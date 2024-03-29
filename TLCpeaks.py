import os
import numpy as np
import pandas as pd
import argparse
import statsmodels.stats.multitest as multi
import tqdm
from multiprocessing import Pool
from scipy.stats import binom_test, combine_pvalues, binom
# dev libraries : 
import timeit

# define function useful for argparse
def boolean_string(s):
    if s not in {'False', 'True'}:
        raise ValueError('Not a valid boolean string')
    return s == 'True'

### Use Argparse to get input variables ### 
### Parse script parameters ###
parser = argparse.ArgumentParser(description='For detailed help consult docs at : https://alexdray86.github.io/TLCpeaks')
parser.add_argument('-i', '--input_bam', type=str,
                    help='Input BAM file to be used to call peaks. [REQUIRED]')
parser.add_argument('-o', '--out_file', type=str,
                    help='Output file with results. [REQUIRED]')
parser.add_argument('-d', '--in_dir', type=str,default=None,
                    help='Input directory containing multiple BAM files. If provided, --input_bam option is ignored. [OPTIONAL]')
parser.add_argument('--out_dir', type=str,default=None,
                    help='Outpud directory if multiple files are provided as input. If provided, --out_file option is ignored. [OPTIONAL]')
parser.add_argument('-t', '--tmp_folder', type=str, default="tmp_folder",
                    help='temporary folder where intermediate results are writen. Typically, one bam file per chromosome and its respective bed file will be written there. This takes a lot of memory space and needs to be removed afterwards. Default is tmp_folder/ [OPTIONAL]')
parser.add_argument('-n', '--n_cpu', type=int, default=4,
                    help='Number of CPU to use for multi-processing. Default is 4. [OPTIONAL]')
parser.add_argument('-w', '--del_weight', type=float, default=0.9,
                    help='Weight to be put on the deletion versus the starting position in the combination of p-values. This weight should be in [0, 1], and the weight assigned to starting position will we (1-w). By default, w=0.9')
parser.add_argument('--min_read_length', type=int, default=17,
        help='Min read length to consider for extracting start positions. Default is 17bp. [OPTIONAL]')
parser.add_argument('--max_read_length', type=int, default=68,
        help='Max read length to consider for extracting start positions. Default is 68bp. [OPTIONAL]')
parser.add_argument('--window_size', type=int, default=200,
        help='Window size to compute statistics on. By default, 200bp are considered, which means +/- 100bp around position of interest. [OPTIONAL]')
parser.add_argument('--adjust_pval', type=boolean_string, default=True, 
        help='Whether adjusted p-value should be used to filter the results. If False, the non-adjusted p-value is used to filter significant peaks. In any cases, the asjtment is done and printed in the output. [OPTIONAL]')
parser.add_argument('-p', '--pval_cutoff', type=float, default=0.05,
                    help='Cut-off on the p-values used to call a peak significant. By default, it will be set at 0.05. If option --adjust_pval is False, then this cutoff will be used on the non-adjusted p-value.')
parser.add_argument('--width_peaks', type=int, default=20,
        help='Once the cross-linking site has been found, this will determine how much we enlarge it to make a peak around it. Default is 20bp peaks. [OPTIONAL]')
parser.add_argument('--method', type=str, default='P_d_and_s',
                    help='Method to compute P-value - define which p-value will be adjusted and used to filter peaks and rank them. *pval_combined* is the old method, all the other ones are the new ones. Options are: pval_combined, P_d_cond_s, P_d_or_s, and P_d_and_s (I will explain next week exactly what they are) Enjoy! Default: P_d_and_s [OPTIONAL]')
parser.add_argument('--mode', type=str, default='ss_and_del',
                   help="['ss_only', 'ss_and_del', 'del_only'] Mode to compute statistical test. 'ss_only' mode uses only starting site to compute the statistics. 'ss_and_del' mode uses both starting site and deletion to compute statistics. 'del_only' mode uses only deletions to compute statistics.")

#'pval_combined','P_d', 'P_s', 'P_d_cond_s', 'P_s_cond_d', 'P_d_or_s', 'P_d_and_s'
args = parser.parse_args()

### Define constants ###
BAM_FILE = args.input_bam
OUT_FILE = args.out_file
MIN_READ_LENGTH = args.min_read_length
MAX_READ_LENGTH = args.max_read_length
WINDOW_SIZE = args.window_size
P_E = 1/(WINDOW_SIZE+1)
WIDTH_PEAK = args.width_peaks
TMP_FOLDER = args.tmp_folder
TMP_CHROM_FILE = TMP_FOLDER + "/list_chromosomes.txt"
N_CPU = args.n_cpu
DEL_WEIGHT = args.del_weight
STR_WEIGHT = 1.0 - DEL_WEIGHT
ADJUST_PVAL = args.adjust_pval
PVAL_CUTOFF = args.pval_cutoff
IN_DIR = args.in_dir
OUT_DIR = args.out_dir
METHOD = args.method
MODE = args.mode
if IN_DIR is None:
    SINGLE_FILE = True
    if BAM_FILE is None or OUT_FILE is None:
        raise Exception('Both input BAM file and output BED file need to be specified. Exiting.')
else:
    SINGLE_FILE = False
    if OUT_DIR is None:
        raise Exception('Input folder with multiple BAM provided. An output folder should be provided too. Exiting.')

# Check inputs
if DEL_WEIGHT < 0.0 or DEL_WEIGHT > 1.0:
    raise Exception('Weight for deletion should be included in [0, 1]. Exiting.')
if (MODE != 'ss_only') and (MODE != 'ss_and_del') and (MODE != 'del_only'):
    raise Exception('Wrong mode provided. Exiting.')

### Class Definition ###
class TLCpeaks(object):
    def __init__(self, chrom, bam_file, tmp_folder, P_D, P_S):
        """Launch peak calling analysis on one chromosome. Outputs a pandas dataframe with all significant peaks

        Args: 
            chrom (string): The chromosome of interest to analyze
            bam_file (string): Path to the bam file containing all chromosomes
            tmp_folder (string): Path to the temp folder where intermediate results are written
            P_D (double): probability of having a deletion among all events (should be precomputed from the bam file)
            P_S (double): probability of having a start site among all events (should be precomputed from the bam file)
        """
        self.chrom      = chrom
        self.bam_file   = bam_file
        self.tmp_folder = tmp_folder
        self.bam_file_1chr = ''
        self.bed_file_1chr = ''
        self.unique_positions = np.empty(0)
        self.unique_positions_pos = np.empty(0)
        self.unique_positions_neg = np.empty(0)
        self.valid_pos_strand = True
        self.valid_neg_strand = True
        self.np_del_pos = np.empty(0) # here
        self.np_del_neg = np.empty(0)
        self.np_str_pos = np.empty(0)
        self.np_str_neg = np.empty(0) 
        self.res_pos_pos, self.res_del_pos, self.res_str_pos = [], [], []
        self.P_d_pos, self.P_s_pos, self.P_d_cond_s_pos, self.P_s_cond_d_pos, self.P_d_or_s_pos, self.P_d_and_s_pos = [], [], [], [], [], []
        self.P_d_neg, self.P_s_neg, self.P_d_cond_s_neg, self.P_s_cond_d_neg, self.P_d_or_s_neg, self.P_d_and_s_neg = [], [], [], [], [], []
        self.res_pos_neg, self.res_del_neg, self.res_str_neg = [], [], []
        self.p_d = P_D
        self.p_s = P_S
        self.pd_res = pd.DataFrame()

    def compute_joint_prob(self, N_d, K_d, N_s, K_s):
        """
        p_d: probability of having a deletion among all deletion+starting site events
        p_s: probability of having a starting site among all deletion+starting site events
        p_e: probability of having an event at position 0 in the window
        N_d: number of deletion in window
        K_d: number of deletion at position 0
        N_s: number of starting site in window
        K_s: number of starting site at position 0
        """
        p_d_e = self.p_d*P_E #joint prob of deletion and event at position 0
        p_s_e = self.p_s*P_E #joint prob of starting site and event at position 0
        N_e = N_d + N_s #total number of event : deletion+starting site

        M_pmf = np.zeros((K_d+1,K_s+1))
        Pd_pmf = np.zeros(K_d+1)
        Ps_pmf = np.zeros(K_s+1)
        # Compute probabilities of having d=0:Kd+1 deletions
        for d in range(K_d+1):
            Pd_pmf[d] = binom.pmf(k=d, n=N_e, p=p_d_e)
        # Compute probabilities of having s=0:Ks+1 starting sites
        for s in range(K_s+1):
            Ps_pmf[s] = binom.pmf(k=s, n=N_e, p=p_s_e)
        # Compute joint probabilities of deletion and starting sites
        for d in range(K_d+1):
            for s in range(K_s+1):
                M_pmf[d,s] = Pd_pmf[d]*Ps_pmf[s]

        # Probability to have more than Kd deletions
        P_d = 1-np.sum(Pd_pmf[0:K_d-1])
        # Probability to have more than Kd deletions
        P_s = 1-np.sum(Ps_pmf[0:K_s-1])
        # Probability to have P(X >= k_x OR Y >= k_y)
        P_or = 1-np.sum(M_pmf[0:K_d-1,0:K_s-1].flatten())
        # Probability to have P(X >= k_x AND Y >= k_y)
        P_and = P_or - (np.sum(Pd_pmf[0:K_d-1]) + np.sum(Ps_pmf[0:K_s-1]) - 2*np.sum(M_pmf[0:K_d-1,0:K_s-1].flatten()))
        # Probability of Deletion conditional on starting sites P(X >= k_x | Y >= k_y)
        if P_s > 0:
            P_d_cond_s = P_and / P_s
        else:
            P_d_cond_s = 0.0
        # Probability of Starting site conditional on Deletions P(Y >= k_y | X >= k_x)
        if P_d > 0:
            P_s_cond_d = P_and / P_d
        else:
            P_s_cond_d = 0.0

        return np.clip([P_d, P_s, P_or, P_and, P_d_cond_s, P_s_cond_d], 0,1)


    def part1_gen_bam_and_bed(self):
    ### PART 1 - Generate bam and bed file from one chromosome ###
        #print('PART 1 - Generate bam and bed file from one chromosome')
        # Generating bam for a single chromosome 
        self.bam_file_1chr = "{0}/{1}.bam".format(self.tmp_folder, self.chrom)
        cmd_ = "samtools view {0} {1} -b > {2}".format(self.bam_file, self.chrom, self.bam_file_1chr)
        os.system(cmd_)
        # Generate a bed file from the bam file, containing the CIGAR information
        # This require bedtools to be installed
        self.bed_file_1chr = self.tmp_folder + "/" + self.chrom + ".bed"
        cmd_str = 'bamtobed -cigar -i {0} > {1}'.format(self.bam_file_1chr, self.bed_file_1chr)
        os.system(cmd_str)
        
    def part2_read_bed_file_record_info(self):
    ### PART 2 - Read bed file and record all information ### 
    # Here we iterate over bed file to record all deletion or start position of each reach reads
        #print('PART 2 - Read bed file and record all information')
        n_selected_reads = 0
        recorded_positions = {}
        list_delet_pos = [] ; list_delet_neg = []
        list_start_pos = [] ; list_start_neg = []
        N_d = 0 ; N_s = 0
        with open(self.bed_file_1chr) as f:
            lines = f.readlines()
            for line in lines:
                fields = line.split("\t")
                chrom = fields[0]
                pos = fields[1]
                end = fields[2]
                strand = fields[5]
                cigar = fields[6]
                read_length = int(end) - int(pos)
                if MODE == 'ss_and_del':
                    if '1D' in cigar : # or '1I' in cigar
                    # option 1 : there is a 1bp deletion event 
                        pos_del = self.find_1d_deletion_position(cigar)
                        N_d+=1
                        if strand == '+':
                            list_delet_pos.append(int(pos) + int(pos_del))
                        else:
                            list_delet_neg.append(int(pos) + int(pos_del))
                    elif '2D' in cigar or '3D' in cigar: # or '3D' in cigar
                    # option 2 : there is a 2bp/3bp deletion event 
                        pos_del = self.find_2d3d_deletion_position(cigar)
                        N_d+=1
                        for pos_d in pos_del:
                            if strand == '+':
                                list_delet_pos.append(int(pos) + int(pos_d))
                            else:
                                list_delet_neg.append(int(pos) + int(pos_d))
                    else:
                        # First filter on read length : 
                        if read_length > 18 and read_length <= 68:
                            N_s+=1
                            if strand == '+':
                                list_start_pos.append(int(pos)-1)
                            elif strand == '-':
                                list_start_neg.append(int(end)+1)
                elif MODE == 'del_only':
                    if '1D' in cigar : # or '1I' in cigar
                    # option 1 : there is a 1bp deletion event 
                        pos_del = self.find_1d_deletion_position(cigar)
                        N_d+=1
                        if strand == '+':
                            list_delet_pos.append(int(pos) + int(pos_del))
                        else:
                            list_delet_neg.append(int(pos) + int(pos_del))
                    elif '2D' in cigar or '3D' in cigar: # or '3D' in cigar
                    # option 2 : there is a 2bp/3bp deletion event 
                        pos_del = self.find_2d3d_deletion_position(cigar)
                        N_d+=1
                        for pos_d in pos_del:
                            if strand == '+':
                                list_delet_pos.append(int(pos) + int(pos_d))
                            else:
                                list_delet_neg.append(int(pos) + int(pos_d))
                elif MODE == 'ss_only':
                    if read_length > 18 and read_length <= 68:
                        N_s+=1
                        if strand == '+':
                            list_start_pos.append(int(pos)-1)
                        elif strand == '-':
                            list_start_neg.append(int(end)+1)

        # Generate unique positions 
        self.unique_positions, c = np.unique(np.concatenate((np.array(list_delet_pos), np.array(list_start_pos))), 
                                                             return_counts=True)
        self.unique_positions_pos = self.unique_positions[c > 1].copy()

        self.unique_positions, c = np.unique(np.concatenate((np.array(list_delet_neg), np.array(list_start_neg))),
                                                             return_counts=True)
        self.unique_positions_neg = self.unique_positions[c > 1].copy()
        
        # Store numpy arrays for delete / start positions 
        # Generating numpy array from lists for better performance 
        self.np_del_pos = np.sort(np.array(list_delet_pos)) ; self.np_del_neg = np.sort(np.array(list_delet_neg))
        self.np_str_pos = np.sort(np.array(list_start_pos)) ; self.np_str_neg = np.sort(np.array(list_start_neg))

        if len(self.unique_positions_pos) < 1:
            self.valid_pos_strand = False
        if len(self.unique_positions_neg) < 1:
            self.valid_neg_strand = False
        return N_d, N_s
            
    def part3_analysis(self):
    ### PART 3 - Analyse all unique positions and generate p-values ###
        #print('PART 3 - Analyse all unique positions and generate p-values')
        prev_position = -1 ; count_line = 0
        del_pos_idx_start, del_pos_idx_end, del_neg_idx_start, del_neg_idx_end = 0, 0, 0, 0
        str_pos_idx_start, str_pos_idx_end, str_neg_idx_start, str_neg_idx_end = 0, 0, 0, 0        
        del_pos_idx_start_summit, del_pos_idx_end_summit, del_neg_idx_start_summit, del_neg_idx_end_summit = 0, 0, 0, 0
        str_pos_idx_start_summit, str_pos_idx_end_summit, str_neg_idx_start_summit, str_neg_idx_end_summit = 0, 0, 0, 0
        
        # Analysis of positive strand
        if self.valid_pos_strand:
            for position in self.unique_positions_pos:
                position = int(position)
                del_pos_idx_start_summit, del_pos_idx_end_summit = self.find_position(del_pos_idx_start_summit, del_pos_idx_end_summit, position, self.np_del_pos)
                str_pos_idx_start_summit, str_pos_idx_end_summit = self.find_position(str_pos_idx_start_summit, str_pos_idx_end_summit, position, self.np_str_pos)
                n_del_pos = len(self.np_del_pos[del_pos_idx_start_summit:del_pos_idx_end_summit])
                n_str_pos = len(self.np_str_pos[str_pos_idx_start_summit:str_pos_idx_end_summit])

                if n_del_pos > 1 or n_str_pos > 1:
                # Here we look at the position w.r.t the positive strand and make statistics if n_event > 1
                    # update indexes 
                    del_pos_idx_start, del_pos_idx_end = self.update_index(del_pos_idx_start, del_pos_idx_end, position, self.np_del_pos, WINDOW_SIZE)
                    window_del_pos = self.np_del_pos[del_pos_idx_start:del_pos_idx_end]
                    N_d = len(window_del_pos)
                    K_d = n_del_pos
                    pval_del = binom_test(n_del_pos, n=len(window_del_pos), p=P_E, alternative='greater')
                    str_pos_idx_start, str_pos_idx_end = self.update_index(str_pos_idx_start, str_pos_idx_end, position, self.np_str_pos, WINDOW_SIZE)
                    window_str_pos = self.np_str_pos[str_pos_idx_start:str_pos_idx_end]
                    N_s = len(window_str_pos)
                    K_s = n_str_pos
                    pval_str = binom_test(n_str_pos, n=len(window_str_pos), p=P_E, alternative='greater')
                    
                    P_d, P_s, P_or, P_and, P_d_cond_s, P_s_cond_d = self.compute_joint_prob(N_d, K_d, N_s, K_s)

                    self.res_pos_pos.append(position) ; self.res_del_pos.append(pval_del) ; self.res_str_pos.append(pval_str)
                    self.P_d_pos.append(P_d), self.P_s_pos.append(P_s), self.P_d_cond_s_pos.append(P_d_cond_s), self.P_s_cond_d_pos.append(P_s_cond_d)
                    self.P_d_or_s_pos.append(P_or), self.P_d_and_s_pos.append(P_and)

                count_line += 1
 
        # Analysis of negative strand
        if self.valid_neg_strand:
            prev_position = -1 ; count_line = 0
            for position in self.unique_positions_neg:
                position = int(position)
                del_neg_idx_start_summit, del_neg_idx_end_summit = self.find_position(del_neg_idx_start_summit, del_neg_idx_end_summit, position, self.np_del_neg)
                str_neg_idx_start_summit, str_neg_idx_end_summit = self.find_position(str_neg_idx_start_summit, str_neg_idx_end_summit, position, self.np_str_neg)

                n_del_neg = len(self.np_del_neg[del_neg_idx_start_summit:del_neg_idx_end_summit])
                n_str_neg = len(self.np_str_neg[str_neg_idx_start_summit:str_neg_idx_end_summit])

                if n_del_neg > 1 or n_str_neg > 1:
                # Here we look at the position w.r.t the negative strand and make statistics if n_event > 1
                    del_neg_idx_start, del_neg_idx_end = self.update_index(del_neg_idx_start, del_neg_idx_end, position, self.np_del_neg, WINDOW_SIZE)        
                    window_del_neg = self.np_del_neg[del_neg_idx_start:del_neg_idx_end]
                    N_d = len(window_del_neg)
                    K_d = n_del_neg
                    str_neg_idx_start, str_neg_idx_end = self.update_index(str_neg_idx_start, str_neg_idx_end, position, self.np_str_neg, WINDOW_SIZE)        
                    window_str_neg = self.np_str_neg[str_neg_idx_start:str_neg_idx_end]
                    N_s = len(window_str_neg)
                    K_s = n_str_neg
                    pval_del = binom_test(n_del_neg, n=len(window_del_neg), p=P_E, alternative='greater')
                    pval_str = binom_test(n_str_neg, n=len(window_str_neg), p=P_E, alternative='greater')

                    P_d, P_s, P_or, P_and, P_d_cond_s, P_s_cond_d = self.compute_joint_prob(N_d, K_d, N_s, K_s)

                    self.res_pos_neg.append(position) ; self.res_del_neg.append(pval_del) ; self.res_str_neg.append(pval_str)
                    self.P_d_neg.append(P_d), self.P_s_neg.append(P_s), self.P_d_cond_s_neg.append(P_d_cond_s), self.P_s_cond_d_neg.append(P_s_cond_d)
                    self.P_d_or_s_neg.append(P_or), self.P_d_and_s_neg.append(P_and)
                count_line += 1
                        
    def part4_assemble_results(self):
    ### PART 4 - Assemble results and keep significant peaks ###
        #print('PART 4 - Assemble results and keep significant peaks')
        comb_pval_pos = []
        for i in range(len(self.res_del_pos)):
            comb_pval = combine_pvalues(np.array([self.res_del_pos[i], self.res_str_pos[i]]), 
                                             method = 'stouffer', 
                                             weights=np.array([DEL_WEIGHT, STR_WEIGHT]))[1]
            if np.isnan(comb_pval):
                comb_pval = min(self.res_del_pos[i], self.res_str_pos[i])
            comb_pval_pos.append(comb_pval)

        comb_pval_neg = []
        for i in range(len(self.res_del_neg)):
            comb_pval = combine_pvalues(np.array([self.res_del_neg[i], self.res_str_neg[i]]), 
                                             method = 'stouffer', 
                                             weights=np.array([DEL_WEIGHT, STR_WEIGHT]))[1]
            if np.isnan(comb_pval):
                comb_pval = min(self.res_del_neg[i], self.res_str_neg[i])
            comb_pval_neg.append(comb_pval)

        pd_res_pos = pd.DataFrame(np.array([np.array([self.chrom for x in range(len(self.res_pos_pos))]), 
                                            np.array(self.res_pos_pos)-int(WIDTH_PEAK/2), np.array(self.res_pos_pos)+int(WIDTH_PEAK/2)+1, 
                                            np.array(['+' for x in range(len(self.res_pos_pos))]),
                                            np.array(self.res_del_pos), np.array(self.res_str_pos), np.array(comb_pval_pos),
                                            np.array(self.P_d_pos), np.array(self.P_s_pos), np.array(self.P_d_cond_s_pos), 
                                            np.array(self.P_s_cond_d_pos), np.array(self.P_d_or_s_pos), np.array(self.P_d_and_s_pos)]),
                              index = ['chr', 'start', 'end', 'strand','pval_del', 'pval_start', 'pval_combined',
                                       'P_d', 'P_s', 'P_d_cond_s', 'P_s_cond_d', 'P_d_or_s', 'P_d_and_s']).T
        #self.P_d_pos, self.P_s_pos, self.P_d_cond_s_pos, self.P_s_cond_d_pos, self.P_d_or_s_pos, self.P_d_and_s_pos
        pd_res_neg = pd.DataFrame(np.array([np.array([self.chrom for x in range(len(self.res_pos_neg))]), 
                                            np.array(self.res_pos_neg)-int(WIDTH_PEAK/2), np.array(self.res_pos_neg)+int(WIDTH_PEAK/2)+1, 
                                            np.array(['-' for x in range(len(self.res_pos_neg))]),
                                            np.array(self.res_del_neg), np.array(self.res_str_neg), np.array(comb_pval_neg),
                                            np.array(self.P_d_neg), np.array(self.P_s_neg), np.array(self.P_d_cond_s_neg), 
                                            np.array(self.P_s_cond_d_neg), np.array(self.P_d_or_s_neg), np.array(self.P_d_and_s_neg)]),
                              index = ['chr', 'start', 'end', 'strand','pval_del', 'pval_start', 'pval_combined',
                                       'P_d', 'P_s', 'P_d_cond_s', 'P_s_cond_d', 'P_d_or_s', 'P_d_and_s']).T

        if self.valid_pos_strand and self.valid_neg_strand:
            self.pd_res = pd.concat((pd_res_pos, pd_res_neg))
            self.pd_res.reset_index(inplace=True, drop=True)
            self.pd_res.sort_values('start', inplace=True)
            self.pd_res['pval_combined'] = pd.to_numeric(self.pd_res['pval_combined'])
        elif self.valid_pos_strand and not self.valid_neg_strand:
            self.pd_res = pd_res_pos
        elif self.valid_neg_strand and not self.valid_pos_strand:
            self.pd_res = pd_res_neg

        if self.valid_pos_strand or self.valid_neg_strand:
            self.pd_res.reset_index(inplace=True, drop=True)
            self.pd_res.sort_values('start', inplace=True)
            self.pd_res['pval_combined'] = pd.to_numeric(self.pd_res['pval_combined'])

    def find_1d_deletion_position(self, cigar_string):
    # This function will read the CIGAR string, and 
    # find deletion's position w.r.t start position
    # (assuming CIGAR string has a single deletion)
        is_num = False
        num_start = 0
        num_end = 0
        pos_in_cigar = 0 ; count_bp = 0
        for c in cigar_string:
            if str.isnumeric(c):
                if is_num is False:
                    is_num = True
                    num_start = pos_in_cigar
            elif c == '\n':
                break
            else:
                if c == 'D' or c == 'I':
                    count_bp += 1
                    break
                else:
                    num_end = pos_in_cigar
                    count_bp += int(cigar_string[num_start:num_end])
                    is_num = False
            pos_in_cigar += 1
        return(count_bp)

    def find_2d3d_deletion_position(self, cigar_string):
    # This function will read the CIGAR string, and 
    # find deletion's position w.r.t start position
    # (assuming CIGAR string has a single deletion)
        is_num = False
        num_start = 0
        num_end = 0
        pos_in_cigar = 0 ; count_bp = 0
        multi_del = [] 
        for c in cigar_string:
            if str.isnumeric(c):
                if is_num is False:
                    is_num = True
                    num_start = pos_in_cigar
            elif c == '\n':
                break
            else:
                if c == 'D':
                    num_end = pos_in_cigar
                    n_del = int(cigar_string[num_start:num_end])
                    for n in range(n_del):
                        multi_del.append(count_bp + n)
                    count_bp += 1
                    break
                else:
                    num_end = pos_in_cigar
                    count_bp += int(cigar_string[num_start:num_end])
                    is_num = False
            pos_in_cigar += 1
        return(multi_del)

    def find_position(self, idx_start, idx_end, position, vect):
        if len(vect) == 0:
            return idx_start,idx_end
        else:
            new_idx_start, new_idx_end = max(min(idx_start, len(vect)-1),0), max(min(idx_end, len(vect)-1),0)
            while(vect[new_idx_start] != position and vect[new_idx_start] < position) and new_idx_start < len(vect)-1:
                new_idx_start += 1
            while(vect[new_idx_end] != position and vect[new_idx_end] < position) and new_idx_end < len(vect)-1:
                new_idx_end += 1
            if vect[new_idx_end] == position and new_idx_end < len(vect)-1:
                while(vect[new_idx_end] == position and new_idx_end < len(vect)-1):
                    new_idx_end += 1
            if new_idx_start == new_idx_end and new_idx_end < len(vect)-1:
                new_idx_end += 1
            if new_idx_end >= len(vect)-1:
                new_idx_end = len(vect) - 1
            if vect[new_idx_start] != position:
                new_idx_start = idx_start
                new_idx_end = idx_start
            return new_idx_start, new_idx_end

    def update_index(self, idx_start, idx_end, position, vect, WINDOW_SIZE):
        if len(vect) == 0:
            return idx_start,idx_end
        else:
            low_bound = position - int(WINDOW_SIZE/2)
            up_bound  = position + int(WINDOW_SIZE/2)
            new_idx_start, new_idx_end = max(min(idx_start, len(vect)-1),0), max(min(idx_end, len(vect)-1),0)
            while(vect[new_idx_start] < low_bound and new_idx_start < len(vect)-1):
                new_idx_start += 1
            while(vect[new_idx_end] < up_bound and new_idx_end < len(vect)-1):
                new_idx_end += 1
            if new_idx_start == new_idx_end and new_idx_end < len(vect)-1:
                new_idx_end += 1
            if new_idx_end >= len(vect)-1:
                new_idx_end = len(vect) - 1
            return new_idx_start, new_idx_end

def remove_folder_and_content(folder):
    for file_name in os.listdir(folder):
        file = path + file_name
        if os.path.isfile(file):
            os.remove(file) 
    os.rmdir(folder)
    
def create_folder_if_not_exists(folder):
    if not os.path.exists(folder):
        os.makedirs(folder)


def generate_chromosome_list(bam_file, tmp_chrom_file):
    cmd_ = "samtools idxstats {0} | cut -f 1 > {1}".format(bam_file, TMP_CHROM_FILE)
    os.system(cmd_)
    chroms = np.array(pd.read_csv(TMP_CHROM_FILE, header=None)[0])
    chroms = chroms[chroms != '*']
    #chroms = chroms[np.array(['random' not in x and 'chrUn' not in x for x in chroms])]
    return(chroms)

def launch_one_chromosome(input_list):
    this_chromosome = input_list[0]
    this_file = input_list[1]
    P_D = input_list[2]
    P_S = input_list[3]
    ChrPrs_object = TLCpeaks(this_chromosome, this_file, TMP_FOLDER, P_D, P_S)
    ChrPrs_object.part1_gen_bam_and_bed()
    ChrPrs_object.part2_read_bed_file_record_info()
    ChrPrs_object.part3_analysis()
    ChrPrs_object.part4_assemble_results()
    return(ChrPrs_object.pd_res)

def launch_one_bam_file(this_bam, out_file):
    start_all = timeit.default_timer()
        
    bai_file = this_bam + '.bai'
    if not os.path.exists(bai_file):
        raise ValueError('Error: no index file was found for ' + this_bam + ', exiting...')

    # Prepare folders and data
    print('Prepare folders and data')
    create_folder_if_not_exists(TMP_FOLDER)
    chroms = generate_chromosome_list(this_bam, TMP_CHROM_FILE)
    chroms = np.flip(chroms)
    list_input = [[x, this_bam] for x in chroms]
   
    # Launch peak calling with multi-processing per chromosome
    print('Launch Deletion/Starting Site event counting with multi-processing per chromosome')
    list_results = []
    with Pool(N_CPU) as p:
        list_results = list(tqdm.tqdm(p.imap(launch_one_chromosome_counts,
                                             list_input),
                              total = len(list_input),
                              position=0, leave=True))
    print(pd.concat(list_results))
    pd_Ns = pd.concat(list_results)
    N_d = np.sum(pd_Ns['N_d'])
    N_s = np.sum(pd_Ns['N_s'])
    P_D = N_d / (N_d + N_s)
    P_S = N_s / (N_d + N_s)
    print("total number of Deletion to consider : {}".format(N_d))
    print("Probability of having a Deletion among all events considered : {}".format(P_D))
    print("total number of starting sites to consider : {}".format(N_s))
    print("Probability of having a Staring site among all events considered : {}".format(P_S))

    # Launch peak calling with multi-processing per chromosome
    print('Launch peak calling with multi-processing per chromosome')
    list_input = [[x, this_bam, P_D, P_S] for x in chroms]

    list_results = []
    with Pool(N_CPU) as p:
        list_results = list(tqdm.tqdm(p.imap(launch_one_chromosome,
                                             list_input),
                              total = len(list_input),
                              position=0, leave=True))

    # Combine results
    print('Combine results and adjust p-values for multiple-testing')
    pd_general_results = pd.concat(list_results)
    pd_general_results[METHOD] = pd_general_results[METHOD].astype(float)
    pd_general_results['pval_combined_adj'] = multi.multipletests(pd_general_results[METHOD])[1]
    if ADJUST_PVAL:
        pd_general_results_sub = pd_general_results.iloc[np.array(pd_general_results['pval_combined_adj'] < PVAL_CUTOFF)].copy()
        pd_general_results_sub.sort_values('pval_combined_adj', inplace=True)
    else:
        pd_general_results_sub = pd_general_results.iloc[np.array(pd_general_results[METHOD] < PVAL_CUTOFF)].copy()
        pd_general_results_sub.sort_values(METHOD, inplace=True)

    pd_general_results_sub = pd_general_results_sub[['chr', 'start', 'end', 'pval_del', 'pval_start','strand', METHOD, 'pval_combined_adj']]
    pd_general_results_sub.to_csv(out_file, sep="\t", index=False, header=False)
    
    stop_all = timeit.default_timer()
    print('Running time of TLCpeaks for {0} : {1} min'.format(this_bam, str( round( (stop_all - start_all)/60, 3)))) 


def launch_one_chromosome_counts(input_list):
    this_chromosome = input_list[0]
    this_file = input_list[1]
    ChrPrs_object = TLCpeaks(this_chromosome, this_file, TMP_FOLDER, 0, 0)
    ChrPrs_object.part1_gen_bam_and_bed()
    N_d, N_s = ChrPrs_object.part2_read_bed_file_record_info()
    return pd.DataFrame(np.array([[N_d, N_s]]), columns=['N_d','N_s'])

### MAIN ###
if __name__ == "__main__":
    
    if SINGLE_FILE:
        start_all = timeit.default_timer()
        
        bai_file = BAM_FILE + '.bai'
        if not os.path.exists(bai_file):
            raise ValueError('Error: no index file was found for ' + BAM_FILE + ', exiting...')

        # Prepare folders and data
        print('Prepare folders and data')
        create_folder_if_not_exists(TMP_FOLDER)
        chroms = generate_chromosome_list(BAM_FILE, TMP_CHROM_FILE)
        chroms = np.flip(chroms)
        list_input = [[x, BAM_FILE] for x in chroms]
       
        # Launch peak calling with multi-processing per chromosome
        print('Launch Deletion/Starting Site event counting with multi-processing per chromosome')
        list_results = []
        with Pool(N_CPU) as p:
            list_results = list(tqdm.tqdm(p.imap(launch_one_chromosome_counts,
                                                 list_input),
                                  total = len(list_input),
                                  position=0, leave=True))
        print(pd.concat(list_results))
        pd_Ns = pd.concat(list_results)
        N_d = np.sum(pd_Ns['N_d'])
        N_s = np.sum(pd_Ns['N_s'])
        P_D = N_d / (N_d + N_s)
        P_S = N_s / (N_d + N_s)
        print("total number of Deletion to consider : {}".format(N_d))
        print("Probability of having a Deletion among all events considered : {}".format(P_D))
        print("total number of starting sites to consider : {}".format(N_s))
        print("Probability of having a Staring site among all events considered : {}".format(P_S))

        # Launch peak calling with multi-processing per chromosome
        print('Launch peak calling with multi-processing per chromosome')
        list_input = [[x, BAM_FILE, P_D, P_S] for x in chroms]
 
        list_results = []
        with Pool(N_CPU) as p:
            list_results = list(tqdm.tqdm(p.imap(launch_one_chromosome,
                                                 list_input),
                                  total = len(list_input),
                                  position=0, leave=True))

        # Combine results
        print('Combine results and adjust p-values for multiple-testing')
        pd_general_results = pd.concat(list_results)
        pd_general_results[METHOD] = pd_general_results[METHOD].astype(float)
        pd_general_results['pval_combined_adj'] = multi.multipletests(pd_general_results[METHOD])[1]
        if ADJUST_PVAL:
            pd_general_results_sub = pd_general_results.iloc[np.array(pd_general_results['pval_combined_adj'] < PVAL_CUTOFF)].copy()
            pd_general_results_sub.sort_values('pval_combined_adj', inplace=True)
        else:
            pd_general_results_sub = pd_general_results.iloc[np.array(pd_general_results[METHOD] < PVAL_CUTOFF)].copy()
            pd_general_results_sub.sort_values(METHOD, inplace=True)

        pd_general_results_sub = pd_general_results_sub[['chr', 'start', 'end', 'pval_del', 'pval_start','strand', METHOD, 'pval_combined_adj']]
        pd_general_results_sub.to_csv(OUT_FILE, sep="\t", index=False, header=False)
        
        stop_all = timeit.default_timer()
        print('Running time of TLCpeaks for {0} : {1} min'.format(BAM_FILE, str( round( (stop_all - start_all)/60, 3)))) 
    else:
        create_folder_if_not_exists(OUT_DIR)
        list_files = os.listdir(IN_DIR)
        list_files = [x for x in list_files if '.bam' in x and '.bai' not in x]
        for this_bam in list_files:
            bai_file = IN_DIR + "/" + this_bam + '.bai'
            bam_path = IN_DIR + "/" + this_bam
            if os.path.exists(bai_file):
                print('working with ' + this_bam)
                out_file = OUT_DIR + "/" + this_bam.rsplit('.', 1)[0] + '.bed'
                launch_one_bam_file(bam_path, out_file)
            else:
                print('Skipping ' + bam_path + ', no index .bai found !')


