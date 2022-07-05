import os
import numpy as np
import pandas as pd
import argparse
import statsmodels.stats.multitest as multi
import tqdm
from multiprocessing import Pool
from scipy.stats import binom_test, combine_pvalues
# dev libraries : 
import timeit

### Use Argparse to get input variables ### 
### Parse script parameters ###
parser = argparse.ArgumentParser(description='For detailed help consult docs at : https://alexdray86.github.io/TLCpeaks')
parser.add_argument('-i', '--input_bam', type=str,
                    help='Input BAM file to be used to call peaks. [REQUIRED]')
parser.add_argument('-o', '--out_file', type=str,
                    help='Output file with results. [REQUIRED]')
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
def boolean_string(s):
    if s not in {'False', 'True'}:
        raise ValueError('Not a valid boolean string')
    return s == 'True'
parser.add_argument('--adjust_pval', type=boolean_string, default=True, 
        help='Whether adjusted p-value should be used to filter the results. If False, the non-adjusted p-value is used to filter significant peaks. In any cases, the asjtment is done and printed in the output. [OPTIONAL]')


args = parser.parse_args()

### Define constants ###
BAM_FILE = args.input_bam
OUT_FILE = args.out_file
MIN_READ_LENGTH = args.min_read_length
MAX_READ_LENGTH = args.max_read_length
WINDOW_SIZE = args.window_size
TMP_FOLDER = args.tmp_folder
TMP_CHROM_FILE = TMP_FOLDER + "/list_chromosomes.txt"
N_CPU = args.n_cpu
DEL_WEIGHT = args.del_weight
STR_WEIGHT = 1.0 - DEL_WEIGHT
ADJUST_PVAL = args.adjust_pval
print(ADJUST_PVAL)

# Check inputs
if DEL_WEIGHT < 0.0 or DEL_WEIGHT > 1.0:
    raise Exception('Weight for deletion should be included in [0, 1]. Exiting.')
if BAM_FILE is None or OUT_FILE is None:
    raise Exception('Both input BAM file and output BED file need to be specified. Exiting.')

### Class Definition ###

class ChromosomeParser(object):
    def __init__(self, chrom, bam_file, tmp_folder):
        """Launch peak calling analysis on one chromosome. Outputs a pandas dataframe with all significant peaks

        Args: 
            chrom (string): The chromosome of interest to analyze
            bam_file (string): Path to the bam file containing all chromosomes
            tmp_folder (string): Path to the temp folder where intermediate results are written
        """
        self.chrom      = chrom
        self.bam_file   = bam_file
        self.tmp_folder = tmp_folder
        self.bam_file_1chr = ''
        self.bed_file_1chr = ''
        self.unique_positions = np.empty(0)
        self.unique_positions_gt1 = np.empty(0)
        self.valid_chromosome = True
        self.np_del_pos = np.empty(0) # here
        self.np_del_neg = np.empty(0)
        self.np_str_pos = np.empty(0)
        self.np_str_neg = np.empty(0) 
        self.res_pos_pos, self.res_del_pos, self.res_str_pos = [], [], []
        self.res_pos_neg, self.res_del_neg, self.res_str_neg = [], [], []
        self.pd_res = pd.DataFrame()
        
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
                if '1D' in cigar : # or '1I' in cigar
                # option 1 : there is a 1bp deletion event 
                    pos_del = self.find_1d_deletion_position(cigar)
                    if strand == '+':
                        list_delet_pos.append(int(pos) + int(pos_del))
                    else:
                        list_delet_neg.append(int(pos) + int(pos_del))
                elif '2D' in cigar or '3D' in cigar: # or '3D' in cigar
                # option 2 : there is a 2bp/3bp deletion event 
                    pos_del = self.find_2d3d_deletion_position(cigar)
                    for pos_d in pos_del:
                        if strand == '+':
                            list_delet_pos.append(int(pos) + int(pos_d))
                        else:
                            list_delet_neg.append(int(pos) + int(pos_d))
                else:
                    # First filter on read length : 
                    if read_length > 18 and read_length <= 68:
                        if strand == '+':
                            list_start_pos.append(int(pos))
                        elif strand == '-':
                            list_start_neg.append(int(end))
        # Generate unique positions 
        self.unique_positions, c = np.unique(np.concatenate((np.unique(list_delet_pos), np.unique(list_delet_neg),
                                                          np.unique(list_start_pos), np.unique(list_start_neg))),
                                                          return_counts=True)

        self.unique_positions_gt1 = self.unique_positions[c > 1].copy()

        # Store numpy arrays for delete / start positions 
        # Generating numpy array from lists for better performance 
        self.np_del_pos = np.sort(np.array(list_delet_pos)) ; self.np_del_neg = np.sort(np.array(list_delet_neg))
        self.np_str_pos = np.sort(np.array(list_start_pos)) ; self.np_str_neg = np.sort(np.array(list_start_neg))

        if len(self.unique_positions_gt1) < 1:
            self.valid_chromosome = False
            
            
    def part3_analysis(self):
    ### PART 3 - Analyse all unique positions and generate p-values ###
        #print('PART 3 - Analyse all unique positions and generate p-values')
        proba = 1/(WINDOW_SIZE+1) # probability to get 1-bp by chance out of 201 bp window
        prev_position = -1 ; count_line = 0
        del_pos_idx_start, del_pos_idx_end, del_neg_idx_start, del_neg_idx_end = 0, 0, 0, 0
        str_pos_idx_start, str_pos_idx_end, str_neg_idx_start, str_neg_idx_end = 0, 0, 0, 0        
        del_pos_idx_start_summit, del_pos_idx_end_summit, del_neg_idx_start_summit, del_neg_idx_end_summit = 0, 0, 0, 0
        str_pos_idx_start_summit, str_pos_idx_end_summit, str_neg_idx_start_summit, str_neg_idx_end_summit = 0, 0, 0, 0
        for position in self.unique_positions_gt1:
            perc_adv = round(100 * count_line / len(self.unique_positions_gt1))
            position = int(position)
            del_pos_idx_start_summit, del_pos_idx_end_summit = self.find_position(del_pos_idx_start_summit, del_pos_idx_end_summit, position, self.np_del_pos)
            del_neg_idx_start_summit, del_neg_idx_end_summit = self.find_position(del_neg_idx_start_summit, del_neg_idx_end_summit, position, self.np_del_neg)
            str_pos_idx_start_summit, str_pos_idx_end_summit = self.find_position(str_pos_idx_start_summit, str_pos_idx_end_summit, position, self.np_str_pos)
            str_neg_idx_start_summit, str_neg_idx_end_summit = self.find_position(str_neg_idx_start_summit, str_neg_idx_end_summit, position, self.np_str_neg)

            n_del_pos = len(self.np_del_pos[del_pos_idx_start_summit:del_pos_idx_end_summit])
            n_del_neg = len(self.np_del_neg[del_neg_idx_start_summit:del_neg_idx_end_summit])
            n_str_pos = len(self.np_str_pos[str_pos_idx_start_summit:str_pos_idx_end_summit])
            n_str_neg = len(self.np_str_neg[str_neg_idx_start_summit:str_neg_idx_end_summit])

            if n_del_pos > 1 or n_str_pos > 1:
            # Here we look at the position w.r.t the positive strand and make statistics if n_event > 1
                # update indexes 
                del_pos_idx_start, del_pos_idx_end = self.update_index(del_pos_idx_start, del_pos_idx_end, position, self.np_del_pos, WINDOW_SIZE)
                window_del_pos = self.np_del_pos[del_pos_idx_start:del_pos_idx_end]
                pval_del = binom_test(n_del_pos, n=len(window_del_pos), p=proba, alternative='greater')
                str_pos_idx_start, str_pos_idx_end = self.update_index(str_pos_idx_start, str_pos_idx_end, position, self.np_str_pos, WINDOW_SIZE)
                window_str_pos = self.np_str_pos[str_pos_idx_start:str_pos_idx_end]
                pval_str = binom_test(n_str_pos, n=len(window_str_pos), p=proba, alternative='greater')
                self.res_pos_pos.append(position) ; self.res_del_pos.append(pval_del) ; self.res_str_pos.append(pval_str)

            if n_del_neg > 1 or n_str_neg > 1:
            # Here we look at the position w.r.t the negative strand and make statistics if n_event > 1
                del_neg_idx_start, del_neg_idx_end = self.update_index(del_neg_idx_start, del_neg_idx_end, position, self.np_del_neg, WINDOW_SIZE)        
                window_del_neg = self.np_del_neg[del_neg_idx_start:del_neg_idx_end]
                str_neg_idx_start, str_neg_idx_end = self.update_index(str_neg_idx_start, str_neg_idx_end, position, self.np_str_neg, WINDOW_SIZE)        
                window_str_neg = self.np_str_neg[str_neg_idx_start:str_neg_idx_end]
                pval_del = binom_test(n_del_neg, n=len(window_del_neg), p=proba, alternative='greater')
                pval_str = binom_test(n_str_neg, n=len(window_str_neg), p=proba, alternative='greater')
                self.res_pos_neg.append(position) ; self.res_del_neg.append(pval_del) ; self.res_str_neg.append(pval_str)
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
                                            np.array(self.res_pos_pos)-10, np.array(self.res_pos_pos)+10, 
                                            np.array(['+' for x in range(len(self.res_pos_pos))]),
                                            np.array(self.res_del_pos), np.array(self.res_str_pos), np.array(comb_pval_pos)]),
                              index = ['chr', 'start', 'end', 'strand','pval_del', 'pval_start', 'pval_combined']).T

        pd_res_neg = pd.DataFrame(np.array([np.array([self.chrom for x in range(len(self.res_pos_neg))]), 
                                            np.array(self.res_pos_neg)-10, np.array(self.res_pos_neg)+10, 
                                            np.array(['-' for x in range(len(self.res_pos_neg))]),
                                            np.array(self.res_del_neg), np.array(self.res_str_neg), np.array(comb_pval_neg)]),
                              index = ['chr', 'start', 'end', 'strand', 'pval_del', 'pval_start', 'pval_combined']).T

        self.pd_res = pd.concat((pd_res_pos, pd_res_neg))

        if self.pd_res.shape[0] >=2:     
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
    cmd_ = "samtools idxstats {0} | cut -f 1 > {1}".format(BAM_FILE, TMP_CHROM_FILE)
    os.system(cmd_)
    chroms = np.array(pd.read_csv(TMP_CHROM_FILE, header=None)[0])
    chroms = chroms[chroms != '*']
    chroms = chroms[np.array(['random' not in x and 'chrUn' not in x for x in chroms])]
    return(chroms)

def launch_one_chromosome(this_chromosome):
    ChrPrs_object = ChromosomeParser(this_chromosome, BAM_FILE, TMP_FOLDER)
    ChrPrs_object.part1_gen_bam_and_bed()
    ChrPrs_object.part2_read_bed_file_record_info()
    if ChrPrs_object.valid_chromosome:
        ChrPrs_object.part3_analysis()
        ChrPrs_object.part4_assemble_results()
        return(ChrPrs_object.pd_res)


### MAIN ###
if __name__ == "__main__":
    start_all = timeit.default_timer()

    # Prepare folders and data
    print('Prepare folders and data')
    create_folder_if_not_exists(TMP_FOLDER)
    chroms = generate_chromosome_list(BAM_FILE, TMP_CHROM_FILE)
    chrome = np.flip(chroms)

    # Launch peak calling with multi-processing per chromosome
    print('Launch peak calling with multi-processing per chromosome')
    list_results = []
    with Pool(N_CPU) as p:
        list_results = list(tqdm.tqdm(p.imap(launch_one_chromosome,
                                             chroms),
                              total = len(chroms),
                              position=0, leave=True))

    # Combine results
    print('Combine results and adjust p-values for multiple-testing')
    pd_general_results = pd.concat(list_results)
    pd_general_results['pval_combined_adj'] = multi.multipletests(pd_general_results['pval_combined'])[1]
    if ADJUST_PVAL:
        pd_general_results_sub = pd_general_results.iloc[np.array(pd_general_results['pval_combined_adj'] < 0.05)].copy()
        pd_general_results_sub.sort_values('pval_combined_adj', inplace=True)
    else:
        pd_general_results_sub = pd_general_results.iloc[np.array(pd_general_results['pval_combined'] < 0.05)].copy()
        pd_general_results_sub.sort_values('pval_combined', inplace=True)

    pd_general_results_sub = pd_general_results_sub[['chr', 'start', 'end', 'pval_del', 'pval_start','strand', 'pval_combined', 'pval_combined_adj']]
    pd_general_results_sub.to_csv(OUT_FILE, sep="\t", index=False, header=False)
    
    stop_all = timeit.default_timer()
    print('Running time of TLCpeaks for {0} : {1} min'.format(BAM_FILE, str( (stop_all - start_all)/60))) 
