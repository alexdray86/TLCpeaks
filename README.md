_________________________________________________________________________________

## TLCpeaks

We present a novel peak calling algorithm called TLCpeaks, which is specifically designed for identifying RNA-binding sites in transcriptome-wide crosslinking and immunoprecipitation (CLIP) data. CLIP is a technique that utilizes covalent crosslinking of RNA-binding proteins (RBPs) to their target RNA to identify binding sites. Current peak callers for CLIP data are limited in their reproducibility and accuracy in predicting high-confidence peaks, and involve high computational burden. TLCpeaks takes advantage of short nucleotide deletions, which are a common feature of current CLIP protocols at various frequencies, together with starting sites of truncated reads. We show that peaks called by TLCpeaks have a higher motif content, better robustness between replicates, and are closer to motif sequences than peaks produced by standard methods PureCLIP and CLIPper. Additionally, our algorithm is highly efficient, taking only about 9 minutes per sample, which is about 20 times faster than PureCLIP and 50 times faster than CLIPper. Our results suggest that TLCpeaks can infer RNA-binding sites with unprecedented speed and accuracy by working in synergy with the TLC-CLIP experimental method. Our algorithm has the potential to enhance the reliability and efficiency of identifying RNA-binding sites, which can greatly contribute to the understanding of gene regulation and cellular processes.

_________________________________________________________________________________

### INSTALL TLCpeaks ###

#### First Git Clone the repository on your computer with : 
git clone https://github.com/alexdray86/TLCpeaks.git


_________________________________________________________________________________

### Create and launch Python Environment

#### Go inside the directory, and create the python3 environment :
python3 -m venv env_TLCpeaks

#### Now, source the environment :
source env_TLCpeaks/bin/activate

#### Then, relevant libraries are installed in the environment using the provided file requirement.txt :
python3 -m pip install -r requirements.txt

_________________________________________________________________________________

### LAUNCH AND USE TLCpeaks ## 

#### In its most basic form, TLCpeaks can be used with 
python3 TLCpeaks.py -i INPUT_BAM -o OUTPUT_BED

The script therefore use a BAM file as input and output a bed file
containing significant peak called

#### To get help and use the options, just launch : 
python3 TLCpeaks.py --help

