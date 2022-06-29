### INSTALL TLCpeaks ### 

# First Git Clone the repository on your computer with : 
git clone https://github.com/alexdray86/TLCpeaks.git

# Go inside the directory, and create the python3 environment :
python3 -m venv env_TLCpeaks

# Now, source the environment :
source env_TLCpeaks/bin/activate

# Then, relevant libraries are installed in the environment using the provided file requirement.txt :
python3 -m pip install -r requirements.txt


### LAUNCH AND USE TLCpeaks ### 

# In its most basic form, TLCpeaks can be used with 
python3 TLCpeaks.py -i INPUT_BAM -o OUTPUT_BED

# The script therefore use a BAM file as input and output a bed file
containing significant peak called

# To get help and use the options, just launch : 
python3 TLCpeaks.py --help

# For detailed decumentation, go to 
XXX to add XXX

