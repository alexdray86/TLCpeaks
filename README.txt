### INSTALL TLCpeaks ### 

# First, the python3 environment should be set up :
python3 -m venv env_TLCpeaks

# The environment can be sourced to be used :
source env_TLCpeaks/bin/activate

# Then, relevant libraries are install in the environment using requirement.txt :
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

