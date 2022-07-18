# Import packages and modules
import sys
# sys.path.append('scripts/')
import pandas as pd 
import core_functions as cfs

modules = ['pandas', 'numpy']
for module in modules:
    try:
        print(module + ' ' + sys.modules[module].__version__)
    except:
        print(module + ' has no __version__ attribute')

# Read meta-information file
input_filepath = input("Please enter input filepath here: ")
input_file = pd.read_csv(input_filepath)
print(input_file.head(2))

# check for NaN values i.e. blank rows
if input_file.isnull().values.any(): 
    input_file = cfs.clean_input_file(input_file)
    
# check if inputs in meta-information file are correctly formatted
cfs.check_input_file(input_file)

# User inputs read_count_filter
global read_count_percent_cutoff
read_count_percent_cutoff_input = input("Please enter percent to use for filtering low read counts here. "
                                  "Common cutoffs are 1-2%. ")
read_count_percent_cutoff = float(read_count_percent_cutoff_input)

# User inputs path to correlation input file, or leaves blank if none
corr_corr_inputpath = input("Please enter correlation input filepath here: ")
if corr_corr_inputpath == '':
    corr_input = None
    print('No correlation input file.')
else:
    corr_input = pd.read_csv(corr_corr_inputpath)
    corr_input.head()

# User inputs string used to name files ('BEV' or 'NGBEV')
global bev_string_id
bev_string_id = input('Please enter either \'BEV\' or \'NGBEV\' to indicate which string is used when naming your CRISPResso files.')
if ((bev_string_id != 'BEV') and (bev_string_id != 'NGBEV')):
    raise Exception('Invalid input. Please enter either \'BEV\' or \'NGBEV\' to specify which string is used in CRISPResso file names. Be careful not to add any extra spaces.')
    
    
# User inputs path to CRISPResso files
global CRISPResso_filepath 
CRISPResso_filepath = input("Please enter CRISPResso filepath here: ")
CRISPResso_filepath = cfs.check_folder_filepath(CRISPResso_filepath)
print(CRISPResso_filepath)

# User inputs filepath to store allele_freq output tables 
global output_filepath 
# output_filepath = '../../RemovingPrimerInput_NotebookUpdate/allele_freq'
output_filepath = input("Please enter output folder filepath here: ")
output_filepath = cfs.check_folder_filepath(output_filepath)
print(output_filepath)

# Run function
cfs.run(input_file,corr_input, bev_string_id, CRISPResso_filepath, read_count_percent_cutoff, output_filepath)

## HEATMAP CODE ##

# User inputs column labels and heatmap filename
global time1
time1 = input("Please input the label for the first %Reads column (e.g. D7, low)")
global time2
time2 = input("Please input the label for the second %Reads column (e.g. D21, high)")
global heatmap_filename
heatmap_filename = input("Please enter the filename for the heatmap.")

# Run heatmap function
cfs.heatmaps(df = input_file,
         vmin = -2,
         vmax = 2,
         filepath = output_filepath,
         time1 = time1,
         time2 = time2,
         filename = heatmap_filename
        )


