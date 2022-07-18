# Import functions and modules

import pandas as pd
import matplotlib as mpl

import core_functions as cfs

mpl.rc('pdf', fonttype=42)
mpl.rcParams['font.sans-serif'] = "Arial"
mpl.rcParams['font.family'] = "sans-serif"

# Read meta-information file
input_filepath = input("Please enter input filepath here: ")
input_file = pd.read_csv(input_filepath)
print(input_file.head(2))

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

# User inputs type of base editor to determine x-axis labels (e.g. ABE -> Position of A)
global be_type
be_type_input = input("Please specify the type of base editor used in the input samples by entering 'A' for A base editor or 'C' for C base editor: ")

# Make sure a base editor is selected and not default value
if (be_type_input != 'A') and (be_type_input != 'C'):
    raise Exception('Invalid input. Please enter either A or C to specify base editor.')

else:
    # Add 'BE' to be_type
    be_type = be_type_input + 'BE'

#make input df with columns BEV, offset, rev_com, left_lim, right_lim, primer, width, height
input_df = pd.DataFrame()

#check for NaN values i.e. blank rows
if input_file.isnull().values.any(): 
    input_file = cfs.clean_input_file(input_file)

allele_freq_input_df = input_file


BEV_test_df = allele_freq_input_df[['sgRNA_sequence', 'BEV_test', 'primer']].copy()
BEV_test_df['type'] = 'test'
BEV_test_df = BEV_test_df.copy().rename(columns = {'sgRNA_sequence': 'guide_seq', 'BEV_test':'BEV'})

BEV_ref_df = allele_freq_input_df[['sgRNA_sequence', 'BEV_ref', 'primer']].copy()
BEV_ref_df['type'] = 'ref'
BEV_ref_df = BEV_ref_df.copy().rename(columns = {'sgRNA_sequence': 'guide_seq', 'BEV_ref':'BEV'})

input_df = pd.concat([BEV_ref_df, BEV_test_df]).reset_index(drop=True)

# Set x-axis limits
input_df['left_lim'] = -25
input_df['right_lim'] = 25

# Set plot dimensions
input_df['width'] = 6 
input_df['height'] = 6 

print(input_df.head(2))


# Run functions row by row
for i,row in input_df.iterrows():
    print(row['BEV'])
    bev_list = row['BEV'].split(';')
    output_name = '_'.join(bev_list)#+'_'+row['primer']
    #bev_df = get_bev_df(bev_list,row['rev_com'],output_name,row['primer'], row['guide_seq'])
    bev_df = cfs.get_bev_df(bev_list,output_name,row['guide_seq'], be_type, CRISPResso_filepath, 
                          bev_string_id, output_filepath)
    cfs.make_nuc_per_plot(bev_df,row['left_lim'],row['right_lim'],bev_list, row['width'],row['height'],
                    be_type, output_name, output_filepath)
    
    #break
