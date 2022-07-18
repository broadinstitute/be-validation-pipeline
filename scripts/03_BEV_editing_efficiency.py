import pandas as pd
import matplotlib as mpl
import seaborn as sns

import core_functions as cfs


mpl.rc('pdf', fonttype=42)
mpl.rcParams['font.sans-serif'] = "Arial"
mpl.rcParams['font.family'] = "sans-serif"
sns.set_context('paper')
sns.set_style('ticks')
boxprops = {'edgecolor': 'k', 'linewidth': 0.5, 'facecolor': 'w'}
lineprops = {'color': 'k', 'linewidth': 0.5}
stripplot_kwargs = dict({'linewidth': 0.5, 'size': 3, 'alpha': 0.8})
boxplot_kwargs = dict({'boxprops': {'linewidth':0,'facecolor':'w'}, 'medianprops': {'linewidth':1},
                       'whiskerprops': {'linewidth':0}, 'capprops': {'linewidth':0},
                       'width': 0.8,'whis':(10,90)})

# Read meta-information file
global old_input_file
# input_filepath = 'AudreyData/Metainfo_input_f2_r1_ABE.csv'
input_filepath = input("Please enter input filepath here: ")
old_input_file = pd.read_csv(input_filepath)

# User inputs filepath where nucleotide_percentage outputs are stored
global nuc_per_filepath
# nuc_per_filepath = 'AudreyData/Validation_CRISPResso_results/nucleotide_percentage/'
nuc_per_filepath = input("Please enter filepath to folder with nucleotide percentage outputs here: ")
nuc_per_filepath = cfs.check_folder_filepath(nuc_per_filepath)
print(nuc_per_filepath)

# User inputs filepath to store editing efficency outputs 
global output_filepath
# output_filepath = 'AudreyData/Validation_CRISPResso_results/'
output_filepath = input("Please enter output folder filepath here: ")
output_filepath = cfs.check_folder_filepath(output_filepath)
print(output_filepath)

# User inputs name of gene being targeted
global experiment
experiment = input("Please enter gene name here: ")

# User inputs type of base editor
global be_type

be_type_input = input("Please specify the type of base editor used in the input samples by entering 'A' for A base editor or 'C' for C base editor: ")

# Make sure a base editor is selected and not default value

if (be_type_input != 'A') and (be_type_input != 'C'):
    raise Exception('Invalid input. Please enter either A or C to specify base editor.')

else:
    be_type = be_type_input + 'BE'
    
# Generate editing efficiency table
df = cfs.run_editing_eff(old_input_file, nuc_per_filepath, output_filepath, experiment, be_type)
save_filepath = output_filepath + 'editing_efficiency_input.csv'
df.to_csv(save_filepath)

# Filter by editing window for plotting purposes
df = df.copy().drop('sg', axis=1).drop_duplicates()
if be_type == 'ABE':
    df.loc[(df['Position'] >= 4) & (df['Position'] <= 8) & (df['Nucleotide'] == 'G'),'Position'].value_counts()
elif be_type == 'CBE':
    df.loc[(df['Position'] >= 4) & (df['Position'] <= 8) & (df['Nucleotide'] == 'T'),'Position'].value_counts()
    
nuc = input("Please enter target nucleotide here: ")

fig,ax = cfs.make_editing_eff_plot(df, nuc,[experiment], be_type, output_filepath)
