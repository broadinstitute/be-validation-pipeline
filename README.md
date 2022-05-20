# Base Editor Validation (BEV) Pipeline

Welcome to the BEV pipeline! This analysis pipeline begins after you have received sequencing files back after conducting a base editor validation experiment. 

### Before beginning, please download all the files in this repo by clicking on the green "Code" button, then selecting "Download ZIP." This zip file will then be saved to your Downloads folder, which you can then move to wherever you'd like to conduct your analysis. 

### Please read and refer to the Sequence_Orientation_Documentation.ipynb in the "docs" folder to ensure that all sequence inputs are in the correct orientation.

## Part 1: Run BEV on GPP LIMS

[CRISPResso2](https://github.com/pinellolab/CRISPResso2) is our preferred tool for analysis of validation outcomes. To begin, please follow the steps described [here](https://gpp-rnd.github.io/be-validation-pipeline/). 

## Part 2: Run BEV notebooks

There are three Jupyter notebooks that use files outputted by CRISPResso to generate useful visualizations. They must be run in the following order:
1. Allele frequencies notebook: 
    * Generates translations for resulting alleles 
    * Visualizes their log-fold changes (relative to an early time point) and their frequencies as heat maps
2. Nucleotide percentage notebook:
    * Calculates percentage frequency of each nucleotide at each target nucleotide position
    * Visualizes percentage frequencies as bar charts broken down by nucleotide at each position 
3. Editing efficiency notebook:
    * Uses outputs from nucleotide percentage notebook to visualize overall editing efficiency (by nucleotide position)

### Step 1: Create Metainformation File

This file will be used in each of the three BEV notebooks. **Please refer to the Sequence_Orientation_Documentation.ipynb 
in the docs folder to make sure sequence inputs are formatted correctly.** 
<br/><br/>
**Columns**: 

* **sg** : a short unique identifier for each row
* **sgRNA_sequence** : sequence of sgRNA as designed (please see 
  [Sequence_Orientation_Documentation.ipynb](docs/Sequence_Orientation_Documentation.ipynb) 
  for more information)
* **translation_ref_seq**: this sequence can be retrieved by following the steps outlined in 
  [Sequence_Orientation_Documentation.ipynb](docs/Sequence_Orientation_Documentation.ipynb) 
    * reference sequence outputted by CRISPResso formatted such that any intronic sequences are lower-case, 
      exons are upper-case, and UTRs are indicated by square brackets (if applicable) <u> must be sequence on strand 
      that is being translated; may not necessarily be the same strand as the sgRNA sequence</u> 
    * Ex. <font color='grey'>tgtcttttctatgatctctttag</font><font color='green'>GGGTGACCCAGTCTATT</font>
* **primer** : name of primer pair (joined by '\_') used to amplify genomic locus as mentioned in sample name; 
  if the two samples being compared in that row have different primers, leave this column blank
    * Ex 1. <font color='purple'>F_C12</font><font color = 'blue'><b>_</b></font><font color='green'>R_C12</font>
    * Ex 2. F3A_R2B
* **frame** : frame for translation (manually determined for each sg / primer pair); position of first coding nucleotide in reference sequence within codon; frame can be 1, 2, 3
    * Ex. given reference sequence: tgtcttttctatgatctctttag<font color='green'>**G**</font>G|GTG|ACC|CAG|TCT|ATT 
        since the first coding nucleotide of the reference sequence (<font color='green'><b>G</b></font>) is the 2nd nucleotide in its codon 
        (\_<font color='green'><b>G</b></font>G) &rightarrow; frame = 2
* **first_codon** : first codon for translation (note: this codon may be incomplete in the translation_ref_seq input)
    * Ex. given reference sequence: tgtcttttctatgatctctttag**GG**|GTG|ACC|CAG|TCT|ATT 
        the first codon input would be **T**GG where the **T** was in the previous exon
* **last_codon** : last codon for translation (note: this codon may be incomplete in the translation_ref_seq input)
    * Ex. given reference sequence: tgtcttttctatgatctctttagGG|GTG|ACC|CAG|TCT|**ATT** 
        the last codon input would be ATT 
* **rev_com** : "True" if the guide sequence/amplicon sequence are on the opposite strand from the strand being translated
* **BEV_ref** : reference sample(s) for log-fold change (LFC) calculation (i.e. early time point, empty vector, etc.); if multiple BEV numbers are given, they should be separated by ';', and they will be treated as replicates that will be averaged
* **BEV_test** : test sample(s) for LFC calculation; if multiple BEV numbers are given, they should be separated by ';', and they will be treated as replicates that will be averaged

**Example input:**


| sg      | sgRNA_sequence       | translation_ref_seq                                  | BEV_start | BEV_end | primer        | frame | first_codon| last_codon| rev_com | BEV_ref | BEV_test |
| ------- | -------------------- | ---------------------------------------- |  -------: |  -----: | ------------- |  ----|----|---: | ------: | ------- | -------- |
| 397   | GTCACCCCTAAAGAGATCAT | tgtcttttctatgatctctttagGGGTGACCCAGTCTATT | 7         | 12      |F_C12_R_C12 |  2    |TGG|ATT| True    | 5;6     | 9;10     |


### Step 3: Create correlation input file

This file will be used to calculate Pearson correlations between replicates.

**Columns**: 

* **sg** : sg identifier 
* **reps_for_correlation** : semicolon-separated BEV numbers of which to calculate the pairwise Pearson correlation of the log-normalized read counts

**Example input:**
    
| sg      | reps_for_correlation |
| ------- | -------------------: | 
| 397     | 7;8 | 
| 397     | 9;10 | 
| 397     | 11;12 | 


### Step 4: Run allele frequencies notebook

Open 01_BEV_allele_frequencies.ipynb in the "notebooks" folder to begin.

### Step 5: Run nucleotide percentage notebook

Open 02_BEV_nucleotide_percentage_plots.ipynb in the "notebooks" folder to begin.

### Step 6: Run editing efficiency notebook

Open 03_BEV_editing_efficiency.ipynb the "notebooks" folder to begin.
