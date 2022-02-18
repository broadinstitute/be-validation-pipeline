# Base Editor Validation (BEV) Pipeline

Welcome to the BEV pipeline! This analysis pipeline begins after you have received sequencing files back after conducting a base editor validation experiment. 

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

This file will be used in each of the three BEV notebooks. **Please refer to [Sequence_Orientation_Documentation.html](docs/Sequence_Orientation_Documentation.html) to ensure that all sequence inputs are in the correct orientation. 
<br/><br/>
**Columns**: 

* **sg** : sg identifier 
* **sgRNA_sequence** : sequence of sgRNA as designed 
* **translation_ref_seq**: reference sequence outputted by CRISPResso formatted such that any intronic sequences are lower-case, exons are upper-case, and UTRs are indicated by square brackets (if applicable) <u> must be sequence on strand that is being translated; may not necessarily be the same strand as the sgRNA sequence</u> 
    * Ex. <font color='grey'>tgtcttttctatgatctctttag</font><font color='green'>GGGTGACCCAGTCTATT</font>
* **BEV_start** : BEV number for first sample in sg
* **BEV_end** : BEV number for last sample in sg
* **primer** : name of primer pair (joined by '\_') used to amplify genomic locus as mentioned in sample name
    * Ex. <font color='purple'>F_C12</font><font color = 'blue'><b>_</b></font><font color='green'>R_C12</font>
* **frame** : frame for translation (manually determined for each sg / primer pair); position of first coding nucleotide in reference sequence within codon; frame can be 1, 2, 3
    * Ex. given reference sequence: tgtcttttctatgatctctttag<font color='green'>**G**</font>G|GTG|ACC|CAG|TCT|ATT 
        since the first coding nucleotide of the reference sequence (<font color='green'><b>G</b></font>) is the 2nd nucleotide in its codon 
        (\_<font color='green'><b>G</b></font>G) &rightarrow; frame = 2
* **first_codon** : first codon for translation 
* **last_codon** : last codon for translation 
* **rev_com** : samples for which reference sequence is on reverse strand 
* **BEV_ref** : reference sample(s) for log-fold change (LFC) calculation (i.e. early time point, empty vector, etc.); if multiple BEV numbers are given, they should be separated by ';', and they will be treated as replicates that will be averaged
* **BEV_test** : test sample(s) for LFC calculation; if multiple BEV numbers are given, they should be separated by ';', and they will be treated as replicates that will be averaged

**Example input:**


| sg      | sgRNA_sequence       | translation_ref_seq                                  | BEV_start | BEV_end | primer        | frame | first_codon| last_codon| rev_com | BEV_ref | BEV_test |
| ------- | -------------------- | ---------------------------------------- |  -------: |  -----: | ------------- |  ----|----|---: | ------: | ------- | -------- |
| 397   | GTCACCCCTAAAGAGATCAT | tgtcttttctatgatctctttagGGGTGACCCAGTCTATT | 7         | 12      |F_C12_R_C12 |  2    |TGG|ATT| True    | 5;6     | 9;10     |


### Step 2: Run allele frequencies notebook

Open [01_BEV_allele_frequencies.ipynb](notebooks/01_BEV_allele_frequencies.ipynb) to begin.

### Step 3: Run nucleotide percentage notebook

Open [02_BEV_nucleotide_percentage_plots.ipynb](notebooks/02_BEV_nucleotide_percentage_plots.ipynb) to begin.

### Step 4: Run editing efficiency notebook

Open [03_BEV_editing_efficiency.ipynb](notebooks/03_BEV_editing_efficiency.ipynb) to begin.
