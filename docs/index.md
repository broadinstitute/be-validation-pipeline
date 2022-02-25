Steps to set up an analysis pipeline to analyze base editor validation experiments using [CRISPResso2](https://github.com/pinellolab/CRISPResso2): <br/>
    1. Using the [BEV tool](https://gpplims.broadinstitute.org/screening/bev/create) on GPP LIMS (available to GPP only) ([Jump to steps](#running-bev-tool-on-gpp-lims)) <br/>
    2. On your local machine ([Jump to steps](#running-crispresso2-on-local-machine))

## Running BEV tool on GPP LIMS

**Step 1: Sample name**

Assign a descriptive name for the samples, without spaces and beginning with your initials (similar to PoolQ). 

**Step 2: Construct file and Barcode file**

* If sequencing available through WalkUp sequencing, you can use the Walk-Up Sample Finder to access sequencing files by inputting the sample ID and lane number.
* If sequencing files are no longer available through WalkUp, the backed up files can also be accessed by providing a file path 

**Step 3: Conditions file**

A .csv file with two columns without headers: <br/>
    1. Barcode sequences <br/>
    2. Descriptive condition names without spaces or special characters (e.g. A549_RDA569_RepA_D7_Dropout) because these sample names will be used to name the demultiplexed FASTQ files. 

The BEV tool demultiplexes files for you. 

**Step 4: Batch file**

Batch file for a CRISPResso2 run is a .txt file with various columns indicating different parameters for each sample 
in the run. The batch file for running CRISPResso2 on BEV is similar to the one used for running CRISPResso2 on a local machine, with some differences (bolded). **Please refer to the "[Sequence_Orientation_Documentation.html](Sequence_Orientation_Documentation.html)" to input the appropriate amplicon and guide sequences.** 

The different columns in a batch file are as described below:
1. `name`: Generate a name of the format, <font color="#8b0000"> BEV_samplenumber_primerpair</font>
    for each sample. For example, <font color="blue">the name for sample 1 would be BEV_001_F1_A1_R1_A1, where 001 is 
    the sample number and F1_A1 is the forward primer and R1_A1 is the reverse primer</font>.
2. **`identifier`**: Condition name from conditions file 
3. `amplicon_seq`: Relevant amplicon sequence for each sample. It is recommended to include the forward and the 
   reverse primers in the amplicon sequence. For more information on how to obtain the right amplicon sequence, refer 
   to  **"[Sequence_Orientation_Documentation.html](Sequence_Orientation_Documentation.html)"**.
4. `guide_seq`: Sequence of sgRNA for each sample.
5. **Other parameters**: The following parameters can be included as columns in the batch file if they are 
    different from the default values of CRISPResso2 which can be found in the [documentation](https://github.com/pinellolab/CRISPResso2). 
    running CRISPResso. For more information on what each of these parameters do, please read the documentation 
    [here](https://github.com/pinellolab/CRISPResso2).
    - **`-w` or `--quantification_window_size` or `--window_around_sgrna`**: (default: 1) For base editors, we 
        recommend using <font color="blue"> 20 </font>. Please bear in mind that a wider window results in more reads 
        with sequencing errors being classified as edited. For knockout experiments, the recommendation is to use the 
        default value.
    - **`-wc` or `--quantification_window_center` or `--cleavage_offset`**: (default: -3) For base editors, we 
        recommend using the center of the guide, i.e. <font color="blue"> -10 </font>. The default value, i.e. the cut 
        position, can be used for SpCas9 knockout experiments. This parameter can be varied accordingly based on the 
        cut position of the nuclease used. 
    - **`--exclude_bp_from_left`**: Due to the presence of stagger sequences in our vectors, reads need to be 
        trimmed and the provided trimming adapters cannot be used. To trim the reads on the left end this parameter can 
        be toggled appropriately depending on the position of the sgRNA in the amplicon sequence.
    - **`--exclude_bp_from_right`**: This argument can be used to trim reads that have low sequencing quality at the ends.

    **Sample batch file shown below:**    

    **name**|**fastq_r1**|**amplicon_seq**|**guide_seq**|**w**|**wc**|**exclude_bp_from_left**|**exclude_bp_from_right**|**plot_window_size**|
    :-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|
    BEV_417_F3_R2|A549_RDA569_RepA_D7_Dropout_|GCAGAAAGTCAGTCCCATGGAATTTTCGCTTCCCACAGGTCTCTGCTAGGGGGCTGGGGTTGGGGTGGGGGTGGTGGGCCTGCCCTTCCAATGGATCCACTCACAGTTTCCATAGGTCTGAAAATGTTTCCTGACTCAGAGGGGGCTCGACGCTAGGATCTGACTGCGGCTCCTCCATGGCAGTGACCCGGAAGGCAGTCTGGCTGCTGCAAGAGGAAAAGTGGGGATCCA|GCTCCTCCATGGCAGTGACC|20|-10|8|8|40|


**Step 4: Running CRISPResso2**

Enter your email address for the notification email, then click the button that says "Run CRISPResso2." Similar to PoolQ, you will get a notification if your BEV run dies or succeeds. 
   
   **Few additional tips/tricks:**

   + For experiments with both knockout and base editor samples, the quantification window size and window center 
    parameters can be toggled appropriately to run all samples together. 
  + Note: you should download the CRISPResso output from BEV, and look through the RUNNING_LOG file in the "crispresso_output" folder to make sure all your samples were aligned and analyzed properly.
+ If you are encountering an error where CRISPResso2 is not able to align any reads to your reference sequence, 
    first check your reference sequence and then consider changing `--default_min_aln_score` argument to <font color='blue'>50</font>. The default value is 60.

**Step 5: Analyzing CRISPResso2 results**

Please use the **"01_BEV_allele_frequencies.ipynb"**, **"02_BEV_nucleotide_percentage_plots.ipynb"** and **"03_BEV_editing_efficiency.ipynb"**. 
notebooks in the "notebooks" folder [here](https://github.com/mhegde/be-validation-pipeline) to further analyze your data.  


## Running CRISPResso2 on local machine

**Step 1: Downloading & demultiplexing files**

* If demultiplexed files were requested from WalkUp sequencing, the received files can be used directly with CRISPResso2. 
If not, please demultiplex the files such that all the reads for a sample are in a single fastq file. 
* Please make sure these demultiplexed files are downloaded onto your local machine before running CRISPResso2. 
* Create a folder to save all the files relevant to the CRISPResso2 run. Within the main folder another folder can be
created to save the input fastq files as shown. 

**Step 2: Installing Docker to run CRISPResso2**

Running CRISPResso2 via Docker is the easiest way to use it. The way to do this is explained 
[here](https://github.com/pinellolab/CRISPResso2#docker).

**Step 3: Setting up a batch file for CRISPResso2**

Batch file for a CRISPResso2 run is a .txt file with various columns indicating different parameters for each sample 
in the run. This batch file should also be placed in the folder created for the CRISPResso2 run. 
**Please refer to the "[Sequence_Orientation_Documentation.html](Sequence_Orientation_Documentation.html)" to input the appropriate amplicon and guide sequences.** 

The different columns in a batch file are as described below:
1. `name`: Generate a name of the format, <font color="#8b0000"> BEV_samplenumber_primerpair</font>
    for each sample. For example, <font color="blue">the name for sample 1 would be BEV_001_F1_A1_R1_A1, where 001 is 
    the sample number and F1_A1 is the forward primer and R1_A1 is the reverse primer</font>.
2. `fastq_r1`: Path to relevant demultiplexed FASTQ input files. It is recommended that you have all FASTQ
    inputs in a single folder.
3. `amplicon_seq`: Relevant amplicon sequence for each sample. It is recommended to include the forward and the 
   reverse primers in the amplicon sequence. For more information on how to obtain the right amplicon sequence, refer 
   to  **"[Sequence_Orientation_Documentation.html](Sequence_Orientation_Documentation.html)"**.
4. `guide_seq`: Sequence of sgRNA for each sample.
5. **Other parameters**: The following parameters can be included as columns in the batch file if they are 
    different for each sample. If they are the same for every sample, they can just be included in the syntax while 
    running CRISPResso. For more information on what each of these parameters do, please read the documentation 
    [here](https://github.com/pinellolab/CRISPResso2).
    - **`-w` or `--quantification_window_size` or `--window_around_sgrna`**: (default: 1) For base editors, we 
        recommend using <font color="blue"> 20 </font>. Please bear in mind that a wider window results in more reads 
        with sequencing errors being classified as edited. For knockout experiments, the recommendation is to use the 
        default value.
    - **`-wc` or `--quantification_window_center` or `--cleavage_offset`**: (default: -3) For base editors, we 
        recommend using the center of the guide, i.e. <font color="blue"> -10 </font>. The default value, i.e. the cut 
        position, can be used for SpCas9 knockout experiments. This parameter can be varied accordingly based on the 
        cut position of the nuclease used. 
    - **`--exclude_bp_from_left`**: Due to the presence of stagger sequences in our vectors, reads need to be 
        trimmed and the provided trimming adapters cannot be used. To trim the reads on the left end this parameter can 
        be toggled appropriately depending on the position of the sgRNA in the amplicon sequence.
    - **`--exclude_bp_from_right`**: This argument can be used to trim reads that have low sequencing quality at the ends.

    **Sample batch file shown below:**    

    **name**|**fastq_r1**|**amplicon_seq**|**guide_seq**|**w**|**wc**
    :-----:|:-----:|:-----:|:-----:|:-----:|:-----:
    BEV_010_F1_A1_R1_A1|validation-inputs/BEV/Plate1/BEV_010_F1_A1_R1_A1.construct.fastq.gz|GCTATTTAGTGTTATCCAAGGAACATCTTCAGTATCTCTAGGATTCTCTGAGCATGGCAGTTTCTGCTTAT|GGAACATCTTCAGTATCTCT|20|-10
    BEV_016_F1_A2_R1\_A2|validation-inputs/BEV/Plate1/BEV_016_F1_A2_R1_A2.construct.fastq.gz|TTATATACCTTTTGGTTATATCATTCTTACATAAAGGACACTGTGAAGGCCCTTTCTTCTGGTTGAGAA|GTTATATCATTCTTACATAA|1|-3

    The sample "BEV_010_F1_A1_R1_A1" is a base editor sample and "BEV_016_F1_A2_R1_A2" is an SpCas9 knockout sample as 
    indicated by the "w" and "wc" columns.

**Step 4: Running CRISPResso2**

To run CRISPResso2,

1. Make sure Docker is running on your computer. If it is, you should see the Docker logo on your menu bar which 
    when clicked on should say "Docker Desktop is running". 
2. Make sure you are running CRISPResso2 in the folder with the batch file. This can be done using the cd command. 
   The outputs will be generated in the same folder as well. 
3. If you are running CRISPResso2 exclusively on base editing samples, please use the --base_edit argument as indicated 
   below.
4. Finally, open your terminal and type the following command:
   
    `docker run -v ${PWD}:/DATA -w /DATA -i pinellolab/crispresso2 CRISPRessoBatch --batch_settings [batch file name] 
    --skip_failed --base_edit`
   
    **Few additional tips/tricks:**

    + For experiments with both knockout and base editor samples, the quantification window size and window center 
    parameters can be toggled appropriately to run all samples together.
    + --skip_failed is a good argument to use while running CRISPResso2 in batch mode. This argument makes sure the job does
    not quit if one sample fails. The failed sample can be identified in the RUNNING_LOG file.
    + One common error you may encounter is **"CRISPResso batch #x was killed by your system. Please decrease the 
    number of processes (-p) and run again."** Using  the `--suppress-plots` or `--suppress-report` arguments might help fix 
    this error.  If not, please try increasing the number of resources such as "CPUs" and "Memory" in your Docker 
    preferences.
    + If you are encountering an error where CRISPResso2 is not able to align any reads to your reference sequence, 
    first check your reference sequence and then consider changing `--default_min_aln_score` argument to <font color='blue'>50</font>. The default value is 60.

**Step 5: Analyzing CRISPResso2 results**

Please use the **"01_BEV_allele_frequencies.ipynb"**, **"02_BEV_nucleotide_percentage_plots.ipynb"** and **"03_BEV_editing_efficiency.ipynb"**.
notebooks in the "notebooks" folder [here](https://github.com/mhegde/be-validation-pipeline) to further analyze your data.  

