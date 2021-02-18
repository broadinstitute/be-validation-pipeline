# Base editor validation pipeline

## Base editor screen validation analysis pipeline using CRISPResso2 
This repository helps set up an analysis pipeline to analyze base editor validation experiments using CRISPResso2.
<p>
<b>Step 1: Downloading & demultiplexing files</b><br>Once sequencing files are received from WalkUp, please use the 
"Setup.ipynb" notebook to download files, demultiplex files if necessary, make a backup and run FASTQC (quality check) 
on the files.
</p>

<p>
<b>Step 2: Accessing sequencing files for CRISPResso2</b><br>Once the setup notebook has been run the sequencing files 
should be available at <b><a href='https://drive.google.com/drive/folders/1uMXOLjvfY9TNlhwj0fVcTp6k2-heQe0c'>"GPP Cloud 
>R&D>Sequencing_Backup"</a></b>. They can be downloaded from here on to your local machine for running CRISPResso2. 
The FASTQC outputs can be found here <b><a href='https://drive.google.com/drive/folders/1NhOYq3_P2Jr3aj_K0KbU1f-9iSWBbKso'>
"GPP Cloud >R&D>FASTQC_outputs"</a></b>.
</p>

<p>
<b>Step 3: Setting up a batch file for CRISPResso2</b><br>Batch file for a CRISPResso2 run is a .txt file with various 
columns indicating different parameters for each sample in the run. The different columns are as described below: 
<ol>
    <li>Generate a name of the format, "BEV_sample_number_primerpair" 
for each sample. For example, the name for sample 1 would be BEV_001_F1_A1_R1_A1, where 001 is the sample number and 
F1_A1 is the forward primer and R1_A1 is the reverse primer.



