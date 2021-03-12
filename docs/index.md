This repository helps set up an analysis pipeline to analyze base editor validation experiments using CRISPResso2 
(https://github.com/pinellolab/CRISPResso2).
<p>
<b>Step 1: Downloading & demultiplexing files</b><br>Once sequencing files are received from WalkUp, please use the 
"Setup.ipynb" notebook in the "notebooks" folder <a href="https://github.com/mhegde/be-validation-pipeline">here</a> on 
GitHub to download files, demultiplex files if necessary, make a backup and run FASTQC (quality check) on the files.
</p>

<p>
<b>Step 2: Accessing sequencing files for CRISPResso2</b><br>Once the setup notebook has been run the sequencing files 
should be available at <b><a href='https://drive.google.com/drive/folders/1uMXOLjvfY9TNlhwj0fVcTp6k2-heQe0c'>"GPP Cloud 
>R&D>Sequencing_Backup"</a></b>. They can be downloaded from here on to your local machine for running CRISPResso2. The
input files should all be in the folder where you will run CRISPResso2 from. It is recommended to create a new folder 
for every new run. The FASTQC outputs can be found here <b><a href='https://drive.google.com/drive/folders/1NhOYq3_P2Jr3aj_K0KbU1f-9iSWBbKso'>
"GPP Cloud >R&D>FASTQC_outputs"</a></b>.
</p>

<p>
<b>Step 3: Installing Docker to run CRISPResso2</b><br>Installing CRISPResso via Docker is the easiest way to use
it. The way to do this is explained <a href="https://github.com/pinellolab/CRISPResso2#docker" target="_blank">here</a>.
</p>

<p>
<b>Step 4: Setting up a batch file for CRISPResso2</b><br>Batch file for a CRISPResso2 run is a .txt file with various 
columns indicating different parameters for each sample in the run. This batch file should also be placed in the folder 
created for the CRISPResso2 run. The different columns in a batch file are as described below: 
<ol>
    <li><b>name</b>: Generate a name of the format, <font color="#8b0000"> BEV_samplenumber_primerpair</font>
    for each sample. For example, <font color="blue">the name for sample 1 would be BEV_001_F1_A1_R1_A1, where 001 is 
    the sample number and F1_A1 is the forward primer and R1_A1 is the reverse primer</font>.</li>
    <li><b>fastq_r1</b>: Path to relevant demultiplexed FASTQ input files. It is recommended that you have all FASTQ
    inputs in a single folder.</li>
    <li><b>amplicon_seq</b>: Relevant reference sequence for each sample. It is recommended to include the forward 
    and the reverse primers in the reference sequence.</li>
    <li><b>guide_seq</b>: Sequence of sgRNA for each sample.</li>
    <li><b>Other parameters</b>: The following parameters can be included as columns in the batch file if they are 
    different for each sample. If they are the same for every sample, they can just be included in the syntax while 
    running CRISPResso. For more information on what each of these parameters do, please read the documentation 
    <a href="https://github.com/pinellolab/CRISPResso2" target="_blank">here</a>. </li>
    <ul>
        <li><b>-w or --quantification_window_size or --window_around_sgrna </b>: (default: 1) For base editors, we 
        recommend using <font color="blue"> 20 </font>. Please bear in mind that a wider window results in more reads 
        with sequencing errors being classified as edited. For knockout experiments, the recommendation is to use the 
        default value.</li>
        <li><b>-wc or --quantification_window_center or --cleavage_offset</b>: (default: -3) For base editors, we 
        recommend using the center of the guide, i.e. <font color="blue"> -10 </font>. The default value, i.e. the cut 
        position, can be used for SpCas9 knockout experiments. This parameter can be varied accordingly based on the 
        cut position of the nuclease used.</li>
        <li><b>--exclude_bp_from_left</b>: Due to the presence of stagger sequences in our vectors, reads need to be 
        trimmed and the provided trimming adapters cannot be used. This parameter can be toggled appropriately depending
        on the position of the sgRNA in the amplicon sequence to trim the reads. </li>
        <li><b>--exclude_bp_from_right</b>: This argument can be used to trim reads that have low sequencing quality at 
        at the ends.</li>
    </ul>
</ol>
<p><b>Sample batch file shown below:</b></p>
<table>
    <tbody>
        <tr>
            <th><b>name</b></th>
            <th><b>fastq_r1</b></th>
            <th><b>amplicon_seq</b></th>
            <th><b>guide_seq</b></th>
            <th><b>w</b></th>
            <th><b>wc</b></th>
        </tr>
        <tr>
            <td>BEV_010_F1_A1_R1_A1</td>
            <td>validation-inputs/BEV/Plate1/BEV_010_F1_A1_R1_A1.construct.fastq.gz</td>
            <td>GCTATTTAGTGTTATCCAAGGAACATCTTCAGTATCTCTAGGATTCTCTGAGCATGGCAGTTTCTGCTTAT</td>
            <td>GGAACATCTTCAGTATCTCT</td>
            <td>20</td>
            <td>-10</td>
        </tr>
        <tr>
            <td>BEV_016_F1_A2_R1_A2</td>
            <td>validation-inputs/BEV/Plate1/BEV_016_F1_A2_R1_A2.construct.fastq.gz</td>
            <td>TTATATACCTTTTGGTTATATCATTCTTACATAAAGGACACTGTGAAGGCCCTTTCTTCTGGTTGAGAA</td>
            <td>GTTATATCATTCTTACATAA</td>
            <td>1</td>
            <td>-3</td>
        </tr>
    </tbody>
</table>
<p>The sample "BEV_010_F1_A1_R1_A1" is a base editor sample and "BEV_016_F1_A2_R1_A2" is an SpCas9 knockout sample as 
indicated by the "w" and "wc" columns.</p>

<p>
<b>Step 5: Running CRISPResso2</b><br>To run CRISPResso2,
<ol>
    <li>Make sure Docker is running on your computer. If it is, you should see the Docker logo on your menu bar which 
    when clicked on should say "Docker Desktop is running".</li> 
    <li>Make sure you are running CRISPResso2 in the folder with the batch file. This can be done using the cd command.
    The outputs will be generated in this folder as well. </li>
    <li>Finally, open your terminal and type the following command:<br>
    <b>docker run -v ${PWD}:/DATA -w /DATA -i pinellolab/crispresso2 CRISPRessoBatch --batch_settings [batch file name] 
    --skip_failed --suppress-plots</b></li>
</ol>

<p>
<b>Few additional tips/tricks:</b>
<ul>
    <li>For experiments with both knockout and base editor samples, the quantification window size and window center 
    parameters can be toggled appropriately to run all samples together.</li>
    <li>--skip_failed is a good argument to use while running batch mode. This argument makes sure the job does not quit
    if one sample fails. The failed sample can be identified in the RUNNING_LOG file. </li>
    <li>One common error you may encounter is <b>"CRISPResso batch #x was killed by your system. Please decrease the 
    number of processes (-p) and run again."</b> Using  the --suppress-plots or --suppress-report arguments might fix this 
    error.  If not, please try increasing the number of resources such as "CPUs" and "Memory" in your Docker 
    preferences.</li>
</ul>

<p>
<b>Step 6: Analyzing CRISPResso2 results</b><br> Please use the "BEV_allele_frequencies_v4.ipynb" notebook in the 
"notebooks" folder <a href="https://github.com/mhegde/be-validation-pipeline">here</a> to further analyze your data.  
</p>