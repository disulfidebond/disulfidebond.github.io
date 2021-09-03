# Overview
This site describes the steps in the Lab Workflow for the Jane Churpek Lab.

* Alignment Step
* Preprocess Step
* Variant Calling: HaplotypeCaller
* Variant Calling: Mutect2

# Preprocess Overview
To recap, the Lab Workflow has the following steps:

* Alignment

* Preprocessing/QAQC

* Variant Calling

The final step, Variant Filtering and Classification, takes place strictly in VarSeq, and does not have any formal workflows in R, although supplementary R scripts have been created to help with analyzing the VarSeq data.

This R notebook describes how to run the Preprocessing Lab Workflow for one or multiple samples, as well as instructions on running the necessary setup steps. It is meant to be a guide and reference, and **not** simply run from start to finish.

For the Preprocessing step of the Lab Workflow, R is used as a wrapper. This means R itself does not run the analysis, but instead it facilitates running open source software tools like Picard and GATK that normally require Bash. The R function `system()` will primarily be used here, along with additional code that monitors for when the analysis is complete and if any errors appear.

## Preprocess Setup

Before beginning, the following setup tasks need to be completed:

* Copy required files

  * The R working directory must have copies of the BAM files you would like to preprocess.

  * The R working directory must have copies of the BAM index files you would like to preprocess.

  * The R working directory must have a copy of the genomic fasta file that was used during the alignment step.

  * The R working directory must have the Read Group files for each BAM file being preprocessed. These files end in '.rg.txt' and can be found in the `/data/workspace_home/workspace_shared/` directory.

  * You have the known_indel files from 1000 genomes copied to the R working directory. These also can be found in the `/data/workspace_home/workspace_shared/` directory. See the Copy Known Indels section below for more information. 

* Create or copy files

  * You must either create or copy a fasta dictionary file from the genomic fasta file that was used during the alignment step. Note that this index file can be re-used as many times as needed, and only needs to be created once.

  * You must either create or copy a fasta index file from the genomic fasta file that was used during the alignment step. Note that this index file can be re-used as many times as needed, and only needs to be created once.

* Other

  * You must create a list of the BAM files you wish to preprocess. Only have at most 10 BAM files at a time in a list file. See the Additional Steps section below for more information.

  * The executable for gatk is installed and on the $PATH

  * The java .jar file for picard is installed and on the $PATH

  * The picard and GATK binaries should be installed and on the $PATH already. If they are not, contact John.
  
  * Copy the file `create_pawfile.preprocess.sh` to your R working directory

### Create Fasta Index File
To create a fasta index file, enter the name of the genome fasta file in the code chunk below, then run the code chunk.

```{r}

fastaFile = c('')
createFastaIdx = paste0('samtools faidx ', fastaFile)
system(createFastaIdx)
print('Completed setup task: create fasta index')
```

### Create Fasta Sequence Dictionary
To create a fasta sequence dictionary, run the next code chunk.
```{r}
createFastaDict = paste0('java -jar /data/workspace_home/workspace_shared/picard.jar CreateSequenceDictionary R=', fastaFile)
```

### Known Indels Files
There are two sets of known indel files: a hg38 version, and a hg19 liftover version.
A liftover in bioinformatics refers to a process whereby one reference genome is converted (or 'lifted over') to another one. In this case, the hg38 version of the human genome was as a template for mapping genomic coordinates onto the hg19 genome.


Both sets of files are in the `/data/workspace_home/workspace_shared/` directory

hg38:

* resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf

* resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf

hg19:

* resources_broad_hg38_v0_Homo_sapiens_assemblyhg19.crossMap.known_indels.fixed.sorted.vcf

* resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg19.crossMap.fixed.sorted.vcf

**Important**: Be sure to copy the index files (ending in '.idx') as well as the VCF files.

### Additional Steps
As mentioned above, you should only run the preprocessing steps for at most 10 samples at a time. The script you'll use to create the parallel processing takes as input a text list of BAM filenames, as well as the suffix for the read group files. For example, if you're using the TCGA file tumor-0119_1401, then the suffix would be '.rg.txt', which will point the parallel processor to the read group file tumor-0119_1401.rg.txt 

## Creating the Preprocess parallel script
Creating the parallel script is very similar to the Alignment workflow. In the code chunk below, enter the file name for the list of bam files, and enter the suffix for the read group identifiers, and indicate whether you are using the hg19 or hg38 reference genome.

Next, you'll see the names of three files in the output of RStudio:

* step1 is the read group correction step

* step2 is the mark duplicates step

* step3 is the BQSR step

```{r}
bamFileList = c('bamfileList.txt')
rgSuffix = c('.rg.txt')
refGenomeVersion = c('hg19')
runPAW_preprocessing = paste0('bash create_pawfile.preprocess.sh ', bamFileList, ' ', rgSuffix, ' ', refGenomeVersion, ' ', fastaFile)
system(runPAW_preprocessing)

```

## Running the Preprocess Script
Finally, enter (or copy and paste) the names for each step of the workflow in the variable names below, then start the parallel script by running the next code chunk

```{r}
step1Var = c('')
runStep1Var = paste0('bash ', step1Var)
step2Var = c('')
runStep2Var = paste0('bash ', step2Var)
step3Var = c('')
runStep3Var = paste0('bash ', step3Var)

system(runStep1Var)
Sys.sleep(10)
system(runStep2Var)
Sys.sleep(10)
system(runStep3Var)
print('parallel jobs completed')

```

# Variant Calling with HaplotypeCaller Overview
To recap, the Lab Workflow has the following steps:

* Alignment

* Preprocessing/QAQC

* Variant Calling

The final step, Variant Filtering and Classification, takes place strictly in VarSeq, and does not have any formal workflows in R, although supplementary R scripts have been created to help with analyzing the VarSeq data.

This R notebook describes how to run the Variant Calling Lab Workflow for one or multiple samples, as well as instructions on running the necessary setup steps. It is meant to be a guide and reference, and **not** simply run from start to finish.

For the Variant Calling step of the Lab Workflow, R is used as a wrapper. This means R itself does not run the analysis, but instead it facilitates running the open source software tool HaplotypeCaller from the Broad Institute, which normally requires Bash. The R function `system()` will primarily be used here, along with additional code that monitors for when the analysis is complete and if any errors appear.

## HaplotypeCaller Setup

Before beginning, the following setup tasks need to be completed:

* Copy required files

  * The R working directory must have copies of the BAM files you would like to analyze.

  * The R working directory must have copies of the BAM index files you would like to analyze.

  * The R working directory must have a copy of the genomic fasta file that was used during the alignment step.

* Create or copy files

  * You must either create or copy a fasta dictionary file from the genomic fasta file that was used during the alignment step. Note that this index file can be re-used as many times as needed, and only needs to be created once.

  * You must either create or copy a fasta index file from the genomic fasta file that was used during the alignment step. Note that this index file can be re-used as many times as needed, and only needs to be created once.

* Other

  * You must create a list of the BAM files you wish to analyze. Only have at most 10 BAM files at a time in a list file. See the Additional Steps section below for more information.

  * The executable for HaplotypeCaller is installed and on the $PATH. If it is not, contact John.
  
  * Copy the file `create_hcRunFile.sh` to your R working directory

## Create Fasta Index File
See [the section above](https://github.com/disulfidebond/disulfidebond.github.io/blob/gh-pages/docs/welcome.md#create-fasta-index-file) if you need to create a Fasta Index File. Otherwise, simply enter the name of the fasta reference file for the variable `fastaFile`

```{r}

fastaFile = c('')

```

## Create Fasta Sequence Dictionary
See [the section above](https://github.com/disulfidebond/disulfidebond.github.io/blob/gh-pages/docs/welcome.md#create-fasta-sequence-dictionary) if you need to create a Fasta Sequence Dictionary, otherwise skip this code chunk.

```{r}
createFastaDict = paste0('java -jar /data/workspace_home/workspace_shared/picard.jar CreateSequenceDictionary R=', fastaFile)
```


### Additional Steps for HaplotypeCaller 
HaplotypeCaller has several modes for running. This workflow will use the germline single sample workflow.
The same bamfile for preprocessing can be used for HaplotypeCaller, but be sure to update the bam file names to be the BQSR bam files. 

```{r}
# required format for bam list file.
# do not include '#' in the actual file

# normal_bamfile1.bam
# normal_bamfile2.bam

```


## Creating the parallel script for HaplotypeCaller 
Creating the parallel script is very similar to the Preprocessing workflow. In the code chunk below, enter the file name for the list of bam files, and indicate whether you are using the hg19 or hg38 reference genome.



```{r}
bamFileList = c('')
refGenomeVersion = c('hg38')
create_hc_runFile = paste0('bash create_hcRunFile.sh ', bamFileList, ' ', refGenomeVersion, ' ', fastaFile)

system(create_hc_runFile)

```

## Running the HaplotypeCaller Script
Finally, enter (or copy and paste) the name for the HaplotypeCaller runFile below, then start the parallel workflow by running the next code chunk.

```{r}
step1Var = c('')
runStep1Var = paste0('bash ', step1Var)
system(runStep1Var)
print('parallel jobs completed')

```

# Overview for Mutect2
To recap, the Lab Workflow has the following steps:

* Alignment

* Preprocessing/QAQC

* Variant Calling

The final step, Variant Filtering and Classification, takes place strictly in VarSeq, and does not have any formal workflows in R, although supplementary R scripts have been created to help with analyzing the VarSeq data.

This R notebook describes how to run the Variant Calling Lab Workflow for one or multiple samples, as well as instructions on running the necessary setup steps. It is meant to be a guide and reference, and **not** simply run from start to finish.

For the Variant Calling step of the Lab Workflow, R is used as a wrapper. This means R itself does not run the analysis, but instead it facilitates running the open source software tool Mutect2 from the Broad Institute, which normally requires Bash. The R function `system()` will primarily be used here, along with additional code that monitors for when the analysis is complete and if any errors appear.

# Setup

Before beginning, the following setup tasks need to be completed:

* Copy required files

  * The R working directory must have copies of the BAM files you would like to analyze.

  * The R working directory must have copies of the BAM index files you would like to analyze.

  * The R working directory must have a copy of the genomic fasta file that was used during the alignment step.

  * You have the known_indel files from 1000 genomes copied to the R working directory. These also can be found in the `/data/workspace_home/workspace_shared/` directory.  

* Create or copy files

  * You must either create or copy a fasta dictionary file from the genomic fasta file that was used during the alignment step. Note that this index file can be re-used as many times as needed, and only needs to be created once.

  * You must either create or copy a fasta index file from the genomic fasta file that was used during the alignment step. Note that this index file can be re-used as many times as needed, and only needs to be created once.

* Other

  * You must create a list of the BAM files you wish to analyze. Only have at most 10 BAM files at a time in a list file. See the Additional Steps section below for more information.

  * The executable for mutect2 is installed and on the $PATH. If it is not, contact John.
  
  * Copy the file `create_mutect2RunFile.sh` to your R working directory

## Create Fasta Index File
To create a fasta index file, enter the name of the genome fasta file in the code chunk below, then run the code chunk.

```{r}

fastaFile = c('')
createFastaIdx = paste0('samtools faidx ', fastaFile)
system(createFastaIdx)
print('Completed setup task: create fasta index')
```

## Create Fasta Sequence Dictionary
To create a fasta sequence dictionary, run the next code chunk.
```{r}
createFastaDict = paste0('java -jar /data/workspace_home/workspace_shared/picard.jar CreateSequenceDictionary R=', fastaFile)
```


## Additional Steps
Mutect2 can be run in two modes: tumor-only, where only tumor files are analyzed, and tumor-normal, where the normal sample is used as a basis to remove any germline variants. In the code chunks below, if you wish to run Mutect2 in tumor-normal mode, then nothing additional needs to be done, but if you wish to run Mutect2 in tumor-only mode, then comment out the 'Tumor-Only section' and uncomment the 'Tumor-Normal section.

Note also that the formatting for the bamfile is slightly different than preprocessing. The input text file must have one line per sample, with the bam file names and bam sample names.

```{r}
# required format for bam list file.
# the sample name can be retrieved from the bam file via samtools view -H bamfile.bam | grep '@RG'
# do not include '#' in the actual file
# An example tumor-normal bam file list is provided, followed by an example tumor-only bam file list


# normal_bamfile1.bam,sample1_normal_name,tumor_bamfile1.bam,sample1_tumor_sample_name
# normal_bamfile2.bam,sample2_normal_name,tumor_bamfile2.bam,sample2_tumor_sample_name

# ,,tumor_bamfile1.bam,
# ,,tumor_bamfile2.bam,

```


# Creating the parallel script
Creating the parallel script is very similar to the Preprocessing workflow. In the code chunk below, enter the file name for the list of bam files, and indicate whether you are using the hg19 or hg38 reference genome.



```{r}
bamFileList = c('')
refGenomeVersion = c('hg38')
create_Mutect2_runFiles = paste0('bash create_mutect2RunFile.tumor_normal.sh ', bamFileList, ' ', refGenomeVersion, ' ', fastaFile)
# If running tumor-only workflow, un-comment the next line and comment-out the previous line
# create_Mutect2_runFiles = paste0('bash create_mutect2RunFile.tumor_only.sh ', bamFileList, ' ', refGenomeVersion, ' ', fastaFile)
system(create_Mutect2_runFiles)

```

# Running the Script
Finally, enter (or copy and paste) the names for the Mutect2 steps of the workflow below. There will be two files, step1Var for the Mutect2 step and step2Var for the Filtering step. Finally, start the parallel workflow by running the next code chunk.

```{r}
step1Var = c('')
step2Var = c('')
runStep1Var = paste0('bash ', step1Var)
runStep2Var = paste0('bash ', step2Var)
system(runStep1Var)
Sys.sleep(2)
system(runStep2Var)
print('parallel jobs completed')

```


