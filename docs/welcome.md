# Alignment Steps
To recap, the Lab Workflow has the following steps:

* Alignment

* Preprocessing/QAQC

* Variant Calling

The final step, Variant Filtering and Classification, takes place strictly in VarSeq, and does not have any formal workflows in R, although supplementary R scripts have been created to help with analyzing the VarSeq data.

This R notebook describes how to run the Alignment Lab Workflow for single samples, as well as multiple samples. It is meant to be a guide and reference, and **not** simply run from start to finish.

For the Alignment step of the Lab Workflow, R is used as a wrapper. This means R itself does not run the analysis, but instead it facilitates running open source software tools like BWA that normally require Bash. The R function `system()` will primarily be used here, along with additional code that monitors for when the analysis is complete and if any errors appear.

## Alignment Single Sample

Before beginning, ensure the following are completed:

* The executable bwa is installed and on the $PATH

* The fastq read files are in the R working directory

* The genome reference file is in the same directory as your fastq files (and where R is running)

* The BWA index has been built, and is in the same directory as the fastq files (and where R is running)

Most likely all of the above have been completed already on Genie, but it's good practice to check ahead of time, similar to how you'd want to assemble all reagents and hardware before running a PCR gel.

The basic usage of BWA is bwa -t6 REFERENCE FASTQ1 FASTQ2 > OUTPUT.BAM
Where:

* bwa is the bwa executable, and -t6 indicates we will use 6 threads

* REFERENCE is the name of the fasta reference file. Usually this will be 'hg19.fa' for hg19 or 'GRCh38.d1.vd1.fa' for hg38

* FASTQ1 is read 1 of the fastq read pairs

* FASTQ2 is read 2 of the fastq read pairs

* '>' indicates we want bwa to write the output to a filename we provide

* OUTPUT.BAM is the filename we provide for the output file

An example could be `bwa -t6 hg19.fa somefastq1.fastq somefastq2.fastq > output.bam`

First, we're going to create variables that will hold the values for the bwa command. In the code block, fill in the correct values for the reference genome, fastq1, fastq2, and the output bam file, then run the code block.

```{r}
refGenome = c('hg19.fa') # this is the reference genome
fastq1 = c('tumor-0119_1401.TCGA-TS-A7PB-01A-11D-A34C-32_gdc_realn.1.fastq.gz')
fastq2 = c('tumor-0119_1401.TCGA-TS-A7PB-01A-11D-A34C-32_gdc_realn.2.fastq.gz')
outputBam = c('tumor-0119_1401.TCGA-TS-A7PB-01A-11D-A34C-32.bam')

# Validation code block
# check to make sure names were entered and the files exist
checkValuesEntered <- function(n) {
  if (n == '') {
    print('Please enter a name for the variable.')
    return(0)
  } else {
    if (!file.exists(n)) {
      print('Warning, unable to locate the file you entered!')
      print('Please check the name and ensure it is in the R working directory, which is:')
      print(getwd())
      return(0)
    } else {
      print('Value checked!')
      return(1)
    }
  }
}
checked = checkValuesEntered(refGenome)
checked = checkValuesEntered(fastq1)
checked = checkValuesEntered(fastq2)
refGenomeName = paste0('refGenome will be ', refGenome)
fastq1Name = paste0('fastq1 will be ', fastq1)
fastq2Name = paste0('fastq2 will be ', fastq2)
outputBamName = paste0('outputBam will be ', outputBam)
print(refGenomeName)
print(fastq1Name)
print(fastq2Name)
print(outputBamName)
# end validation code block
bwaCommand_singleSample = paste0('bwa mem -t4 ', refGenome, ' ', fastq1, ' ', fastq2, ' > ', outputBam)
```

When you're ready to run the BWA aligner, run the next code chunk. Note that you will not be able to run any commands or open any files in RStudio until the alignment finishes, or you stop the alignment.

```{r}
system(bwaCommand_singleSample)
```

## Alignment with Multiple Samples in Parallel

At a basic level, the command for running alignments of multiple samples in parallel is the same as running an alignment for a single sample. The setup for multiple samples in parallel is very similar to running a single sample, with the following modifications:

1. We need to identify and link together read pairs of fastq files for each alignment.

2. R and RStudio can only run a single task at a time, so we need to create/use a script of some kind that can run the alignment of different samples in parallel, which R or RStudio will then start and monitor.

3. We a way to separate out the number of jobs that are run simultaneously so that we do not overload the Genie server.

A bash script has been created named create_paw.bwa.sh, which incorporates all three of these modifications and creates a bash executable script that RStudio can run. Similar to the Preprocess workflow, you will need to provide create_paw.bwa.sh the alignment reference file name, and a text file with a comma-separated list of fastq file names.

The fastq file must have each read pair on one line, separated by a comma. An example is:

>sample1fastq1.fastq,sample1fastq2.fastq

>sample2fastq1.fastq,sample2fastq2.fastq

>sample3fastq1.fastq,sample3fastq2.fastq

>sample4fastq1.fastq,sample4fastq2.fastq

>sample5fastq1.fastq,sample5fastq2.fastq

Note that since we'll be using 4 threads (cores) for each alignment, running more than 6 alignments at a time is not recommended. If your list of fastq files has more than 6 entries, then you'll receive a warning telling you this.

In the following code chunk, fill in the reference file name and the fastq list file name with the corresponding variables.

```{r}

fastaName = c('hg19.fa')
fastqList = c('fastqList.txt') # name of text file with list of fastq names

# Validation code block
# check to make sure names were entered and the files exist
checkValuesEntered <- function(n) {
  if (n == '') {
    print('Please enter a name for the variable.')
    return(0)
  } else {
    if (!file.exists(n)) {
      print('Warning, unable to locate the file you entered!')
      print('Please check the name and ensure it is in the R working directory, which is:')
      print(getwd())
      return(0)
    } else {
      print('Value checked!')
      return(1)
    }
  }
}
checked = checkValuesEntered(fastaName)
checked = checkValuesEntered(fastqList)
# end validation code block
alignmentCommand = paste0('bash create_pawfile.bwa.sh ', fastaName, ' ' , fastqList)
system(alignmentCommand)

```

Once you've created the executable file, you should see the name of it appear. Copy and paste this filename into the next code chunk, then run the code chunk to start the parallel alignment step. 


```{r}
scriptName = c('bwa_align_0901_1459')
runParallelBWA = paste0('bash ', scriptName)
if (scriptName == '') {
  print('please enter the name of the executable file')
} else {
  system(runParallelBWA)
}
```

A minor difference between running parallel workflows and single workflows is you won't see any output in RStudio when the parallel workflow runs, except a message indicating the alignments have completed. Instead, a logfile has been created for each alignment, which will contain all of the alignment log info. It's unlikely you'll encounter errors, but if you do, then check the log corresponding to the sample with errors.


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


