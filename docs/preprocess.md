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
