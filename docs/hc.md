# Variant Calling with HaplotypeCaller Overview
To recap, the Lab Workflow has the following steps:

* [Alignment](https://disulfidebond.github.io/alignment)
* [Preprocessing/QAQC](https://disulfidebond.github.io/preprocess)
* [Variant Calling: HaplotypeCaller](https://disulfidebond.github.io/hc)
* [Variant Calling: Mutect2](https://disulfidebond.github.io/mutect2)

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
