# Alignment Steps
To recap, the Lab Workflow has the following steps:

[Breadcrumb link](https://disulfidebond.github.io/alignment)

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

