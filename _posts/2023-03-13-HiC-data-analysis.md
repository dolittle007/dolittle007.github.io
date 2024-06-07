---
layout: post
title: "Hi-C data analysis"
date: 2022-03-13
category: tutorial
tags: [Hi-C, juicer, pairtools, mcool, analysis]
---

An Introduction to Hi-C data analysis.

<!--more-->

### Software preparation

```bash
conda install -c bioconda pairtools # install pairtools
conda install -c jrhawley mustache-hic
conda install -c bioconda hicexplorer
conda install -c bioconda pairix
conda install -c bioconda samtools
conda install -c anaconda openjdk # install java


# Install juicer tools [version 1.22.01](https://github.com/aidenlab/juicer/wiki/Download)
wget -c https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.22.01.jar
```



### fastq to pairs bam file

#### Step1: Alignment

Now that you have a genome file, index file and a reference fasta file you are all set to align your Micro-C library to the reference. Please note the specific settings that are needed to map mates independently and for optimal results with our proximity library reads.
| Parameter   |      Alignment function      |
|:------------|:----------------------------|
| mem    | set the bwa to use the BWA-MEM algorithm, a fast and accurate alignment algorithm optimized for sequences in the range of 70bp to 1Mbp |
| -5    |    for split alignment, take the alignment with the smallest coordinate (5’ end) as primary, the mapq assignment of the primary alignment is calculated independent of the 3’ alignment   |
| -S    | skip mate rescue |
|-P     | skip pairing; mate rescue performed unless -S also in use|
|-T0     |The T flag set the minimum mapping quality of alignments to output, at this stage we want all the alignments to be recorded and thus T is set up to 0, (this will allow us to gather full stats of the library, at later stage we will filter the alignments by mapping quality|
|-t|number of threads, default is 1. Set the numbers of threads to not more than the number of cores that you have on your machine (If you don’d know the number of cores, used the command lscpu and multiply Thread(s) per core x Core(s) per socket x Socket(s))|
|*.fasta or *.fa|Path to a reference file, ending with .fa or .fasta, e,g, hg38.fasta|
|*.fastq or *.fastq.gz|Path to two fastq files; path to read 1 fastq file, followed by fastq file of read 2 (usually labeled as R1 and R2, respectively). Files can be in their compressed format (.fastq.gz) or uncompressed (.fastq). In case your library sequence is divided to multiple fastq files, you can use a process substitution < with the cat command (see example below)|
|-o|sam file name to use for output results [stdout]. You can choose to skip the -o flag if you are piping the output to the next command using ‘\|’ |

Command:
```bash
bwa mem -5SP -T0 -t<threads> <ref.fasta> <HiC_R1.fastq> <HiC_R2.fastq> -o <aligned.sam>
```
Example (one pair of FASTQ file):
```bash
bwa mem -5SP -T0 -t16 hg38.fasta HiC_R1.fastq HiC_R2.fastq -o aligned.sam
```
Example (multiple pairs of FASTQ file):
```bash
bwa mem -5SP -T0 -t16 hg38.fasta <(zcat file1.R1.fastq.gz file2.R1.fastq.gz file3.R1.fastq.gz) <(zcat file1.R2.fastq.gz file2.R2.fastq.gz file3.R2.fastq.gz) -o aligned.sam
```

#### Step2: Recording valid ligation events

We use the __parse__ module of the __pairtools__ pipeline to find ligation junctions in Micro-C (and other proximity ligation) libraries. When a ligation event is identified in the alignment file the pairtools pipeline will record the outer-most (5’) aligned base pair and the strand of each one of the paired reads into __.pairsam__ file (pairsam format captures SAM entries together with the Hi-C pair information). In addition, it will also asign a pair type for each event. e.g. if both reads aligned uniquely to only one region in the genome, the type UU (Unique-Unique) will be assigned to the pair. The following steps are necessary to identify the high quality valid pairs over low quality events (e.g. due to low mapping quality):

__pairtools parse__ options

| Parameter   |Value|      Function                |
|:------------|:----|:-----------------------------|
|min-mapq | 40 | Mapq threshold for defining an alignment as a multi-mapping alignment. Alignment with mapq <40 will be marked as type M (multi)|
|walks-policy|5unique|Walks is the term used to describe multiple ligations events, resulting three alignments (instead of two) for a read pair. However, there are cases in which three alignment in read pairs are the result of one ligation event, pairtool parse can rescue this event. walks-policy is the policy for reporting un-rescuable walk. 5unique is used to report the 5’-most unique alignment on each side, if present (one or both sides may map to different locations on the genome, producing more than two alignments per DNA molecule)|

Command
```bash
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in <cores>\
--nproc-out <cores> --chroms-path <ref.genome> <aligned.sam> > <parsed.pairsam>
```
Example
```bash
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 8 --nproc-out 8 --chroms-path hg38.genome aligned.sam >  parsed.pairsam
```
At the parsing step, pairs will be flipped such that regardless of read1 and read2, pairs are always recorded with first side of the pair having the lower genomic coordinates.

#### Step3: Sorting the pairsam file
The parsed pairs are then sorted using pairtools sort
__pairtools sort__ options:
|Parameter|Function|
|:--------|:-------|
|--tmpdir|rovide a full path to a temp directory. A good rule of thumb is to have a space available for this directory at a volume of x3 of the overall volume of the fastq.gz files. Using a temp directory will help avoid memory issues|
|-nproc|Number of processes to split the sorting work|

Command
```bash
pairtools sort --nproc <cores> --tmpdir=<path/to/tmpdir> <parsed.pairsam> > <sorted.pairsam>
```
Example
```bash
pairtools sort --nproc 16 --tmpdir=/home/ubuntu/ebs/temp/  parsed.pairsam > sorted.pairsam
```
Please note that an absolute path for the temp directory is required for pairtools sort, e.g. path of the structure ~/ebs/temp/ or ./temp/ will not work, instead, something of this sort is needed /home/user/ebs/temp/


#### Step4: merging the pairsam file (optional)

Sequencing technical replicates (replicates within an experiment)
Sequencing replicates are merged after the alignment but before duplicate removal step (__pairsam dedup__), since reads resulting from a single PCR duplication event may exist in both replicates.
Merging is performed on (intermediate) pairsam files using __pairsam merge__.

Command
```bash
pairtools merge *.sorted.pairsam -o merged.pairsam
```
Example
```bash
pairtools merge *.sorted.pairsam -o merged.pairsam
```
#### Step5: Removig PCR duplicates

__pairtools dedup__ detects molecules that could be formed via PCR duplication and tags them as “DD” pair type. These pairs should be excluded from downstream analysis. Use the pairtools dedup command with the –output-stats option to save the dup stats into a text file.
__pairtools dedup__ options:
|Parameter|Function|
|:--------|:-------|
|--mark-dups|If specified, duplicate pairs are marked as DD in “pair_type” and as a duplicate in the sam entries|
|--output-stats|Output file for duplicate statistics. Please note that if a file with the same name already exists, it will be opened in the append mode|

Command:
```bash
pairtools dedup --nproc-in <cores> --nproc-out <cores> --mark-dups --output-stats <stats.txt> \
--output <dedup.pairsam> <sorted.pairsam>
```
Example:
```bash
pairtools dedup --nproc-in 8 --nproc-out 8 --mark-dups --output-stats stats.txt --output dedup.pairsam sorted.pairsam
```

#### Step6: Filtering
Sometimes you may need certain types of pairs based on their properties, such as mapq, pair type, distance or orientation. For all these manipulations, there is pairtools select which requires a file and pythonic condition as an input:

Command:
```bash
pairsamtools select <condition> -o <filtered_pairsam> <deduped_pairsam> 
```

Example:
```bash
pairsamtools select '(pair_type == "UU") or (pair_type == "UR") or (pair_type == "RU")' -o filtered.pairsam.gz deduped.pairsam.gz 
```
#### Step7: Generate .pairs and bam files

The __pairtools split__ command is used to split the final __.pairsam__ into two files: __.sam__ (or __.bam__) and __.pairs__ (__.pairsam__ has two extra columns containing the alignments from which the Hi-C pair was extracted, these two columns are not included in __.pairs__ files)



| Parameter | Function |
| :------ | :----------- |
| --output-pairs | Output pairs file. If the path ends with .gz or .lz4 the output is pbgzip-/lz4c-compressed. If you wish to pipe the command and output the pairs fils to stdout use - instead of file name |
| --output-sam | Output sam file. If the file name extension is .bam, the output will be written in bam format. If you wish to pipe the command, use - instead of a file name. please note that in this case the sam format will be used (and can be later converted to bam file e.g. with the command samtools view -bS -@16 -o temp.bam |

Command:
```bash
pairtools split --nproc-in <cores> --nproc-out <cores> --output-pairs <mapped.pairs> \
--output-sam <unsorted.bam> <filtered.pairsam>
```
Example:
```bash
pairtools split --nproc-in 8 --nproc-out 8 --output-pairs mapped.pairs --output-sam unsorted.bam dedup.pairsam
```
The __.pairs__ file can be used for generating contact matrix.



### References

* [4DN Hi-C data processing pipeline](https://data.4dnucleome.org/resources/data-analysis/hi_c-processing-pipeline)
* [Hi-C data analysis Bootcamp 2018](https://hms-dbmi.github.io/hic-data-analysis-bootcamp/#1)
* [Pairtools walkthrough](https://pairtools.readthedocs.io/en/latest/examples/pairtools_walkthrough.html)
* [Micro-C data analysis](https://micro-c.readthedocs.io/en/latest/fastq_to_bam.html)
* [Hi-C processing workflow based on 4DN consortium pipeline](https://gist.github.com/ggirelli/bbc52daad8f16564777785fbe98ec6ed)
* 
ut_prefix.stats.txt | pairtools select '(pair_type == "UU") or (pair_type == "UR") or (pair_type == "RU")' | pairtools split --nproc-in $cores --nproc-out $cores --output-pairs $output_prefix.pairs

bgzip --force $output_prefix.pairs
pairix -f $output_prefix.pairs.gz
```
```bash
sam2pair.sh test.sam test

This script will output test.pairs.gz.
```


### References

* [4DN Hi-C data processing pipeline](https://data.4dnucleome.org/resources/data-analysis/hi_c-processing-pipeline)
* [Hi-C data analysis Bootcamp 2018](https://hms-dbmi.github.io/hic-data-analysis-bootcamp/#1)
* [Pairtools walkthrough](https://pairtools.readthedocs.io/en/latest/examples/pairtools_walkthrough.html)
* [Micro-C data analysis](https://micro-c.readthedocs.io/en/latest/fastq_to_bam.html)
* [Hi-C processing workflow based on 4DN consortium pipeline](https://gist.github.com/ggirelli/bbc52daad8f16564777785fbe98ec6ed)
 
