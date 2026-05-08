# Transcriptome_general

Pipeline:
FastQC => Trimmomatic => Trinity Assembly => Trinotate Annotation
Turns out BaseSpace has the run already demultiplexed from the sample sheet. If you need to manually demultiplex, you can download the BaseSpace bcl-convert from typing in the instrument's serial number after searching for a bcl-convert download from Illumina. 

Samples are paired and locate in individual folders for each respective sample as forward-reads.fastq.gz and reverse-reads.fastq.gz (more complicated names irl) 

To find all fastq.gz files in a directory:
```
find . -name "*.fastq.gz"
```

FastQC is already installed on the HPC.
```
module avail fastqc

ml fastqc

mkdir fastqc_raw

find . -name "*.fastq.gz" | xargs fastqc -o fastqc_raw --threads 8
```

This will start recursively looking at each fastq.gz to run fastqc on, with outputs going into the new 'fastqc_raw' folder. 

You generally have to download the .html files to view them properly in a browser.

To remove the string of characters at the end of the individual sample file names (e.g. U001A1-ds.94d8e0cae4a34da2becc444efadb897a), use the following:
```
for f in *-ds.*
do
    mv "$f" "${f%-ds.*}"
done
```

### Trimmomatic

For individual reads (both forward and reverse) of each of the samples

First, make a directory for the trimmed reads and then batch samples (bash loop pattern) :) 

```
mkdir trimmed
```

```
#!/bin/bash
#SBATCH --job-name=trimmomatic_p1-100_individual
#SBATCH --partition=standard
#SBATCH --account=panilab
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=trimmomatic_individual_%j.out
#SBATCH --error=trimmomatic_individual_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bqc8em@virginia.edu

mkdir -p trimmed trimmed/unpaired

find . -name "*_R1_001.fastq.gz" | while read r1
do
    r2=${r1/_R1_/_R2_}

    sample=$(basename "$r1" _R1_001.fastq.gz)

    java -jar /standard/panilab/Becky/shortread-RNAseq/Trimmomatic/trimmomatic-0.39.jar PE \
      -threads 8 \
      "$r1" \
      "$r2" \
      trimmed/${sample}_paired_R1.fastq.gz \
      trimmed/unpaired/${sample}_unpaired_R1.fastq.gz \
      trimmed/${sample}_paired_R2.fastq.gz \
      trimmed/unpaired/${sample}_unpaired_R2.fastq.gz \
      ILLUMINACLIP:/standard/panilab/Becky/shortread-RNAseq/Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 \
      LEADING:3 \
      TRAILING:3 \
      SLIDINGWINDOW:4:15 \
      MINLEN:36
done

```
