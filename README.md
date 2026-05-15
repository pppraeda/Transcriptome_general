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
Needed to pull out the one ghost sample from the rest of the 23 ghost samples, so I lined up the well numbers with the order of samples loaded for sequencing. The pancakes were in well 20/sample 20, so I moved those left and right reads to their own directory (Pancake-trimmed), and moved all of the other left/right reads to Ghost-trimmed. 
To concatenate the 23 pairs of ghost reads:
```
zcat *_paired_R1.fastq.gz | gzip > ghost_left.fq.gz

zcat *_paired_R2.fastq.gz | gzip > ghost_right.fq.gz

```
This will take some time but can be done locally. Ended up with left reads of 2.23 Gb and right reads of 2.34 Gb. 

### Assembly with Trinity
I already have a conda environment set up with all the requirements for Trinity, so i'll be using that for assemblies. 

P1-100 Ghost assembly:

```
#!/bin/bash
#SBATCH --partition=standard
#SBATCH --account=panilab
#SBATCH --job-name=trinity_ghost-p1-100
#SBATCH --output=trinity.out
#SBATCH --error=trinity.err
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --mail-user=bqc8em@virginia.edu
#SBATCH --mail-type=END,FAIL

ml miniforge

conda activate /home/bqc8em/.conda/envs/Trinity


Trinity \
  --seqType fq \
  --left /standard/panilab/Becky/shortread-RNAseq/P1-100/BaseSpace/WorleyLab-442058840/BCLConvert_02_13_2026_02_18_28Z-899399509/trimmed/Ghost-trimmed/ghost_left.fq.gz \
  --right /standard/panilab/Becky/shortread-RNAseq/P1-100/BaseSpace/WorleyLab-442058840/BCLConvert_02_13_2026_02_18_28Z-899399509/trimmed/Ghost-trimmed/ghost_right.fq.gz \
  --CPU 16 \
  --max_memory 100G \
  --normalize_reads
```

Pancake P1-100 assembly:
```
#!/bin/bash
#SBATCH --partition=standard
#SBATCH --account=panilab
#SBATCH --job-name=trinity_pancake
#SBATCH --output=trinity.out
#SBATCH --error=trinity.err
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --mail-user=bqc8em@virginia.edu
#SBATCH --mail-type=END,FAIL

ml miniforge

conda activate /home/bqc8em/.conda/envs/Trinity


Trinity \
  --seqType fq \
  --left /standard/panilab/Becky/shortread-RNAseq/P1-100/BaseSpace/WorleyLab-442058840/BCLConvert_02_13_2026_02_18_28Z-899399509/trimmed/Pancake-trimmed/U020D3_S20_L001_paired_R1.fastq.gz \
  --right /standard/panilab/Becky/shortread-RNAseq/P1-100/BaseSpace/WorleyLab-442058840/BCLConvert_02_13_2026_02_18_28Z-899399509/trimmed/Pancake-trimmed/U020D3_S20_L001_paired_R2.fastq.gz \
  --CPU 16 \
  --max_memory 100G \
  --normalize_reads
```

-----------------------------

### P4-300 Trimming

Navigate to the proper directory containing the data for P4-300. This data had to be moved to an external hard drive then to HPC due to something with BaseSpace (wasn't being easy for some reason). For the Undetermined files that exceeded 10 Gb (maximum upload size for HPC, apparently), I used FileZilla to get them onto the Panilab folder just in case we needed them (AKA so all data is present, even though I will move forward with demultiplexed data and paired reads only). 

```
#!/bin/bash
#SBATCH --job-name=trimmomatic_p4-300_individual
#SBATCH --partition=standard
#SBATCH --account=panilab
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=trimmomatic_individual_P4-300_%j.out
#SBATCH --error=trimmomatic_individual_P4-300_%j.err
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

I moved sample 20 (Pancake) into a new directory called Pancake and did the same with other reads (into Ghost). Assembly scripts are below.


Pancake: 

```
#!/bin/bash
#SBATCH --partition=standard
#SBATCH --account=panilab
#SBATCH --job-name=trinity_pancake-p4-300
#SBATCH --output=trinity.out
#SBATCH --error=trinity.err
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --mail-user=bqc8em@virginia.edu
#SBATCH --mail-type=END,FAIL

ml miniforge

conda activate /home/bqc8em/.conda/envs/Trinity


Trinity \
  --seqType fq \
  --left /standard/panilab/Becky/shortread-RNAseq/P4-300/BaseSpace/WorleyLab-442058840/BCLConvert_03_12_2026_12_03_44Z-906926024/trimmed/Pancake/U020D3_S20_L001_paired_R1.fastq.gz \
  --right /standard/panilab/Becky/shortread-RNAseq/P4-300/BaseSpace/WorleyLab-442058840/BCLConvert_03_12_2026_12_03_44Z-906926024/trimmed/Pancake/U020D3_S20_L001_paired_R2.fastq.gz\
  --CPU 16 \
  --max_memory 100G \
  --normalize_reads
```

For Ghost, I needed to concatenate all left and right reads, respectively, again using the following commands while in the proper directory (.../trimmed/Ghost):

```
zcat *_paired_R1.fastq.gz | gzip > ghost_left.fq.gz

zcat *_paired_R2.fastq.gz | gzip > ghost_right.fq.gz

```

Assembly
```
#!/bin/bash
#SBATCH --partition=standard
#SBATCH --account=panilab
#SBATCH --job-name=trinity_pancake-p4-300
#SBATCH --output=trinity.out
#SBATCH --error=trinity.err
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --mail-user=bqc8em@virginia.edu
#SBATCH --mail-type=END,FAIL

ml miniforge

conda activate /home/bqc8em/.conda/envs/Trinity


Trinity \
  --seqType fq \
  --left /standard/panilab/Becky/shortread-RNAseq/P4-300/BaseSpace/WorleyLab-442058840/BCLConvert_03_12_2026_12_03_44Z-906926024/trimmed/Ghost/ghost_left.fastq.gz \
  --right /standard/panilab/Becky/shortread-RNAseq/P4-300/BaseSpace/WorleyLab-442058840/BCLConvert_03_12_2026_12_03_44Z-906926024/trimmed/Ghost/ghost_right.fastq.gz\
  --CPU 16 \
  --max_memory 128G \
  --normalize_reads
```
---------------------------

### Check BUSCO scores of assemblies: 
```
  busco \
  -i /standard/panilab/Becky/shortread-RNAseq/P1-100/BaseSpace/WorleyLab-442058840/BCLConvert_02_13_2026_02_18_28Z-899399509/trimmed/Ghost-trimmed/trinity_out_dir/Trinity-trim-ghost.fasta \
  -l /scratch/bqc8em/busco_downloads/lineages/eukaryota_odb10 \
  -o P1-100_busco_output \
  -m transcriptome \
  -c 16
```
----------------------------
GENERAL FOR BLASTING ASSEMBLIES


Make BLAST database of assemblies:
```
  makeblastdb \
  -in Trinity.fasta \
  -dbtype nucl \
  -out trinity_db
```
Blast the database:
```
  blastn \
  -query query.fasta \
  -db trinity_db \
  -out results.out \
  -outfmt 6 \
  -evalue 1e-5 \
  -num_threads 16
```


Ghost P1-100: 
Path: /standard/panilab/Becky/shortread-RNAseq/P1-100/BaseSpace/WorleyLab-442058840/BCLConvert_02_13_2026_02_18_28Z-899399509/trimmed/Ghost-trimmed/trinity_out_dir
```
cd /standard/panilab/Becky/shortread-RNAseq/P1-100/BaseSpace/WorleyLab-442058840/BCLConvert_02_13_2026_02_18_28Z-899399509/trimmed/Ghost-trimmed/trinity_out_dir

  makeblastdb \
  -in /standard/panilab/Becky/shortread-RNAseq/P1-100/BaseSpace/WorleyLab-442058840/BCLConvert_02_13_2026_02_18_28Z-899399509/trimmed/Ghost-trimmed/trinity_out_dir/Trinity-trim-ghost.fasta \
  -dbtype nucl \
  -out trinity_ghost_100_db
```

IMPORTANT: Make sure you're using tblastn to blast protein sequences against the transcriptome! 
```
tblastn \
  -query /standard/panilab/Becky/shortread-RNAseq/P1-100/BaseSpace/WorleyLab-442058840/BCLConvert_02_13_2026_02_18_28Z-899399509/trimmed/Ghost-trimmed/trinity_out_dir/BLAST/piwi-1_hmiamia-prot.fasta \
  -db trinity_ghost_100_db \
  -out 	piwi-1_ghost \
  -outfmt 6
```

Interpreting the outputs
Example output (piwi-1 from Hofstenia miamia) 
                                                                        align   mis-    gap     query          subject
Query sequence ID                   Subject (Trinity transcr.)  % ID    length  match   open    start   end    start   end     evalue      bitscore

lcl|KJ658709.1_prot_AID23631.1_1	TRINITY_DN1954_c0_g1_i4	    34.173	834	    505	    19	    1	    813	    2938	506	    6.02e-158	490
lcl|KJ658709.1_prot_AID23631.1_1	TRINITY_DN22208_c0_g1_i2	32.360	788	    497	    18	    42	    813	    719	    3022	2.07e-123	398
lcl|KJ658709.1_prot_AID23631.1_1	TRINITY_DN22208_c0_g1_i3	32.360	788	    497	    18	    42	    813	    719	    3022	2.33e-123	398
lcl|KJ658709.1_prot_AID23631.1_1	TRINITY_DN22208_c0_g1_i1	32.360	788	    497	    18	    42	    813	    719	    3022	4.72e-123	398
lcl|KJ658709.1_prot_AID23631.1_1	TRINITY_DN8162_c0_g1_i2	    46.352	233	    115	    7	    547	    774	    35	    718	    1.46e-54	189
lcl|KJ658709.1_prot_AID23631.1_1	TRINITY_DN19459_c0_g1_i7	23.040	842	    522	    32	    42	    792	    2487	67	    2.12e-36	149
lcl|KJ658709.1_prot_AID23631.1_1	TRINITY_DN19459_c0_g1_i2	23.280	872	    530	    36	    20	    791	    2821	323	    5.55e-36	148
lcl|KJ658709.1_prot_AID23631.1_1	TRINITY_DN19459_c0_g1_i8	23.280	872	    530	    36	    20	    791	    2821	323	    5.68e-36	148
lcl|KJ658709.1_prot_AID23631.1_1	TRINITY_DN19459_c0_g1_i10	23.040	842	    522	    32	    42	    792	    2740	320	    8.36e-36	147
lcl|KJ658709.1_prot_AID23631.1_1	TRINITY_DN3954_c0_g1_i1	    23.826	852	    521	    38	    42	    805	    1542	3977	2.63e-35	146
lcl|KJ658709.1_prot_AID23631.1_1	TRINITY_DN19459_c0_g1_i11	23.245	869	    538	    35	    20	    792	    2574	67	    1.43e-34	143
lcl|KJ658709.1_prot_AID23631.1_1	TRINITY_DN3954_c0_g1_i3	    23.700	827	    507	    37	    42	    784	    1542	3902	8.63e-33	137
lcl|KJ658709.1_prot_AID23631.1_1	TRINITY_DN19459_c0_g1_i5	27.576	330	    190	    12	    504	    791	    1038	70	    2.84e-28	121
lcl|KJ658709.1_prot_AID23631.1_1	TRINITY_DN19459_c0_g1_i4	29.064	203	    134	    4	    598	    791	    675	    70	    2.94e-24	103
lcl|KJ658709.1_prot_AID23631.1_1	TRINITY_DN13661_c0_g2_i1	27.245	323	    206	    12	    504	    809	    977	    45	    4.59e-24	105
lcl|KJ658709.1_prot_AID23631.1_1	TRINITY_DN19459_c0_g1_i15	22.340	752	    454	    33	    20	    679	    2163	22	    2.11e-19	94.7
lcl|KJ658709.1_prot_AID23631.1_1	TRINITY_DN26450_c0_g1_i1	43.243	74	    42	    0	    269	    342	    245	    24	    3.96e-15	72.0
lcl|KJ658709.1_prot_AID23631.1_1	TRINITY_DN3954_c0_g1_i2	    23.775	408	    242	    20	    42	    404	    1542	2693	9.16e-05	47.4


Step by step: 
1) Find the annotated coding sequence on NCBI nucleotide
2) Download protein sequence
3) Rename the file and adjust the extension to .fasta (instead of .txt)
4) Upload .fasta file to the BLAST folder of the appropriate Trinity assembly directory
5) Run the blast command with your subject info changed (in the BLAST directory where the db is located!) 
6) Copy the info from the chosen BLAST hit into the spreadsheet
7) Copy the TRINITY ID into the command below:
   

Finding the transcript in your assembly:
(Do this in the directory your assembly is in) 
```
echo "TRINITY_DN1954_c0_g1_i4" > ids.txt
```
You can change this as needed, otherwise it will overwrite the ids.txt file every time you use this command. The text in " " is just an example here, which happens to be the top hit for piwi-1 (H. miamia). 

After you have a text file (which is unfortunately required for picking out using seqtk):
```
seqtk subseq Trinity-trim-ghost.fasta ids.txt > piwi1.fa
```


------------------------------
