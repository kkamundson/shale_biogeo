# Processing of metagenomic reads
### The following code is provided for the processing of raw metagenomic reads to MAGs for the 209 samples presented in Amundson et al., in prep.

## Raw reads to MAGs 
Raw metagenomic reads were first trimmed using sickle: 
'''
sickle pe -f R1_All.fastq -r R2_All.fastq -t sanger -o R1_All_trimmed.fastq -p R2_All_trimmed.fastq -s R1R2_singles.fastq
'''

Forward and reverse reads were merged, and assmebled with IDBA-UD:
'''
fq2fa --merge --filter  R1_All_trimmed.fastq R2_All_trimmed.fastq R1R2_All_trimmed.fa
idba_ud -r R1R2_All_trimmed.fa  -o idba_assembled_output --num_threads 20
'''

Another assembler, MEGAHIT, was used to assemble reads of samples that assembled poorly:
'''
megahit --k-min 27 --k-max 127 --k-step 10 -1 R1_All_trimmed.fastq -2 R2_All_trimmed.fastq -m 0.9 -t 15 -o megahit_out 
'''

Each assembly was next filtered to only contain contigs >5000bp. This subset assembly was used to map to determine contig coverage and bin metagenome-assembled genomes (MAGs). Coverage was determined by mapping reads to the fitlered assembly using bowtie: 
'''
bowtie2-build scaffolds_5000.fa scaffolds_5000.fa_DB --threads 10

bowtie2 -p 15 -x scaffolds_5000.fa_DB -1 R1_All_trimmed.fastq -2 R2_All_trimmed.fastq --un unmapped_paired.fq --al mapped_paired.fa -S All_mappedtoall_paired.sam > bowtie_log
'''

The sam file generated from bowtie2 was converted to a bam file and sorted using samtools: 
'''
samtools view -@ 5 -bS All_mappedtoall_paired.sam > out.bam
samtools sort -@ 5 -T out.bam.sorted -o out.sorted.bam out.bam
'''

Bins (MAGs) were identified using MetaBAT2 using the filtered contigs file (scaffold_5000.fa) and a depth file created by summarizing the sorted bam file using 'jgi_summarize_bam_contig_depths' which is built into the script 'runMetaBat.sh' provided by MetaBAT2. This was executed by: 
'''
runMetaBat.sh --verysensitive PRB-71_scaffold_5000.fa PRB-71.sorted.bam
'''








