# Processing of metagenomic reads
### The following code is provided for the processing of raw metagenomic reads to MAGs for the 209 samples presented in Amundson et al., in prep.

## Raw reads to MAGs 
Raw metagenomic reads were first trimmed using sickle: 
```
sickle pe -f R1_All.fastq -r R2_All.fastq -t sanger -o R1_All_trimmed.fastq -p R2_All_trimmed.fastq -s R1R2_singles.fastq
```

Forward and reverse reads were merged, and assmebled with IDBA-UD:
```
fq2fa --merge --filter  R1_All_trimmed.fastq R2_All_trimmed.fastq R1R2_All_trimmed.fa
idba_ud -r R1R2_All_trimmed.fa  -o idba_assembled_output --num_threads 20
```

Another assembler, MEGAHIT, was used in select cases to assemble reads of samples that assembled poorly:
```
megahit --k-min 27 --k-max 127 --k-step 10 -1 R1_All_trimmed.fastq -2 R2_All_trimmed.fastq -m 0.9 -t 15 -o megahit_out 
```

Each assembly was next filtered to only contain contigs >5000bp. This subset assembly was used to map to determine contig coverage and bin metagenome-assembled genomes (MAGs). Coverage was determined by mapping reads to the fitlered assembly using bowtie: 
```
bowtie2-build scaffolds_5000.fa scaffolds_5000.fa_DB --threads 10

bowtie2 -p 15 -x scaffolds_5000.fa_DB -1 R1_All_trimmed.fastq -2 R2_All_trimmed.fastq --un unmapped_paired.fq --al mapped_paired.fa -S All_mappedtoall_paired.sam > bowtie_log
```

The sam file generated from bowtie2 was converted to a bam file and sorted using samtools: 
```
samtools view -@ 5 -bS All_mappedtoall_paired.sam > out.bam
samtools sort -@ 5 -T out.bam.sorted -o out.sorted.bam out.bam
```

Bins (MAGs) were identified using MetaBAT2 using the filtered contigs file (scaffold_5000.fa) and a depth file created by summarizing the sorted bam file using 'jgi_summarize_bam_contig_depths' which is built into the script 'runMetaBat.sh' provided by MetaBAT2. This was executed by: 
```
runMetaBat.sh --verysensitive PRB-71_scaffold_5000.fa PRB-71.sorted.bam
```

For some samples that assembled poorly, an iterative or sub-assembly approach was also performed. 
For sub assemblies, trimmed metagenomic reads were subsampled to 10 or 20%, merged, and then assembled with IDBA-UD
```
# Subsample R1 and R2 at 10%
python /opt/scripts/bin/pullseq_random_fastq.py -i R1_All_trimmed.fastq -o R1_All_trimmed_10percent.fastq -s 10
python /opt/scripts/bin/pullseq_random_fastq.py -i R2_All_trimmed.fastq -o R2_All_trimmed_10percent.fastq -s 10

# Merge the reads from the previous step
fq2fa --merge --filter R1_All_trimmed_10percent.fastq R2_All_trimmed_10percent.fastq R1R2_All_trimmed_10percent.fa

# Assemble subsetted metagenomic reads using IDBA-UD
idba_ud -r R1R2_All_trimmed_10percent.fa -o idba_assembled_output_10percent --num_threads 10
```

For iterative assemblies, trimmed metagenomic reads were first mapped to unique MAGs recovered for that basin. Then, reads that did not recruit were used for an iterative assembly. 
```
# Map to concatenated file of the already identified medium/high quality, unique bins from the same basin
bowtie2 -p 10 -x MON-2_3bins_toMapto.fa_DB -1 ../../../R1_All_trimmed.fastq -2 ../../../R2_All_trimmed.fastq --un-conc unmapped_paired.fq --al-conc mapped_paired.fa -S All_mappedtoall_paired.sam > bowtie_log

# Assemble with megahit (or with IDBA, as is provided above)
megahit --k-min 27 --k-max 127 --k-step 10 -1 unmapped_paired.1.fq -2 unmapped_paired.2.fq -m 0.9 -t 15 -o megahit_out --out-prefix MON-2_ir1
```


## Curating the MAG database 
To build a genome-resolved database from the 209 metagenomes, all assemblies (including sub-assemblies and iterative assemblies) were binned. MAGs were first evaluated for completeness and contamination using checkM:
```
checkm lineage_wf -t 5 -x fa scaffold_5000.fa.metabat-bins--verysensitive/ scaffold_5000.fa.metabat-bins--verysensitive/checkm

checkm qa -o 2 --tab_table -f scaffold_5000.fa.metabat-bins--verysensitive/checkm/checkm_summary.txt scaffold_5000.fa.metabat-bins--verysensitive/checkm/lineage.ms scaffold_5000.fa.metabat-bins--verysensitive/checkm
```

Only bins that achieved >50% completion and >10% contamination were retained. MAGs were first compiled at a basin-level and dereplicated at the strain (99% ANI) level. This resulted in 1100 MAGs unique to one basin from 11 distinct basins. To form the final database of MAGs, genomes were again dereplicated with dRep (v.3.0.0) as follows: 
```
dRep dereplicate all_PRB_bins/ -g all_PRB_bins/*.fa -comp 50 -p 10
```

This resulted in 978 unique MAGs from all metagenomes included in this study. 


## Taxonomic assignment and annotating MAGs for genomic traits and potential function
Next, taxonomy was assigned to the final MAG database and genomes were annotated for genomic traits and functional potential.

Taxonmy was assigned to MAGs using GTDB-tk (ref 220).
```
gtdbtk classify_wf -x fa --genome_dir ./ --out_dir gtdb_v2.4_ref220 --cpus 10 --mash_db /home/opt/gtdbtk/data/release220/mash
```

Genomic traits (optimal temperature, salinity, etc.) were predicted from genomeSPOT using MAG amino acid sequences, as identified by prodigal. 
```
# run prodigal on each MAG in a loop
prodigal -i ../"$element".fa -a "$element".faa -d "$element".fna -p meta -m -o "$element".outputdone

# run genomeSPOT on each MAG in a loop
python -m genome_spot.genome_spot --models /home/opt/Miniconda3/miniconda3/envs/GenomeSPOT/models --contigs ../prodigal/"$element".fna --proteins ../prodigal/"$element".faa --output "$element"
```

Data from genomeSPOT was produced for each individual genome as an output tsv file. Predictions for each MAG were combined into one single csv using the following code in R
```
library(dplyr)
library(tidyr)
library(readr)

# Read elements from the .txt file
# this is a list of MAGs, which is the prefix on the predictions.tsv files
elements <- read_lines("list_all_edit.txt")

# List all files
files <- list.files(pattern = "*.predictions.tsv")

# Ensure the length of elements matches the number of files
if (length(elements) != length(files)) {
  stop("The number of elements does not match the number of files.")
}

# Function to read and transform each file
read_and_transform <- function(file, element) {
  df <- read_tsv(file, col_types = cols(
    target = col_character(),
    value = col_character(),  # Read value as character to handle mixed types
    error = col_double(),
    units = col_character(),
    is_novel = col_logical(),
    warning = col_character()
  ))
  
  df <- df %>%
    mutate(element = element) %>%
    pivot_wider(names_from = target, values_from = c(value, error, units, is_novel, warning))
  return(df)
}

# Combine all files
combined_df <- mapply(read_and_transform, files, elements, SIMPLIFY = FALSE) %>%
  bind_rows()

# Convert numeric columns to appropriate types
numeric_cols <- names(combined_df)[grepl("value_|error_", names(combined_df))]
combined_df <- combined_df %>%
  mutate(across(all_of(numeric_cols), ~ as.numeric(.), .names = "converted_{.col}"))

# Write to a single file
write_csv(combined_df, "combined_data2.csv")
```

Each MAG was also annotated for functional potential using DRAM v1.4.4. 
```
source /opt/Miniconda2/miniconda2/bin/activate DRAM1.4.4

DRAM.py annotate -i './*.fa' -o ./ANNOTATIONS --threads 15 --use_fegenie --use_sulfur

DRAM.py distill -i ./ANNOTATIONS/annotations.tsv -o ./ANNOTATIONS/summarize --rrna_path ./ANNOTATIONS/rrnas.tsv --trna_path ./ANNOTATIONS/trnas.tsv
```

## Mapping reads to MAGs to determine relative abundance 
In order to determine relative abundance of each MAG, trimmed and rarefied metagenomic reads were mapped to the unique database of 978 medium and high quality MAGs. The sam files from mapping were converted and sorted using samtools, and were input to coverM to determine coverage of genomes across samples. To be considered present, a MAG had to achieve >1x coverage across 75% of its genome
```
# First, metagenomes >10Gbp were rarefied to 10Gbp to reduce sequencing bias
# $1 list of sample prefixes

for element in $(<$1)
do
gunzip "$element"_R1_All_trimmed.fastq.gz
gunzip "$element"_R2_All_trimmed.fastq.gz
reformat.sh in1="$element"_R1_All_trimmed.fastq in2="$element"_R2_All_trimmed.fastq out1="$element"_R1_All_trimmed_10Gbp.fastq out2="$element"_R2_All_trimmed_10Gbp.fastq -samplebasestarge t=10000000000 -sampleseed=1234
gzip "$element"_R1_All_trimmed.fastq
gzip "$element"_R2_All_trimmed.fastq
done
exit 0

# Next, bbmap was used to map reads against a concatenated file of 978 MAGs. Samtools was used to convert the sam file to a bam file, and sort it. The bbtools reformat script was used to ensure filtering of reads that achieved 95% identity of mapping. 
# $1 list of sample prefixes

for element in $(<$1)
do
gunzip ./"$element"_R1_All_trimmed_10Gbp.fastq.gz
gunzip ./"$element"_R2_All_trimmed_10Gbp.fastq.gz
bbmap.sh in1=./"$element"_R1_All_trimmed_10Gbp.fastq in2=./"$element"_R2_All_trimmed_10Gbp.fastq ambig=all ref=allShaleMAGs_978MQHQ_dRepd_12.16.2022.fa threads=15 out="$element"_mapped_rarified.sam -Xmx48G
samtools view -@ 15 -bS "$element"_mapped_rarified.sam > "$element"_mapped_rarified.bam
reformat.sh in="$element"_mapped_rarified.bam out="$element"_mapped_rarified_95.bam minidfilter=0.95 primaryonly=f
samtools sort -@ 15 -o "$element"_mapped_rarified_95_sorted.bam "$element"_mapped_rarified_95.bam
gzip ./"$element"_R1_All_trimmed_10Gbp.fastq
gzip ./"$element"_R2_All_trimmed_10Gbp.fastq
done
exit 0

# Finally, coverM was run using the bam files produced from mapping with bbmap in the previous step. Three commands were run.
# output reads_per_base
coverm genome --proper-pairs-only --genome-fasta-extension fa --genome-fasta-directory ./bins_978uniqueMAGs --bam-files *_95_sorted.bam --th
reads 20 --min-read-percent-identity-pair 0.95 --min-covered-fraction 0 -m reads_per_base --output-file coverm_reads_per_base.txt &> reads_per_base_stats.txt

# output MAGs with min-covered_fraction >0.75
coverm genome --proper-pairs-only --genome-fasta-extension fa --genome-fasta-directory ./bins_978uniqueMAGs --bam-files *_95_sorted.bam --th
reads 20 --min-read-percent-identity-pair 0.95 --min-covered-fraction 0.75 --output-file coverm_min75.txt &> min75_stats.txt

# ouput trimmed_mean
coverm genome --proper-pairs-only --genome-fasta-extension fa --genome-fasta-directory ./bins_978uniqueMAGs --bam-files *_95_sorted.bam --th
reads 20 --min-read-percent-identity-pair 0.95 -m trimmed_mean --output-file coverm_trimmed_mean.txt &> trimmed_mean_stats.txt
```

The three outputs from coverM were next combined based on the rules stated above using the following R code. Relative abundance of each MAG was determined as the proportion of coverage out of the total for a given sample.  
```
# Load the three matrix files - text files from coverm
file1 <- read.delim("coverm_reads_per_base.txt", header=TRUE, row.names = 1)
file2 <- read.delim("coverm_min75.txt", header=TRUE, row.names = 1)
file3 <- read.delim("coverm_trimmed_mean.txt", header=TRUE, row.names = 1)

# file one is reads per base, which needs to multiplied by 151
file1 <- file1 * 151

# convert .txt to .csv (tab to comma delimited)
#write.csv(file1, "file1.csv")
#write.csv(file2, "file2.csv")
#write.csv(file3, "file3.csv")

# now read in using read.csv
#file1 <- read.csv("file1.csv", header=TRUE, stringsAsFactors = FALSE)
#file2 <- read.csv("file2.csv", header=TRUE, stringsAsFactors = FALSE)
#file3 <- read.csv("file3.csv", header=TRUE, stringsAsFactors = FALSE)

# Check if the first column in row 1 of file2 reads 'unmapped'
# Remove the entire first row of file2
file2 <- file2[-1,]

# check the dimensions again -- will return nothing if they match
if (!all(dim(file1) == dim(file2) & dim(file2) == dim(file3))) {
  stop("The matrix files have different dimensions")
}

# force to recognize these files as matrixes
file1 <- as.matrix(file1)
file2 <- as.matrix(file2)
file3 <- as.matrix(file3)

# remove any NA, Na, or NaN values in the matrixes
file1[is.na(file1) | file1 == "Na" | is.nan(file1)] <- 0
file2[is.na(file2) | file2 == "Na" | is.nan(file2)] <- 0
file3[is.na(file3) | file3 == "Na" | is.nan(file3)] <- 0


# Create a new matrix with the same dimensions as the input matrixes
result <- matrix(0, nrow(file1), ncol(file1))

# Fill in the new matrix with values from file3 that meet the requirements in file1 and file2
# the requirement for file 1 (reads per base) is coverage, here 1x
# the requirement for file 2 (min75) is greater than zero
for (i in 1:nrow(file1)) {
  for (j in 1:ncol(file1)) {
    if (file1[i,j] > 1 & file2[i,j] > 0) {
      result[i,j] <- file3[i,j]
    }
  }
}

# Rename the columns of the resulting matrix to match the first file
colnames(result) <- colnames(file1)

# look at them to decide how to rename in the next step
print(colnames(result))

# Shorten the column names by removing the character string '_mapped_rarified_90_sorted.Reads.per.base'
colnames(result) <- gsub("_mapped_rarified_97_sorted.Reads.per.base", "", colnames(result))

# check again to make sure things look good 
print(colnames(result))

row.names(result) <- row.names(file1)

# Save the resulting matrix to a file
write.csv(result, "coverage_combined-coverM_1xcov-75genome.csv", row.names=TRUE)

# get column (sample) names and write them to a csv for making an NMDS
col_names <- colnames(result)
write.csv(col_names, "column_names.csv", row.names=FALSE)

# now convert to relativeabundance 
col_sums <-colSums(result)
relabund <- t((t(result) / col_sums) * 100)

write.csv(relabund, "relabund_combined-coverM_1xcov-75genome.csv", row.names = TRUE)
```





