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

Another assembler, MEGAHIT, was used to assemble reads of samples that assembled poorly:
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






