# Load packages
library(tidyverse)
library(ggplot2)

setwd("/Users/kkamundson/Documents/csu/projects/shale_biogeo_project/tool_creation/iteration_natcom")

#### READ IN FILES ####
# read in MAG table - this will always be the same (unless updated in the future by me)
mags <- read.delim("shale_MAGs_978_v2.txt", header = TRUE)

# read in user feature table 
feat <- read.delim("../iteration1/feature_table_w_tax.txt", header = TRUE, skip = 1)


################ FORMAT MAG TABLE ##############################################
# This is truncating the MAG tax string at different levels for a 'lookup' later on 
# Technically only needs to be done once per environment load 
# GENUS - truncate taxonomic string for genus level 
mags <- mags %>%
  mutate(MAG_genus = sub("(;s__.*)$", "", MAG_FULL_tax))

# FAMILY - truncate the taxonomy string 
mags <- mags %>%
  mutate(MAG_family = sub("(;g__.*)$", "", MAG_FULL_tax))

# ORDER - truncate the taxonomy string
mags <- mags %>%
  mutate(MAG_order = sub("(;f__.*)$", "", MAG_FULL_tax))

# CLASS - truncate the taxonomy string
mags <- mags %>%
  mutate(MAG_class = sub("(;o__.*)$", "", MAG_FULL_tax))

#write.csv(mags, "MAGs_trucated_tax_strings.csv")



################### FORMAT ASV TABLE ###########################################
# rename the first column (ASVs) as such 
colnames(feat)[1] <- "ASV"

# extract the columns 'ASV' and 'taxonomy' from the feature table & rename
taxa <- feat %>% select(ASV, taxonomy)
taxa <- as.data.frame(taxa)
colnames(taxa)[2] <- "ASV_FULL_tax"

taxa <- taxa %>%
  separate(ASV_FULL_tax, into = c("d","p","c","o","f","g","s"), sep = ";", remove = FALSE)

# GENUS - truncate the taxonomy string to everything after 'genus' by removing all species
taxa <- taxa %>%
  mutate(ASV_genus = sub("(;s__.*)$", "", ASV_FULL_tax))

# FAMILY - truncate the taxonomy string to family level
taxa <- taxa %>%
  mutate(ASV_family = sub("(;g__.*)$", "", ASV_FULL_tax))

# ORDER - truncate the taxonomy string to order level
taxa <- taxa %>%
  mutate(ASV_order = sub("(;f__.*)$", "", ASV_FULL_tax))

# CLASS - truncate the taxonomy string to class level
taxa <- taxa %>%
  mutate(ASV_class = sub("(;o__.*)$", "", ASV_FULL_tax))

# CHECK - count total number of rows in taxa and input 'feat' and checking that they match 
nrow(taxa)
nrow(feat)

# SAVE - count the total number of ASVs in a feature table
TOTAL_FEAT <- nrow(feat)

#write.csv(taxa, "taxa_out_test.csv")


################################################################################
######################## START OF COMPARING ####################################
################################################################################


########## COUNTING ASVs THAT ARE INSUFFICIENTLY CLASSIFIED ####################
# Counting the number of ASVs that are Unassigned at the domain (d) level
insuf_ASV_UNASSIGNED <- sum(taxa$d == "Unassigned", na.rm = TRUE)

# Counting the number of ASVs that do not classify past the domain level 
insuf_ASV_justDOMAIN <- sum(is.na(taxa$p))

# Counting the number of ASVs that do not classify past the phylum level 
insuf_ASV_justPHYLUM <- sum(is.na(taxa$c))

# Counting the number of ASVs that do not classify past the class level 
insuf_ASV_justCLASS <- sum(is.na(taxa$o))

# Make a new column 'match_level' and mark all of the above as having insufficient data to continue forward 
merged_data <- taxa %>%
  mutate(match_level = ifelse(is.na(o), "insufficient ASV classification", NA))

# Count the number of ASVs with insufficient classification to move forward
count_ASV_insuf <- sum(merged_data$match_level == "insufficient ASV classification", na.rm = TRUE)



#################### FULL TAX MATCH ############################################
# Assuming 'taxa' is the dataframe with the 'ASV_FULL_tax' column
# Assuming 'mags' is the dataframe with the 'MAG_FULL_tax', 'MAG', and 'compl' columns

# Select rows with the highest completion value for each 'MAG_FULL_tax'
# Need to do this as there might be more than one MAG with identical taxanomic classification - so in those cases it should choose the highest completion MAG
mags_filtered_FULL <- mags %>%
  group_by(MAG_FULL_tax) %>%
  slice_max(order_by = compl, with_ties = FALSE) %>%
  ungroup()

# Count the number of unique MAGs at full classification in dataset
MAGs_unique_FULL <- nrow(mags_filtered_FULL)

# Filter rows with NA in merged_data 
na_rows_full <- merged_data %>% filter(is.na(match_level))

# left join to look up MAG for matching full classificaiton between ASV and MAG
merged_data_full <- left_join(na_rows_full, mags_filtered_FULL %>% select(MAG_FULL_tax, MAG), by = c("ASV_FULL_tax" = "MAG_FULL_tax"))

# Add a new column based on condition
merged_data_full <- merged_data_full %>%
  mutate(match_level = ifelse(!is.na(MAG), "full tax", NA))

# Combine with original dataset
merged_data <- bind_rows(merged_data %>% filter(!is.na(match_level)), merged_data_full)

# CHECK SOME NUMBERS !
## First - merged_data should be the same as 'feat' 
nrow(merged_data)

# Second - how many ASVs match at the full taxanomic string match level?
count_ASV_FULL <- sum(merged_data$match_level == "full tax", na.rm = TRUE)



#################### GENUS LEVEL MATCH #########################################
# Assuming 'taxa' is the dataframe with the 'ASV_genus' column
# Assuming 'mags' is the dataframe with the 'MAG_genus' column
# Assuming 'merged_data' is the dataframe from the previous step

# Filter mags to just unique with highest completion 
mags_filtered_GENUS <- mags %>%
  group_by(MAG_genus) %>%
  slice_max(order_by = compl, with_ties = FALSE) %>%
  ungroup()

# Count the number of unique MAGs at genus level classification in dataset
MAGs_unique_GENUS <- nrow(mags_filtered_GENUS)

# Filter rows with NA in taxa match (MAG/level match columns)
na_rows_genus <- merged_data %>% filter(is.na(match_level))

# Left join on filtered rows
merged_data_genus <- left_join(na_rows_genus, mags_filtered_GENUS %>% select(MAG_genus, MAG), by = c("ASV_genus" = "MAG_genus"))

# Need to rename the column due to the default naming
merged_data_genus <- merged_data_genus %>%
  rename(MAG = MAG.y)

# Fill in new_column with 'genus' for newly added rows
merged_data_genus <- merged_data_genus %>%
  mutate(match_level = ifelse(!is.na(MAG), "genus", NA))

# Combine original merged_data with new_merged_data
merged_data <- bind_rows(merged_data %>% filter(!is.na(match_level)), merged_data_genus)

# Need to remove random MAG.x column it makes
merged_data <- merged_data %>%
  select(-MAG.x)

# Again - CHECK some numbers!
# check merged data
nrow(merged_data)

# Count how many ASVs matched at the genus level
count_ASV_GENUS <- sum(merged_data$match_level == "genus", na.rm = TRUE)




########################### FAMILY LEVEL #######################################
# Assuming 'taxa' is your dataframe with the 'ASV_genus' column
# Assuming 'mags' is your dataframe with the 'MAG_genus' column
# Assuming 'merged_data' is your dataframe from the previous step
# Filter mags to just unique with highest completion 

# Filter mags to just unique with highest completion at family level
mags_filtered_FAMILY <- mags %>%
  group_by(MAG_family) %>%
  slice_max(order_by = compl, with_ties = FALSE) %>%
  ungroup()

# Count the number of unique MAGs at family level classification in dataset
MAGs_unique_FAMILY <- nrow(mags_filtered_FAMILY)

# Filter rows with NA in taxa match (MAG/level match columns)
na_rows_family <- merged_data %>% filter(is.na(match_level))

# Left join on filtered rows
merged_data_family <- left_join(na_rows_family, mags_filtered_FAMILY %>% select(MAG_family, MAG), by = c("ASV_family" = "MAG_family"))

# Need to rename the column due to the default naming
merged_data_family <- merged_data_family %>%
  rename(MAG = MAG.y)

# Fill in new_column with 'family' for newly added rows
merged_data_family <- merged_data_family %>%
  mutate(match_level = ifelse(!is.na(MAG), "family", NA))

# Combine original merged_data with new_merged_data
merged_data <- bind_rows(merged_data %>% filter(!is.na(match_level)), merged_data_family)

# Need to remove random MAG.x column it makes
merged_data <- merged_data %>%
  select(-MAG.x)

# Again - CHECK some numbers!
# check merged data
nrow(merged_data)

# Count how many ASVs matched at the genus level
count_ASV_FAMILY <- sum(merged_data$match_level == "family", na.rm = TRUE)



######################### ORDER LEVEL ##########################################
# Assuming 'taxa' is the dataframe with the 'ASV_order' column
# Assuming 'mags' is the dataframe with the 'MAG_order' column
# Assuming 'merged_data' is the dataframe from the previous step
# Filter mags to just unique with highest completion 

# Filter mags to just unique with highest completion at order level
mags_filtered_ORDER <- mags %>%
  group_by(MAG_order) %>%
  slice_max(order_by = compl, with_ties = FALSE) %>%
  ungroup()

# Count the number of unique MAGs at order level classification in dataset
MAGs_unique_ORDER <- nrow(mags_filtered_ORDER)

# Filter rows with NA in taxa match (MAG/level match columns)
na_rows_order <- merged_data %>% filter(is.na(match_level))

# Left join on filtered rows
merged_data_order <- left_join(na_rows_order, mags_filtered_ORDER %>% select(MAG_order, MAG), by = c("ASV_order" = "MAG_order"))

# Need to rename the column due to the default naming
merged_data_order <- merged_data_order %>%
  rename(MAG = MAG.y)

# Fill in new_column with 'family' for newly added rows
merged_data_order <- merged_data_order %>%
  mutate(match_level = ifelse(!is.na(MAG), "order", NA))

# Combine original merged_data with new_merged_data
merged_data <- bind_rows(merged_data %>% filter(!is.na(match_level)), merged_data_order)

# Need to remove random MAG.x column it makes
merged_data <- merged_data %>%
  select(-MAG.x)

# Again - CHECK some numbers!
# check merged data
nrow(merged_data)

# Count how many ASVs matched at the genus level
count_ASV_ORDER <- sum(merged_data$match_level == "order", na.rm = TRUE)



################################################################################
################################################################################
################################################################################


####################### Filter out samples with low counts #####################
# Calculate the sums of each column from 2 to the second last column
feat_column_sums <- colSums(feat[2:(ncol(feat) - 1)])

# Identify columns that sum to less than 100
samples_to_remove <- names(feat_column_sums[feat_column_sums < 1000])

# Extract columns that sum to less than 100 into a new data matrix
removed_samples <- feat %>%
  select(ASV, all_of(samples_to_remove), taxonomy)

# Count the number of removed samples - subtracting by 2 since it also copies over the ASV and tax columns
removed_samples_numb <- (ncol(removed_samples)) -2

# Remove these columns from the original data matrix
feat_filt <- feat %>%
  select(-all_of(samples_to_remove))



################ Convert feature table to relative abundance ###################
# Function to convert to relative abundance 
convert_to_relative_abundance <- function(column) {
  column_sum <- sum(column)
  relative_abundance <- (column / column_sum) * 100
  return(relative_abundance)
}

# Use function on feature table to convert
feat_filt_relab <- feat_filt %>%
  mutate(across(2:(ncol(.) - 1), convert_to_relative_abundance))
  
#write.csv(feat_filt_relab, "relabund_test_out.csv")


######### Convert feature table to long format and bind with MAG data ##########
# Pivot longer to get the feature table in the right format
feat_filt_relab_long <- feat_filt_relab %>%
  pivot_longer(
    cols = 2:(ncol(feat_filt_relab) - 1),
    names_to = "Sample",
    values_to = "Relabund")

# Expand taxonomy column
feat_filt_relab_long <- feat_filt_relab_long %>%
  separate(taxonomy, into = c("d","p","c","o","f","g","s"), sep = ";", remove = FALSE)

# Match the ASV ID with the MAG it matched to from the comparing section 
feat_filt_relab_long <- feat_filt_relab_long %>%
  left_join(
    merged_data %>% select(ASV, MAG, match_level),
    by = "ASV")

# Match the MAG ID with the MAG metatdata in the 'mags' data matrix 
feat_filt_relab_long <- feat_filt_relab_long %>%
  left_join(
    mags %>% select(MAG, basin, compl, contam, sPROD, METHANO, rhodanase, phsA, dsrB),
    by = "MAG")


###############################################################################
###################### FINAL SAVED OUTPUT TABLES ###############################

# First need to collect the MAG data and add it to the ASV (merged_data) matrix
merged_data_FINAL <- merged_data %>%
  left_join(
    mags %>% select(MAG, basin, compl, contam, sPROD, METHANO, rhodanase, phsA, dsrB),
    by = "MAG")

# Need to establish an order to the data so that it presents in a logical order downstream
match_order <- c("full tax", "genus", "family", "order", "insufficient ASV classification", "NA")

# Clean up final merged data matrix for user output table
merged_data_FINAL <- merged_data_FINAL %>%
  arrange(match_level)

remove_columns <- c("d", "p", "c", "o", "f", "g", "s", "ASV_genus", "ASV_family", "ASV_order", "ASV_class")

merged_data_OUTPUT <- merged_data_FINAL %>%
  select(-all_of(remove_columns))

##### SAVE FINAL DATA FRAME FOR USER ####
write_csv(merged_data_OUTPUT, "ASV_to_MAGs_output.csv")


# Summarize the data to get counts of each match_level
match_level_counts <- merged_data_FINAL %>%
  group_by(match_level) %>%
  summarize(count = n()) %>%
  ungroup()

# Also want to add percentages to this - setting up for ggplot pie chart etc below
match_level_counts <- match_level_counts %>%
  mutate(percentage = count / sum(count) * 100) %>%
  mutate(percentage = round(count / sum(count) * 100, 2)) %>%
  mutate(label = paste0(match_level, "\n(n=", count, ", ", percentage, "%)"))

##### SAVE COUNTS AS OUTPUT ####
write.csv(match_level_counts, "ASV_matching_percents.csv", row.names = F)



################################################################################
################################################################################
############################ FIGURES TIME ! ####################################


##################### PIE CHART ################################################

# Need to order the data first so that it presents in a logical order
match_level_counts <- match_level_counts %>%
  mutate(match_level = factor(match_level, levels = match_order))

# Plot pie chart
p1 <- ggplot(match_level_counts, aes(x = "", y = count, fill = match_level)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c("#F2671F", "#C91B26", "#9C0F5F", "#60047A", "#160A47","#808080")) +
  coord_polar(theta = "y") +
  theme_void() +
  labs(title = "Matching of ASVs by taxonomic level",
       fill = "Match Level") +
  theme(legend.position = "right") +
  geom_text(aes(label = label), color="white", position = position_stack(vjust = 0.5))
p1
#ggsave("p1_PieChart.pdf", plot = p1)



###################### BUBBLE CHART ############################################
# For the next plot I want to show the lowest level of tax classification for an ASV
# Already have taxa string split into the various delimited levels
# So need to find the lowest level - aka where the first NA is

# Define the columns to check
find_lowest_tax <- c("d", "p", "c", "o", "f", "g", "s")

# Function to find the first NA column
first_na_column <- function(row) {
  first_na <- which(is.na(row))
  if (length(first_na) == 0) {
    return("s")  # No NAs in the row
  } else {
    return(names(row)[first_na[1]])
  }}

# Apply the function to each row in the data frame
feat_filt_relab_long$first_na <- apply(feat_filt_relab_long[find_lowest_tax], 1, first_na_column)

# Filter the table to remove really low abundance and NAs in the phylum column
feat_filt_relab_long_p2 <- feat_filt_relab_long %>%
  filter(Relabund >0.05) %>%
  filter(!is.na(p))

# Manipulate the file so that everything is ordered correctly for graphing
feat_filt_relab_long_p2 <- feat_filt_relab_long_p2 %>%
  mutate(first_na = factor(first_na, levels = c("c", "o", "f", "g", "s"))) %>%
  mutate(match_level = factor(match_level, levels = match_order))

# Make bubble plot
p2<-ggplot(feat_filt_relab_long_p2 %>% filter(Relabund >1), aes(x=first_na, y=p, fill = match_level)) +
  geom_jitter(aes(size=Relabund), shape = 21, width = 1, height = 0.5) +
  theme_bw() + 
  scale_fill_manual(values = c("#F2671F", "#C91B26", "#9C0F5F", "#60047A", "#160A47","#808080")) +
  labs(title = "Matching ASVs to MAGs\nASVs sized by relative abundance and colored by taxonomic level of match to MAG",
       fill = "Match level",
       size = "% Relative Abundance") +
  xlab("Lowest level of ASV taxonomic classification") + 
  ylab("Phylum") +
  guides(fill = guide_legend(override.aes = list(size = 5)))
p2  
#ggsave("p2_ASVs_that_match_1perc.pdf", plot = p2)




##################### STACKED BAR CHARTS COMM COMP #############################
# Establishing the color pallette becuase it became an issue in plotting
library(viridis)
plasma_palette <- viridis_pal(option = "plasma")(50)
reversed_palette <- rev(plasma_palette)

# Stack barcharts - ALL - ASV - PHYLUM 
p3 <- ggplot(feat_filt_relab_long, aes(x=Sample, y=Relabund, fill=rev(p))) +
  geom_bar(stat="identity") +
  theme_bw() + 
  labs(title = "Community composition per sample by ASV taxonomy",
       fill = "Phylum") +
  scale_fill_manual(values = reversed_palette) +
  guides(fill = guide_legend(ncol = 2)) +
  xlab("Sample") +
  ylab("% Relative Abundance") +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))
p3
#ggsave("p3_ASV_comcomp_Phylum.pdf", plot = p3)




############### Stack barcharts - MATCHING MAGS REPRESENTATIVES ################
# First need to filter and re-order some things so the plots look pretty
feat_filt_relab_long_ordered <- feat_filt_relab_long %>%
  mutate(match_level = factor(match_level, levels = rev(match_order)))

# Stack barcharts - ALL - ASV - FAMILY - MATCHING MAGS REPRESENTATIVES
p4<-ggplot(feat_filt_relab_long_ordered, aes(x=Sample, y=Relabund, fill=match_level)) +
  geom_bar(stat="identity") +
  theme_bw() + 
  labs(title = "Proportion of ASV community that linked to MAGs",
       fill = "Taxonomic level match\nof ASV-MAG match") +
  scale_fill_manual(values = c("#F2671F", "#C91B26", "#9C0F5F", "#60047A", "#160A47", "#808080"), limits = match_order) +
  xlab("Sample") +
  ylab("% Relative Abundance") +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))
p4
#ggsave("p4_ASV_comcomp_Matches.pdf", plot = p4)




############# BAR CHARTS LINKING TO FUNCTION ###################################
# Stack barcharts - FUNCTION - SULFIDE 
feat_filt_relab_long <- feat_filt_relab_long %>%
  rowwise() %>%
  mutate(S_pathway = case_when(
    sPROD == FALSE ~ "FALSE",
    is.na(sPROD) ~ NA_character_,
    sPROD == TRUE ~ {
      genes <- c()
      if (isTRUE(rhodanase)) genes <- c(genes, "rhodanase")
      if (isTRUE(dsrB)) genes <- c(genes, "dsrB")
      if (isTRUE(phsA)) genes <- c(genes, "phsA")
      if (length(genes) == 0) "TRUE (no genes)" else paste(genes, collapse = ", ")
    }
  )) %>%
  ungroup()

# Arrange positions of categories
feat_filt_relab_long <- feat_filt_relab_long %>%
  mutate(S_pathway = factor(S_pathway, 
                            levels = c(NA, "FALSE", 
                                       sort(unique(S_pathway[!is.na(S_pathway) & S_pathway != "FALSE"])))))

p5 <- ggplot(feat_filt_relab_long, aes(x=Sample, y=Relabund, fill=S_pathway)) +
  geom_bar(stat = "identity") + 
  theme_bw() +
  labs(title = "Proportion of inferred sulfide producers",
       fill = "Sulfide production\npathway") + 
  xlab("Sample") +
  ylab("% Relative Abundance") +
  scale_fill_manual(values = c("#d3d3d3", "#008080", "#b4c7a8", "#70a494", "#ca5828", "#f6ebbc", "#dd8a58")) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))
p5
ggsave("p5_ASV_function_sulfide_NEW.pdf", plot = p5)


# Stack barcharts - FUNCTION - METHANOGENS
p6<-ggplot(feat_filt_relab_long, aes(x=Sample, y=Relabund, fill=METHANO)) +
  geom_bar(stat = "identity") + 
  theme_bw() +
  labs(title = "Proportion of inferred methanogens",
       fill = "Potential\nmethanogen") +
  xlab("Sample") +
  ylab("% Relative Abundance") +
  scale_fill_manual(values = c("#d3d3d3", "#097969", "#434343")) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))
p6
ggsave("p6_ASV_function_methanogens.pdf", plot = p6)



