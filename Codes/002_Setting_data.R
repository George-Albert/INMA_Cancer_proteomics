
###########################################################################
###########################################################################
###                                                                     ###
###        THIS SCRIPT SETS THE READS, METADATA AND FEATURE DATA        ###
###                                                                     ###
###########################################################################
###########################################################################

###########################
##  0. Load dependences  ##
###########################
{
  library(tidyverse)
  library(readxl)
  library(writexl)
  library(openxlsx)
  library(tools)
}

#################
##  Functions  ##
#################
{
  create_dir=function(x){suppressWarnings(dir.create(x,recursive=TRUE))}
  dcols=function(x){data.frame(colnames(x))}
  ul=function(x,n=5){x[1:min(nrow(x),n),1:min(ncol(x),n)]}
  options(width=1000)
}

###################################
##  Setting wd and loading data  ##
###################################
main_wd <- getwd()
# "INMA_cancer/Analysis"
setwd(main_wd)
input_dir   <- "Inputs"
output_dir  <- "Outputs"
raw_dir     <- "001_Raw_data"
quality_dir <- "002_Quality_annotated"
sample_names <- list.files(path = file.path(input_dir,quality_dir),pattern = "xlsx")

### Load each file and keep it in a list
list_of_files <- lapply(sample_names, function(x) read.xlsx(file.path(input_dir,quality_dir,x)))
length(list_of_files)

### Obtain the sample name
var_names_to_list <- file_path_sans_ext(basename(sample_names))
list_of_files <- setNames(list_of_files,var_names_to_list)

##############################
##    1_Au_4h_D240.xlsx     ##
##    10_Au_24h_D240.xlsx   ##
##    11_Au_24h_D240.xlsx   ##
##    12_Au_24h_D238.xlsx   ##
##   13_PEG_24h_D240.xlsx   ##
##   14_PEG_24h_D240.xlsx   ##
##   15_PEG_24h_D240.xlsx   ##
##   16_PEG_24h_D240.xlsx   ##
##  2_Au_4h_D240_bis.xlsx   ##
##    3_Au_4h_D240.xlsx     ##
##    4_Au_4h_D238.xlsx     ##
##  5_PEG_4h_D240_bis.xlsx  ##
##  6_PEG_4h_D240_bis.xlsx  ##
##    7_PEG_4h_D240.xlsx    ##
##    8_PEG_4h_D238.xlsx    ##
##    9_Au_24h_D240.xlsx    ##
##############################

paste("There are",length(sample_names),"files")

### Load each file separately 
# Iterate over each file and load it into the global environment
var_name_vec <- data.frame(Samples.Names=character())

for (file in sample_names) {
  # Read the Excel file
  df <- read.xlsx(file.path(input_dir,quality_dir,file))
  
  # Get the base name of the file without extension to use as variable name
  var_name <- file_path_sans_ext(basename(file))
  
  # Move leading numbers to the end of the variable name
  cleaned_var_name <- sub("^([0-9]+)_(.*)$", "\\2_\\1", var_name)
  var_name_vec[file,] <- cleaned_var_name
  # Assign the data frame to the global environment with the file name
  assign(cleaned_var_name, df, envir = .GlobalEnv)
}

rm(df)
paste("There are",length(ls(pattern = "Au"))+length(ls(pattern = "PEG")),"samples on the global environment")

### Save the new sample names
write.table(var_name_vec,file.path(input_dir,"txt","samples_names.txt"))

### Lets create de reads matrix
columns_to_keep=c("Gene.Name","accession","description","spec.count","EMPAI")

### Here we selected the columns of interest and add  the sample names to the columns
list_of_files <- setNames(list_of_files,var_name_vec$Samples.Names)

list_of_files_selected <- lapply(seq_along(list_of_files), function(i) {

  df <- list_of_files[[i]]
  df <- df[, columns_to_keep]
  colnames(df)[colnames(df) != "accession"] <- paste0(colnames(df)[colnames(df) != "accession"],"_",
                                                      names(list_of_files[i]))
  return(df)
})

### Now we want to combine all the rows by the accession column
full_data_df <- Reduce(function(x, y) merge(x, y, by = "accession", all = TRUE,sort=F), list_of_files_selected)
# full_data_df_1 <- Reduce(function(x, y) merge(x, y, by = "Gene.Name", all = TRUE,sort=F), list_of_files_selected)

# num_rows_list <- lapply(list_of_files_selected, nrow)
# total_rows <- sum(unlist(num_rows_list))
gene_name_columns   <- grep("Gene.Name_", colnames(full_data_df), value = TRUE)
description_columns <- grep("description_", colnames(full_data_df), value = TRUE)

# Make a copy of the original df to no modify it
full_data_df_copy <- full_data_df

# Create a new column "Gene.Name"combining the columns that contain 'Gene.Name'
full_data_df_copy$Gene.Name <- apply(full_data_df_copy[gene_name_columns], 1, function(x) {
  # Return the first no Nan value
  na.omit(x)[1]
})

# Save the original 'Gene.Name' columns and 'accession' in a separate data frame for future review
original_gene_name_data <- full_data_df_copy[, c("accession", gene_name_columns)]

# Remove the original 'Gene.Name' columns from the main data frame
data_df <- full_data_df_copy[, !(names(full_data_df_copy) %in% gene_name_columns)]

reads_spec <- data_df[,c(1,50,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48)]
reads_empai <- data_df[,c(1,50,4,7,10,13,16,19,22,25,28,31,34,37,40,43,46,49)]

colnames(reads_spec) <- gsub("spec.count_","",colnames(reads_spec))
colnames(reads_empai) <- gsub("EMPAI_","",colnames(reads_empai))

### Prepare the reads data
# Find indices of duplicates in Gene.Names

# duplicate_indices <- which(duplicated(reads$Gene.Name) | duplicated(reads$Gene.Name, fromLast = TRUE))
duplicate_indices <-which(duplicated(reads_spec$Gene.Name))

# Replace only the first occurrence of duplicates in Gene.Names with corresponding values in Protein.Names
reads_spec$Gene.Name[duplicate_indices] <- reads_spec$accession[duplicate_indices]
rownames(reads_spec[duplicate_indices,])
rownames(reads_spec) <- reads_spec$Gene.Name
length(which(reads_spec$accession != reads_empai$accession))
rownames(reads_empai) <- rownames(reads_spec)

# Remove the accession and the Gene.Name columns
reads_spec$accession <- NULL
reads_spec$Gene.Name <- NULL

reads_empai$accession <- NULL
reads_empai$Gene.Name <- NULL
# Remove Reverse proteins
reads_spec <- reads_spec %>% filter(!grepl("Reverse_", rownames(.))) 
reads_empai <- reads_empai[rownames(reads_spec),]

# rownames(reads) <- reads$Gene.Name
feature_data <- data_df[,c(1,50,2,5,8,11,14,17,20,23,26,29,32,35,38,41,44,47)]
feature_data$Gene.Name[duplicate_indices] <- feature_data$accession[duplicate_indices]
rownames(feature_data[duplicate_indices,])
rownames(feature_data) <- feature_data$Gene.Name
feature_data <- feature_data[rownames(reads_spec),]
# Create a new column "Description"combining the columns that contain 'Description'
feature_data$Description <- apply(feature_data[description_columns], 1, function(x) {
  # Return the first no Nan value
  na.omit(x)[1]
})

feature_data <- feature_data[,c(1,2,19)]

metadata <- data.frame(Original.ID = rownames(var_name_vec),Sample.Names = var_name_vec$Samples.Names)

# Separate the 'Sample.Name' column into four new columns
metadata <- separate(metadata, Sample.Names, into = c("Particle", "Time", "Dx", "Sample.Order"), sep = "_",remove = F)
index_bis <- which(metadata$Sample.Order == "bis")
# Replace de "bis" string by the  number of the sample
metadata[index_bis,]$Sample.Order <- c(2,5,6)
rownames(metadata) <- metadata$Sample.Names

### Create some metadata from description column of features df

#Protein name
pattern <- "^(.*?)\\sOS=.*$"
protein_name <- gsub(pattern, "\\1", feature_data$Description)
feature_data$Protein.Name <- protein_name

# OS
pattern <-"OS=(.*?)\\sOX="
organism_source <- str_match(feature_data$Description, pattern)[,2]
feature_data$Organism.Source <- organism_source

# Taxonomy
pattern <- "OX=(.*?)\\s(GN=|Pe=)"
taxonomy <- str_match(feature_data$Description, pattern)[,2]
feature_data$Taxonomy <- taxonomy

# Gene Name
pattern <- "GN=(.*?)\\sPe="
gene_name <- str_match(feature_data$Description, pattern)[,2]
feature_data$Gene.Name.description <- gene_name

# Pe (Protein Existance ???)
pattern <- "Pe=(.*?)\\sSV="
pe <- str_match(feature_data$Description, pattern)[,2]
feature_data$Pe <- pe
length(which(feature_data$Pe==3)) ### 24 proteins with low confidence levels

# SV (Sequence Version???)
pattern <- "SV=(\\d+)"
sv <- str_match(feature_data$Description, pattern)[,2]
feature_data$SV <- sv
length(which(feature_data$SV==4)) ### 20 proteins == 4, max is 5 

### Now we can remove the "Description" column
feature_data$Description <- NULL
rownames(feature_data) <- feature_data$Gene.Name

### Check congruency
length(which(rownames(feature_data) != rownames(reads_spec)))
length(which(reads_empai$Gene.Name != reads_spec$Gene.Name))
length(which(colnames(reads_empai) != rownames(metadata)))

metadata$Sample.Order <- as.numeric(metadata$Sample.Order)
metadata <- metadata[order(metadata$Sample.Order),]
metadata$short_setup <- paste0(metadata$Particle,"_",metadata$Time)

#Define colors and fill vectors
color_mapping <- c("Au_4h" = "green",
                   "PEG_4h" = "red",
                   "Au_24h" = "darkgreen",
                   "PEG_24h" = "firebrick4")
metadata$Fill  <- "white"
metadata$Color <- color_mapping[match(metadata$short_setup, unique(metadata$short_setup))]
metadata[9:16,]$Fill <- metadata[9:16,]$Color

reads_spec <-reads_spec[,rownames(metadata)] 
reads_empai <- reads_empai[,colnames(reads_spec)]

write.table(reads_empai,file= file.path(input_dir,"txt","reads_empai.txt"))
write.table(reads_spec,file= file.path(input_dir,"txt","reads_spec.txt"))
write.table(metadata,file= file.path(input_dir,"txt","metadata.txt"))
write.table(feature_data,file= file.path(input_dir,"txt","feature_data.txt"))


#################################################################
##            Plot the quantity of genes per sample            ##
#################################################################

reads <- reads_spec 

# Set all the NAn values to zero
reads[is.na(reads)] <- 0

# create a table to introduce to histogram
non_zero_counts <- reads %>%
  summarise(across(everything(), ~ sum(. > 0))) %>%
  pivot_longer(cols = everything(), names_to = "Condition", values_to = "NonZeroGeneCount")

# Create histogram
hist_plt <- ggplot(non_zero_counts, aes(x = Condition, y = NonZeroGeneCount)) +
  geom_col(fill = "skyblue",color="black") +
  labs(title = "Number of Genes with by Condition",
       x = "Condition",
       y = "Number of Genes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

hist_plt

create_dir(output_dir)

pdf(file = file.path(output_dir,"001_Num_genes_per_sample_all.pdf"))
print(hist_plt)
dev.off()

# Remove the "PEG_4h_D238_8" condition because of the low quanity of genes
reads[,"PEG_4h_D238_8"] <- NULL

# create a table to introduce to histogram
non_zero_counts <- reads %>%
  summarise(across(everything(), ~ sum(. > 0))) %>%
  pivot_longer(cols = everything(), names_to = "Condition", values_to = "NonZeroGeneCount")

# Create histogram of genes per sample before filtering
non_zero_counts$Condition <- factor(non_zero_counts$Condition,levels = non_zero_counts$Condition)
color <- c(rep("green",4),rep("red",3),rep("darkgreen",4),rep("firebrick4",4))

hist_plt_rm <- ggplot(non_zero_counts, aes(x = Condition, y = NonZeroGeneCount,
                                        fill=Condition)) +
  geom_col(color="black") +
  scale_fill_manual(values=color)+
  labs(title = "Number of Genes by Condition after removed the sample_8",
       x = "Condition",
       y = "Number of Genes") +
  ylim(0,500)+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

hist_plt_rm

pdf(file = file.path(output_dir,"002_Num_genes_after_sample8_removed.pdf"))
print(hist_plt_rm)
dev.off()

