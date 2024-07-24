
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
  library(plotrix)
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
# "INMA_cancer/Analyses"
setwd(main_wd)
input_dir   <- "Inputs"
output_dir  <- "Outputs"
raw_dir     <- "006_Raw_data_added"
quality_dir <- "007_Quality_annotated_added"
sample_names <- list.files(path = file.path(input_dir,quality_dir),pattern = "xlsx")

### Load each file and keep it in a list
list_of_files <- lapply(sample_names, function(x) read.xlsx(file.path(input_dir,quality_dir,x)))
length(list_of_files)

### Obtain the sample name
var_names_to_list <- file_path_sans_ext(basename(sample_names))
list_of_files <- setNames(list_of_files,var_names_to_list)

paste("There are",length(sample_names),"files")

### Lets create de reads matrix
columns_to_keep=c("Gene.Name","accession","description","spec.count","EMPAI")

### Here we selected the columns of interest and add  the sample names to the columns

list_of_files_selected <- lapply(seq_along(list_of_files), function(i) {

  df <- list_of_files[[i]]
  df <- df[, columns_to_keep]
  colnames(df)[colnames(df) != "accession"] <- paste0(colnames(df)[colnames(df) != "accession"],"_",
                                                      names(list_of_files[i]))
  return(df)
})

names(list_of_files_selected) <- names(list_of_files)

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

##############################
##  Prepare the Reads Data  ##
##############################

reads_spec <- data_df[,c(grep("Gene.Name|accession",colnames(data_df)),grep("spec.count",colnames(data_df)))]
reads_empai <- data_df[,c(grep("Gene.Name|accession",colnames(data_df)),grep("EMPAI",colnames(data_df)))]

colnames(reads_spec) <- gsub("spec.count_","",colnames(reads_spec))
colnames(reads_empai) <- gsub("EMPAI_","",colnames(reads_empai))
# Find indices of duplicates in Gene.Names

# duplicate_indices <- which(duplicated(reads$Gene.Name) | duplicated(reads$Gene.Name, fromLast = TRUE))
duplicate_indices_by_accession <-which(duplicated(reads_spec$accession))
reads_spec <- reads_spec[-duplicate_indices_by_accession,]
reads_empai <- reads_empai[-duplicate_indices_by_accession,]

duplicate_indices <-which(duplicated(reads_spec$Gene.Name))

# Replace only the first occurrence of duplicates in Gene.Names with corresponding values in Protein.Names
reads_spec$Gene.Name[duplicate_indices] <- reads_spec$accession[duplicate_indices]
reads_empai$Gene.Name[duplicate_indices] <- reads_empai$accession[duplicate_indices]

rownames(reads_spec) <- reads_spec$Gene.Name
rownames(reads_empai) <- reads_empai$Gene.Name

length(which(reads_spec$accession != reads_empai$accession))

# Remove the accession and the Gene.Name columns
reads_spec$accession <- NULL
reads_spec$Gene.Name <- NULL

reads_empai$accession <- NULL
reads_empai$Gene.Name <- NULL
# Remove Reverse proteins
reads_spec <- reads_spec %>% filter(!grepl("Reverse_", rownames(.))) 
reads_empai <- reads_empai[rownames(reads_spec),]

# rownames(reads) <- reads$Gene.Name

###########################
##  Create Feature Data  ##
###########################
feature_data <- data_df[,grep("accession|Gene.Name|description",colnames(data_df))]
feature_data <- feature_data[-duplicate_indices_by_accession,] 
feature_data$Gene.Name[duplicate_indices] <- feature_data$accession[duplicate_indices]
rownames(feature_data) <- feature_data$Gene.Name
feature_data <- feature_data[rownames(reads_spec),]

length(which(rownames(reads_spec)!=rownames(feature_data)))

# Create a new column "Description"combining the columns that contain 'Description'
feature_data$Description <- apply(feature_data[description_columns], 1, function(x) {
  # Return the first no Nan value
  na.omit(x)[1]
})

feature_data <- feature_data[,c("accession","Gene.Name","Description")]

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

#######################
##  Create Metadata  ##
#######################
var_name_vec <- var_names_to_list
metadata <- data.frame(Original.ID = var_name_vec)

# Fix nomenclature
metadata$Sample.Names <- gsub("^Au","Au_",metadata$Original.ID)
metadata$Sample.Names <- gsub("pd","Pd",metadata$Sample.Names )
metadata$Sample.Names <- gsub("5_min","5min",metadata$Sample.Names )

# Create the indices
evs_idx <- grep("^Evs_",metadata$Sample.Names)
cit_idx <- grep("Cit",metadata$Sample.Names)
lys_idx <- grep("Lys",metadata$Sample.Names)
au_idx  <- grep("Au",metadata$Sample.Names)
pd_idx  <- grep("Pd",metadata$Sample.Names)
pt_idx  <- grep("Pt",metadata$Sample.Names)
ctl_idx <- grep("Control",metadata$Sample.Names)
ctl_idx2 <- grep("5min",metadata$Sample.Names)
min_idx <- grep("5min",metadata$Sample.Names)
h4_idx <- grep("_4h",metadata$Sample.Names)
h24_idx <- grep("24h",metadata$Sample.Names)

extract_last_value <- function(x) {
  parts <- strsplit(x, "_")[[1]]
  return(tail(parts, n=1))
}

# Create columns
metadata$Evs     <- "NO"
metadata$Citrate <- "NO"
metadata$Lysine  <- "NO"
metadata$Au      <- "NO"
metadata$Pd      <- "NO"
metadata$Pt      <- "NO"
metadata$Control <- "NO"

# Evaluate the indices in the df
metadata$Evs[evs_idx]     <- "Yes"
metadata$Citrate[cit_idx] <- "Yes"
metadata$Lysine[lys_idx]  <- "Yes"
metadata$Au[au_idx]       <- "Yes"
metadata$Pd[pd_idx]       <- "Yes"
metadata$Pt[pt_idx]       <- "Yes"
metadata$Control[ctl_idx] <- "Yes"
metadata$Control[ctl_idx2]<- "Yes"
metadata$Sample.Order     <- sapply(metadata$Sample.Names,extract_last_value)
metadata$Time             <- "5min"
metadata$Time[min_idx]    <- "5min"
metadata$Time[h4_idx]     <- "4h"
metadata$Time[h24_idx]    <- "24h"

length(which(metadata$Time== "5min"))
# [1] 24
length(which(metadata$Time== "4h"))
# [1] 24
length(which(metadata$Time== "24h"))

metadata$Particle <- "No_particle"

metadata[metadata$Au=="Yes",]$Particle <- "Au"
metadata[metadata$Au=="Yes" & metadata$Lysine =="Yes",]$Particle <- "Au_Lys"
metadata[metadata$Au=="Yes" & metadata$Citrate =="Yes",]$Particle <- "Au_Cit"
metadata[metadata$Pd=="Yes",]$Particle <- "Pd"
metadata[metadata$Pt=="Yes",]$Particle <- "Pt"

metadata$short_setup <- paste0(metadata$Particle,"_Evs_",metadata$Evs,"_",metadata$Time)
metadata$Setup <- paste0(metadata$Particle,"_Evs_",metadata$Evs)

metadata <- metadata %>% 
  arrange(Particle, desc(Evs), desc(Time))

#Define colors and fill vectors

color_mapping <- c("#FFFF33", "gold", "goldenrod3", "grey","green","green", "blue", "blue")

names(color_mapping) <- unique(metadata$Setup)
col_vec <- sapply(color_mapping, function(x) color.id(x)[1])

metadata$Fill  <- "white"
metadata$Color <- col_vec[match(metadata$Setup, unique(metadata$Setup))]
metadata[metadata$Evs =="NO",]$Fill <- metadata[metadata$Evs =="NO",]$Color

# metadata_copy <- metadata



reads_spec <- reads_spec[,metadata$Original.ID]
colnames(reads_spec) <- metadata$Sample.Names

reads_empai <- reads_empai[,metadata$Original.ID]
colnames(reads_empai) <- metadata$Sample.Names

rownames(metadata) <- metadata$Sample.Names

### Check congruency
length(which(rownames(feature_data) != rownames(reads_spec)))
length(which(rownames(reads_empai) != rownames(reads_spec)))
length(which(colnames(reads_spec) != rownames(metadata)))

write.table(reads_empai,file= file.path(input_dir,"txt","reads_empai_added.txt"))
write.table(reads_spec,file= file.path(input_dir,"txt","reads_spec_added.txt"))
write.table(metadata,file= file.path(input_dir,"txt","metadata_added.txt"))
write.table(feature_data,file= file.path(input_dir,"txt","feature_data_added.txt"))


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

non_zero_counts$Color <- metadata$Color
color_base <- non_zero_counts$Color[which(!duplicated(non_zero_counts$Condition))]
names(color_base) <- non_zero_counts$Condition

non_zero_counts <- non_zero_counts %>% 
  mutate(Evs=ifelse(grepl("Evs",Condition),1,0)) %>% 
  arrange(desc(Evs))

non_zero_counts$Condition <- factor(non_zero_counts$Condition,levels = unique(non_zero_counts$Condition))

non_zero_counts_1 <- non_zero_counts[1:36,]
non_zero_counts_2 <- non_zero_counts[37:nrow(non_zero_counts),]

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 73.0   206.0   268.0   346.1   468.0  1668.0 

# Create histogram for all particles
hist_plt <- ggplot(non_zero_counts, aes(x = Condition, y = NonZeroGeneCount,fill = Condition)) +
  geom_col(color="black") +
  scale_fill_manual(values=color_base)+
  labs(title = "Number of Genes by Condition",
       x = "Condition",
       y = "Number of Genes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.text = element_text(size = 8),  # Ajusta el tamaño del texto de la leyenda
        legend.title = element_text(size = 10), # Ajusta el tamaño del título de la leyenda
        legend.key.size = unit(0.5, 'lines'))

hist_plt

create_dir(output_dir)

pdf(file = file.path(output_dir,"003_Num_genes_per_sample_all_added.pdf"),width = 15,height = 8)
print(hist_plt)
dev.off()

# Create histogram for Evs particles
non_zero_counts_1$Condition <- factor(non_zero_counts_1$Condition,levels = non_zero_counts_1$Condition)
color_base <- non_zero_counts_1$Color[which(!duplicated(non_zero_counts_1$Condition))]
names(color_base) <- non_zero_counts_1$Condition

hist_plt_1 <- ggplot(non_zero_counts_1, aes(x = Condition, y = NonZeroGeneCount,fill = Condition)) +
  geom_col(color="black") +
  scale_fill_manual(values=color_base)+
  labs(title = "Number of Genes by Condition (just Evs)",
       x = "Condition",
       y = "Number of Genes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.text = element_text(size = 8),  # Ajusta el tamaño del texto de la leyenda
        legend.title = element_text(size = 10), # Ajusta el tamaño del título de la leyenda
        legend.key.size = unit(0.5, 'lines'))

hist_plt_1

create_dir(output_dir)

pdf(file = file.path(output_dir,"004_Num_genes_per_sample_Evs_added.pdf"),width = 15,height = 8)
print(hist_plt_1)
dev.off()

# Create histogram for naked particles
non_zero_counts_2$Condition <- factor(non_zero_counts_2$Condition,levels = non_zero_counts_2$Condition)
color_base <- non_zero_counts_2$Color[which(!duplicated(non_zero_counts_2$Condition))]
names(color_base) <- non_zero_counts_2$Condition

hist_plt_2 <- ggplot(non_zero_counts_2, aes(x = Condition, y = NonZeroGeneCount,fill = Condition)) +
  geom_col(color="black") +
  scale_fill_manual(values=color_base)+
  labs(title = "Number of Genes by Condition (naked particles)",
       x = "Condition",
       y = "Number of Genes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.text = element_text(size = 8),  # Ajusta el tamaño del texto de la leyenda
        legend.title = element_text(size = 10), # Ajusta el tamaño del título de la leyenda
        legend.key.size = unit(0.5, 'lines'))

hist_plt_2

create_dir(output_dir)

pdf(file = file.path(output_dir,"005_Num_genes_per_sample_particle_naked_added.pdf"),width = 15,height = 8)
print(hist_plt_2)
dev.off()

