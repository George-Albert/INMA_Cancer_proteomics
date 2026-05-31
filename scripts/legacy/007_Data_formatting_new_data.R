###########################
##  0. Load dependences  ##
###########################
{
  library(tidyverse)
  library(readxl)
  library(writexl)
  library(openxlsx)
  library(stringr)
}

#################
##  Functions  ##
#################
{
  create_dir=function(x){suppressWarnings(dir.create(x,recursive=TRUE))}
  dcols=function(x){data.frame(colnames(x))}
  options(width=1000)
  # input_dir <- location of the input data
  # name      <-  name of the file to examine
  # Functio to clean the tables
  cleaning_tables <- function(name, input_dir = "Inputs/006_Raw_data_added"){
    
    tab <- read_excel(file.path(input_dir,name))
    tab <- data.frame(tab)
    
    # str(tab)
    # class(tab)# Inspect the object
    
    ### We identified numeric columns that are defined as characters and are the ones who have the semicolons
    ### We are going to remove the semicolons and change the data type of the columns with numbers from character 
    ### to numeric
    
    ### There are some rows with two values and also we have duplicate values (Check with María and JS)
    semicolon_raws <- grep("(.+);(.+)", tab$NSAF)
    semicolon_raws_empai <- grep("(.+);(.+)", tab$EMPAI)
    
    duplicates_in_raws <- tab[semicolon_raws,]
    duplicates_in_raws_empai <- tab[semicolon_raws_empai,]
    
    ### We removed by NSAF duplicates values
    ### Let`s remove it by now to follow with the cleaning
    tab[semicolon_raws,]$NSAF <- str_replace_all(tab[semicolon_raws,"NSAF"],";(.+)","")
    tab[semicolon_raws_empai,]$EMPAI <- str_replace_all(tab[semicolon_raws_empai,"EMPAI"],";(.+)","")
    
    # str(tab)
   safe<-tab
    for (col in seq_along(tab)){
      # if is character remove the semicolons
      if (is.character(tab[,col])) {
        
        tab[semicolon_raws,col] <- str_replace_all(tab[semicolon_raws,col],";(.+)","")
        tab[,col] <- str_replace_all(tab[,col],";","")
        tab[,col] <- str_replace_all(tab[,col],"E","e")
        tab[,col] <- str_replace_all(tab[,col],",",".")
        
      }
    }
    
    ### create a reg expression to remove the last decimal number from the row in order to get just 
    ### one and remove the duplicate
    
    index <- grep("^[^.]*\\.[^.]*\\.[^.]*$", tab$seq.coverage)
    pattern <-"\\d\\.\\d+$"
    tab[index,"seq.coverage"] <- gsub(pattern,"", tab[index,"seq.coverage"])
    
    ### Besides the semicolon, we found scientific notation in Excel that provokes NAs when 
    ### we convert to numeric
    safe <- tab # Create a copy
    
    ### Convert the columns seq.coverage, NSAF and EMPAI to numeric
    tab$seq.coverage <- as.numeric(tab$seq.coverage)
    tab$NSAF         <- as.numeric(tab$NSAF)
    tab$EMPAI        <- as.numeric(tab$EMPAI)
    
    ### Check the data type
    # str(tab)
    
    ### Nan values in Gene.Name column
    index_gene_name <- which(is.na(tab$Gene.Name))
    
    ### Most of these are contaminants ???
    nan_gene_name <- tab[index_gene_name,]
    
    ### Lets add the same name as the accession column to gene name
    tab[index_gene_name,"Gene.Name"] <- tab[index_gene_name,"accession"]
    
    ### Create a column with the gene name issues found
    tab$Gene.Name.issues <- "No"
    tab[index_gene_name,]$Gene.Name.issues <- "contaminant"
    
    ### Duplicates by Gene.Name
    duplicate_index <-which(duplicated(tab$Gene.Name)) 
    
    if (length(duplicate_index) > 0) {
      tab[duplicate_index,]$Gene.Name.issues <- "Duplicates"
    }
    
    ### Substitute the duplicates 
    ### Find duplicates in specific columns and add a column indicating duplicates
    tab <- tab %>%
      mutate(duplicated_row = duplicated(tab[, c("seq.coverage", "seq.count","spec.count","NSAF","EMPAI")]) | 
               duplicated(tab[, c("seq.coverage", "seq.count","spec.count","NSAF","EMPAI")], fromLast = TRUE))
    
    tab[which(tab$duplicated_row==T),]$duplicated_row <- "Duplicates"
    tab[which(tab$duplicated_row==F),]$duplicated_row <- "No"
    
    ### Duplicates by Accession
    duplicate_accession <- which(duplicated(tab$accession))
    # tab$Accesion.Dupes <- "No"  
    
    if (length(duplicate_accession) > 0) {
      # tab[duplicate_accession,]$Accesion.Dupes <- "Duplicates"
      print("There are accession Duplicates")
    } else{print("There are not accession Dupes")}
    
    {
    ### Create some metadata from description column
    #Protein name
    pattern <- "^(.*?)\\sOS=.*$"
    protein_name <- gsub(pattern, "\\1", tab$description)
    tab$Protein.Name <- protein_name
    
    # OS
    pattern <-"OS=(.*?)\\sOX="
    organism_source <- str_match(tab$description, pattern)[,2]
    tab$Organism.Source <- organism_source
    
    # Taxonomy
    pattern <- "OX=(.*?)\\s(GN=|Pe=)"
    taxonomy <- str_match(tab$description, pattern)[,2]
    tab$Taxonomy <- taxonomy
    
    # Gene Name
    pattern <- "GN=(.*?)\\sPe="
    gene_name <- str_match(tab$description, pattern)[,2]
    tab$Gene.Name.description <- gene_name
    
    # Pe (Protein Existance ???)
    pattern <- "Pe=(.*?)\\sSV="
    pe <- str_match(tab$description, pattern)[,2]
    tab$Pe <- pe
    length(which(tab$Pe==3)) ### 9 proteins with low confidence levels
    
    # SV (Sequence Version???)
    pattern <- "SV=(\\d+)"
    sv <- str_match(tab$description, pattern)[,2]
    tab$SV <- sv
    length(which(tab$SV==4)) ### 14 proteins == 4, max is 5 
    
    ### Remove the contaminants
    remove_index <- grep(pattern = "contaminant",tab$accession)  
    tab <- tab[-remove_index,]
    }
    
    # dim(tab)[1]
    
    # print(summary(tab))
    
    ### Save data already cleaned
    quality_dir <- "Inputs/007_Quality_annotated_added"
    create_dir(file.path(quality_dir))
    write.xlsx(tab, file = file.path(quality_dir,name))
    
    return(tab)
    
  }
  ### Function to verify if a dataframe is empty by rows
  is_not_empty <- function(df) {
    nrow(df) > 0
  }
}

###################################
##  Setting wd and loading data  ##
###################################
main_wd <- getwd()
setwd(main_wd)
input_dir <- "Inputs/006_Raw_data_added"
output_dir <- "Outputs"

### Rename the file names
sample_names <- list.files(path = input_dir,pattern = "xlsx")
sample_names_rename <- gsub(" ", "_",sample_names)
sample_names_rename[1]<- gsub("24_h", "24h",sample_names_rename[1])

file.rename(file.path(input_dir,sample_names),file.path(input_dir,sample_names_rename))

sample_names <- list.files(path = input_dir,pattern = "xlsx")
##################################
##  1. Annotate quality issues  ##
##################################

# Initialize empties data frame and lists to store data we want to save
dim_results_df      <- data.frame(Sample.Name = character(), N.genes = integer(), stringsAsFactors = FALSE)
summary_list        <- list()
dupe_gene_name_list <- list()
dupe_row_list       <- list()
### Apply workflow to all files
for(sample in sample_names){
  
  print(sample)
  tab=cleaning_tables(sample)
  
  # Store the dimensions of the first table in the data frame
  dim_results_df <- rbind(dim_results_df, data.frame(Sample.Name = sample, N.genes = dim(tab)[1]))
  
  # Extract the summary of each table
  summary_output <- summary(tab)
  summary_matrix <- matrix(summary_output, nrow=nrow(summary_output), ncol = ncol(summary_output))  
  colnames(summary_matrix) <- colnames(summary_output)
  summary_list[[sample]] <- as.data.frame(summary_matrix)
  
  # Extract the dupes by gene name
  dupe_gene_name <- tab[which(tab$Gene.Name.issues == "Duplicates"),]
  dupe_gene_name_list[[sample]] <- dupe_gene_name
  # Extract the dupes by "seq.coverage", "seq.count","spec.count","NSAF","EMPAI"
  dupe_row <- tab[which(tab$duplicated_row == "Duplicates"),]
  dupe_row_list[[sample]] <- dupe_gene_name
  
}

dupe_row_list_filtered <- Filter(is_not_empty,dupe_row_list)
dupe_gene_name_list_filtered <- Filter(is_not_empty,dupe_gene_name_list)

### Save data up to now
create_dir("Inputs/txt")
create_dir("Inputs/xlsx")

write.table(dim_results_df, file = file.path("Inputs","txt","Num_genes_per_condition_007_added.txt"))
write.xlsx(dim_results_df, file = file.path("Inputs","xlsx","Num_genes_per_condition_007_added.xlsx"))

# write.xlsx(summary_list, file = file.path("Inputs","xlsx","Summary_each_condition_tab.xlsx"))
# capture.output(summary_list, file = file.path("Inputs","txt","Summary_each_condition_tab.txt"))
write.xlsx(dupe_gene_name_list_filtered, file = file.path("Inputs","xlsx","dupe_gene_name_condition_007_added.xlsx"))
# capture.output(dupe_gene_name_list_filtered, file = file.path("Inputs","txt","dupe_gene_name_condition.txt"))
write.xlsx(dupe_row_list_filtered, file = file.path("Inputs","xlsx","dupe_rows_condition_007_added.xlsx"))

