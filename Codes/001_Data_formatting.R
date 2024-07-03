<<<<<<< HEAD

#################################################################################
#################################################################################
###                                                                           ###
###  THIS SCRIPT CLEANS AND PREPARE THE DATA TO BE ABLE TO PROCESS AND CREATE ###
###               THE READS,METADATA AND FEATURE DATA MATRICES                ###
###                                                                           ###
#################################################################################
#################################################################################

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
  cleaning_tables <- function(name, input_dir = "Inputs/001_Raw_data" ){
    
    tab <- read_excel(file.path(input_dir,name))
    tab <- data.frame(tab)
    
    # str(tab)
    # class(tab)# Inspect the object
    
    ### We identified numeric columns that are defined as characters and are the ones who have the semicolons
    ### We are going to remove the semicolons and change the data type of the columns with numbers from character 
    ### to numeric
    safe <- tab # Create a copy
    
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
    
    for (col in seq_along(tab)){
      # if is character remove the semicolons
      if (is.character(tab[,col])) {
        
        tab[semicolon_raws,col] <- str_replace_all(tab[semicolon_raws,col],";(.+)","")
        tab[,col] <- str_replace_all(tab[,col],";","")
        tab[,col] <- str_replace_all(tab[,col],"E","e")
        
      }
    }
    
    ### Besides the semicolon, we found scientific notation in Excel that provokes NAs when 
    ### we convert to numeric
    
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
    
    # dim(tab)[1]
    
    # print(summary(tab))
    
    ### Save data already cleaned
    quality_dir <- "Inputs/002_Quality_annotated"
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
input_dir <- "Inputs/001_Raw_data"
output_dir <- "Outputs"
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

write.table(dim_results_df, file = file.path("Inputs","txt","Num_genes_per_condition.txt"))
write.xlsx(dim_results_df, file = file.path("Inputs","xlsx","Num_genes_per_condition.xlsx"))

write.xlsx(summary_list, file = file.path("Inputs","xlsx","Summary_each_condition_tab.xlsx"))
# capture.output(summary_list, file = file.path("Inputs","txt","Summary_each_condition_tab.txt"))
write.xlsx(dupe_gene_name_list_filtered, file = file.path("Inputs","xlsx","dupe_gene_name_condition.xlsx"))
# capture.output(dupe_gene_name_list_filtered, file = file.path("Inputs","txt","dupe_gene_name_condition.txt"))
write.xlsx(dupe_row_list_filtered, file = file.path("Inputs","xlsx","dupe_rows_condition.xlsx"))
=======
#This script creates a feature dictionary pairing RV_numbers, Gene_names and protein GI_number accessions.
# 0. Load dependences
# 1. ID conversion

# 0. Load dependences
{
    library(readxl)
    library(writexl)
    library(openxlsx)

    library(stringr)
    
    create_dir=function(x){suppressWarnings(dir.create(x,recursive=TRUE))}
    dcols=function(x){data.frame(colnames(x))}
    ul=function(x,n=5){x[1:min(nrow(x),n),1:min(ncol(x),n)]}
    options(width=1000)
}


## 1. Annotate quality issues.

{
    Annotate_quality_issues=function(name)
    {
        
        ## Cargar tabla
        tab=read_excel(paste0("Inputs/1_Raw/",name))
        tab=data.frame(tab)
        
        ## Ctr replace de los ;
        for(i in 1:ncol(tab))
        tab[,i]=str_replace_all(tab[,i],";","")
        
        ##Hay filas que parecen concatenar más de una proteína, de momento las marco como OK="No". En ellas las columnas sseq.coverage, NSAF y EMPAI e vuelven NAs al intentar pasarlas a numeros
        tab$seq.coverage=as.numeric(tab$seq.coverage)
        tab$seq.count=as.numeric(tab$seq.count)
        tab$spec.count=as.numeric(tab$spec.count)
        tab$NSAF=as.numeric(tab$NSAF)
        tab$EMPAI=as.numeric(tab$EMPAI)
        tab$HeavySpec=as.numeric(tab$HeavySpec)
        tab$LightSpec=as.numeric(tab$LightSpec)

        
        set_coverage=(which(is.na(tab$seq.coverage)))
        set_NSAF=(which(is.na(tab$NSAF)))
        set_EMPAI=(which(is.na(tab$EMPAI)))
        set_any_issue=(which(is.na(tab$seq.coverage) | is.na(tab$NSAF) | is.na(tab$EMPAI)))
        #print(paste(length(set_coverage),"NAs en seqcover",length(set_NSAF),"n NSAF y",length(set_EMPAI),"en EMPAI:",length(set_consensus), "en todos."))
        
        ## Por otra parte, en algunas está el accession pero no el gene.name. En ellas copio el accession a gene name de momento.
        set_missing_gene_names=which(is.na(tab$Gene.Name))
        tab$Gene.Name[set_missing_gene_names]=tab$accession[set_missing_gene_names]
        genes_dupe=tab$Gene.Name[which(duplicated(tab$Gene.Name))]
        set_genes_dupes_all=which(tab$Gene.Name %in% genes_dupe)
        set_genes_dupes_all_but_first=which(duplicated(tab$Gene.Name))
        
        #print(paste(length(set_missing_gene_names),"NAs en gene_name."))
        
        
        ## Creo una columna llamada "Ok" para definir tres status de QC: OK="Yes", "No" (filas múltiples) y "Find_gene_name".
        tab$Ok="Yes"
        tab$Ok[set_genes_dupes_all_but_first]="Dupe_gene_names"
        tab$Ok[set_missing_gene_names]="Find_gene_name"
        tab$Ok[set_any_issue]="No"
        
        tab$Ok=factor(tab$Ok)
        
        ## Imprimo el summary para ver que todo ha ido bien
        print(summary(tab))

        ## Vuelvo a repetir la carga y limpiado de ; para que lo que escriba esté igual que la original con la columna de Ok
        safe=tab
        tab=read_excel(paste0("Inputs/1_Raw/",name))
        tab=data.frame(tab)
        
        ## Ctr replace de los ;
        for(i in 1:ncol(tab))
        tab[,i]=str_replace_all(tab[,i],";","")
        tab$Ok=safe$Ok

        
        
        ## Salvo la tabla.
        create_dir("Inputs/2_Quality_annotated")
        write.xlsx(tab, file = paste0("Inputs/2_Quality_annotated/",name))
        return(tab)
    }
    
    sample_names <- list.files(path = "Inputs/1_Raw")[1:12]
    
    for(sample in sample_names)
        tab=Annotate_quality_issues(sample)
}

## 2. Load Quality-annotated sample tables, and select the good quality, interesting parts:
## row-wise: select only rows that are Ok="Yes" or "Find_Gene_name" (discard dupes and multiples for now), and set accession as gene_names of those missing.
## columns-wise: select only the columns that are interesting (in principle)

{
    filter_useful_data=function(name,columns_to_keep=c("Gene.Name","accession","description","spec.count"))
    {
        
        ## Cargar tabla
        tab=read_excel(paste0("Inputs/2_Quality_annotated/",name))
        tab=data.frame(tab)
        tab=tab[which(tab$Ok %in% c("Yes","Find_gene_name")),]
      
        tab$seq.coverage=as.numeric(tab$seq.coverage)
        tab$seq.count=as.numeric(tab$seq.count)
        tab$spec.count=as.numeric(tab$spec.count)
        tab$NSAF=as.numeric(tab$NSAF)
        tab$EMPAI=as.numeric(tab$EMPAI)
        tab$HeavySpec=as.numeric(tab$HeavySpec)
        tab$LightSpec=as.numeric(tab$LightSpec)
        
        set_missing_gene_names=which(is.na(tab$Gene.Name))
        #print(paste(length(set_missing_gene_names),"NAs en gene_name."))
        
        tab$Gene.Name[set_missing_gene_names]=tab$accession[set_missing_gene_names]
                
        tab$Ok=factor(tab$Ok)
        
        ## Imprimo el summary para ver que todo ha ido bien
        print(summary(tab))
        
        return(tab[,columns_to_keep])
    }
    
    ## Esto lo muevo a mano:
    Au_4h_1=filter_useful_data("1_Au_4h_D240.xlsx")
    Au_4h_2=filter_useful_data("2_Au_4h_D240_bis.xlsx")
    Au_4h_3=filter_useful_data("3_Au_4h_D240.xlsx")
    PEG_4h_4=filter_useful_data("4_PEG_4h_D240.xlsx")
    PEG_4h_5=filter_useful_data("5_PEG_4h_D240_bis.xlsx")
    PEG_4h_6=filter_useful_data("6_PEG_4h_D240_bis.xlsx")
    Au_24h_7=filter_useful_data("7_Au_24h_D240.xlsx")
    Au_24h_8=filter_useful_data("8_Au_24h_D240.xlsx")
    Au_24h_9=filter_useful_data("9_Au_24h_D240.xlsx")
    PEG_4h_10=filter_useful_data("10_PEG_24h_D240.xlsx")
    PEG_4h_11=filter_useful_data("11_PEG_24h_D240.xlsx")
    PEG_4h_12=filter_useful_data("12_PEG_24h_D240.xlsx")
    
}

## 3. Merge carefully the matrices to get the data matrix, along with the feature table and the metadata table; and save.

{
    
    careful_merge_pairs=function(a,b,name_a,name_b,target_column){
                
        all=merge(a,b,by="Gene.Name",all=TRUE)
        
        
        ## Copio los accessions y descriptions de las proteinas solo presentes en una tabla, y les pongo las spec counts a cero.
        all$accession.x[which(is.na(all$accession.x))]=all$accession.y[which(is.na(all$accession.x))]
        all$accession.y[which(is.na(all$accession.y))]=all$accession.x[which(is.na(all$accession.y))]

        all$description.x[which(is.na(all$description.x))]=all$description.y[which(is.na(all$description.x))]
        all$description.y[which(is.na(all$description.y))]=all$description.x[which(is.na(all$description.y))]
        
        tgt_x=paste0(target_column,".x")
        tgt_y=paste0(target_column,".y")

        all[,tgt_x][which(is.na(all[,tgt_x]))]=0
        all[,tgt_y][which(is.na(all[,tgt_y]))]=0

        ## Miro si hay discrepancias en accessions y/o descripciones linkadas a un solo gene name
        issues_accession=length(which(all$accession.x!=all$accession.y))
        issues_description=length(which(all$description.x!=all$description.y))
        set_issues=which(all$accession.x!=all$accession.y | (all$description.x!=all$description.y))
        names_x=paste0(colnames(a),".x")
        names_tgt=c(paste0(colnames(a),".x"),paste0(target_column,".y"))
        names_tgt[1]="Gene.Name"
        if(length(set_issues)>0)
        {
            print("Incongruencias al mergear, me quedo con las anotaciones .x")
            print(all[set_issues,])
        }else{print("Todo ok")}
        all=all[,names_tgt]
        colnames(all)[2]="accession"
        colnames(all)[3]="description"
        colnames(all)[4]=name_a
        colnames(all)[5]=name_b
        
        return(all)

    }
    
    Samples_1_2=careful_merge_pairs(Au_4h_1,Au_4h_2,name_a="Au_4h_1",name_b="Au_4h_2",target_column="spec.count")
    Samples_3_4=careful_merge_pairs(Au_4h_3,PEG_4h_4,name_a="Au_4h_3",name_b="PEG_4h_4",target_column="spec.count")
    Samples_5_6=careful_merge_pairs(PEG_4h_5,PEG_4h_6,name_a="PEG_4h_5",name_b="PEG_4h_6",target_column="spec.count")
    Samples_7_8=careful_merge_pairs(Au_24h_7,Au_24h_8,name_a="Au_24h_7",name_b="Au_24h_8",target_column="spec.count")
    Samples_9_10=careful_merge_pairs(Au_24h_9,PEG_4h_10,name_a="Au_24h_9",name_b="PEG_4h_10",target_column="spec.count")
    Samples_11_12=careful_merge_pairs(PEG_4h_11,PEG_4h_12,name_a="PEG_4h_11",name_b="PEG_4h_12",target_column="spec.count")
    #"Todo ok" en todas

   

   careful_merge_groups=function(a,b){
               name_a=colnames(a)[4:ncol(a)]
               name_b=colnames(b)[4:ncol(b)]

       all=merge(a,b,by="Gene.Name",all=TRUE)
       
       ## Copio los accessions y descriptions de las proteinas solo presentes en una tabla, y les pongo las spec counts a cero.
       all$accession.x[which(is.na(all$accession.x))]=all$accession.y[which(is.na(all$accession.x))]
       all$accession.y[which(is.na(all$accession.y))]=all$accession.x[which(is.na(all$accession.y))]

       all$description.x[which(is.na(all$description.x))]=all$description.y[which(is.na(all$description.x))]
       all$description.y[which(is.na(all$description.y))]=all$description.x[which(is.na(all$description.y))]
       
       for(name in c(name_a,name_b))
       all[,name][which(is.na(all[,name]))]=0

       ## Miro si hay discrepancias en accessions y/o descripciones linkadas a un solo gene name
       issues_accession=length(which(all$accession.x!=all$accession.y))
       issues_description=length(which(all$description.x!=all$description.y))
       set_issues=which(all$accession.x!=all$accession.y | (all$description.x!=all$description.y))
       names_tgt=c("Gene.Name","accession.x","description.x",name_a,name_b)
       if(length(set_issues)>0)
       {
           print("Incongruencias al mergear, me quedo con las anotaciones .x")
           print(all[set_issues,])
       }else{print("Todo ok")}
       all=all[,names_tgt]
       colnames(all)[2]="accession"
       colnames(all)[3]="description"
       return(all)

   }
   
   samples_1_to_4=careful_merge_groups(Samples_1_2,Samples_3_4)
   #"Todo ok"
   samples_5_to_8=careful_merge_groups(Samples_5_6,Samples_7_8)
   #"Todo ok"
   samples_9_to_12=careful_merge_groups(Samples_9_10,Samples_11_12)
   
   #"Incongruencias al mergear, me quedo con las anotaciones .x"
   #  Gene.Name                   accession.x                                                                                                    description.x Au_24h_9 PEG_4h_10 accession.y                                                                                                    description.y PEG_4h_11 PEG_4h_12
#   3     ACAP1 Reverse_sp|A5PK26|ACAP1_BOVIN Arf-GAP with coiled-coil, ANK repeat and PH domain-containing protein 1 OS=Bos taurus OX=9913 GN=ACAP1 PE=2 SV=1        0         2      A5PK26 Arf-GAP with coiled-coil, ANK repeat and PH domain-containing protein 1 OS=Bos taurus OX=9913 GN=ACAP1 PE=2 SV=1         2         0


   samples_1_to_8=careful_merge_groups(samples_1_to_4,samples_5_to_8)

   samples_1_to_12=careful_merge_groups(samples_1_to_8,samples_9_to_12)


   check_uniqueness_genes=function(tab){
       genes_dupes=tab$Gene.Name[which(duplicated(tab$Gene.Name))]
       print(tab[which(tab$Gene.Name %in% genes_dupes),])
   }

   check_uniqueness_genes(Samples_1_2)
   check_uniqueness_genes(Samples_3_4)
   check_uniqueness_genes(Samples_5_6)
   check_uniqueness_genes(Samples_7_8)
   check_uniqueness_genes(Samples_9_10)
   check_uniqueness_genes(Samples_11_12)
   check_uniqueness_genes(samples_1_to_4)
   check_uniqueness_genes(samples_5_to_8)
   check_uniqueness_genes(samples_9_to_12)
   check_uniqueness_genes(samples_1_to_8)
   check_uniqueness_genes(samples_1_to_12)

   reads=samples_1_to_12[,4:15]
   rownames(reads)=samples_1_to_12$Gene.Name
   feature_data=samples_1_to_12[,2:3]
   rownames(feature_data)=samples_1_to_12$Gene.Name

   meta_data=data.frame(Particle=rep(c(rep("Au",3),rep("PEG",3)),2),Time=c(rep("4h",6),rep("24h",6)))
   rownames(meta_data)=colnames(reads)
   
   create_dir("Inputs/3_Merged")
   write.table(reads, "Inputs/3_Merged/spectral_reads.txt")
   write.table(feature_data, "Inputs/3_Merged/feature_data.txt")
   write.table(meta_data, "Inputs/3_Merged/meta_data.txt")

   }
>>>>>>> 3367b67161cf93f2b757dcc1d606f2ece18fadcd
