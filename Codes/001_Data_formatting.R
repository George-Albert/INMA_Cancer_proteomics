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
