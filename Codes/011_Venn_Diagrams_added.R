############################################################################
############################################################################
###                                                                      ###
###                       EXPRESSION LEVELS (Venn Diagram) SCRIPT        ###
###                                                                      ###
############################################################################
############################################################################

###########################
##  0. Load dependences  ##
###########################
{
  library(readxl)
  library(tidyverse)
  library(writexl)
  library(openxlsx)
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(stringr)
  library(cowplot)
  library(openxlsx)
  library(preprocessCore)
  library(stringr)
  library(reshape2)
  library(eulerr)
}

#################
##  Functions  ##
################# 
{
  create_dir=function(x){suppressWarnings(dir.create(x,recursive=TRUE))}
  dcols=function(x){data.frame(colnames(x))}
  ul=function(x,n=5){x[1:min(nrow(x),n),1:min(ncol(x),n)]}
  options(width=1000)
  my_name <- function(v1) {
    deparse(substitute(v1))
  }
}

###################################
##  Setting wd and loading data  ##
###################################

main_wd <- getwd()
# "/Users/George_A/Documents/Portafolio/INMA_cancer/Analyses"
setwd(main_wd)
input_dir <-"Inputs" 
output_dir <- "Outputs"

feature_data=read.table(file.path(input_dir,"txt","feature_data_added.txt"))
metadata=read.table(file.path(input_dir,"txt","metadata_filtered_added.txt"))
reads_spec=read.table(file.path(input_dir,"txt","reads_spec_added.txt"))

reads_spec <- reads_spec[,rownames(metadata)]
reads_spec[is.na(reads_spec)] <- 0
length(which(rownames(metadata) != colnames(reads_spec)))

##################################################################
##                     SET INPUT PARAMETERS                     ##
##################################################################
reads <- reads_spec

remove_extreme_vec <- c("variance","mean","both")
rm_ext <- remove_extreme_vec[3]

th_mean_vec <- c(0.2,0.5,1)
th_mean <- th_mean_vec[1]

#########################
##  Expression Levels  ##
#########################

# Log Transformation
exp=log2(reads+1)

exp_norm=normalize.quantiles.robust(as.matrix(exp),copy=FALSE, 
                                    remove.extreme=rm_ext,
                                    n.remove=1,use.median=FALSE,
                                    use.log2=FALSE)

### number of samples and genes
n_samples <- ncol(reads)
n_genes <- nrow (reads)
th <- 0.05

### Design matrix 
metadata$short_setup <- factor(metadata$short_setup,levels = unique(metadata$short_setup))
metadata$Time=factor(metadata$Time,levels=c("5min","4h","24h"))

group <- factor(metadata$short_setup)

design=model.matrix(~0+group,data=metadata$short_setup)
colnames(design) <- levels(group)

fit=lmFit(exp_norm,design)
fit2=eBayes(fit,trend=T, robust=T)

betas <- data.frame(fit2$coefficients)
p.value <- fit2$p.value

conditions <- colnames(betas)

# Set the threshold 
umbrales=quantile(exp_norm,c(0.75,0.8,0.95))
print(umbrales)
vec <- umbrales[3]
venn_list <- list()

for (condition in conditions) {
  
  
  print(condition)

  exp_levels   <- data.frame(row.names=rownames(betas),betas[condition])
  exp_genes=rownames(exp_levels)[which(exp_levels>vec)]
  
  venn_dir <- "010_Venn_diagram_added"
  create_dir(file.path(input_dir,venn_dir))

  write.table(exp_genes, file.path(input_dir,venn_dir,paste0(condition,"_quant_",names(vec),".txt")))
  
  create_dir(file.path(output_dir,"Figures",venn_dir))
  
  venn_list[[condition]] <- exp_genes

  name_euler <- paste0("Euler_diagram_quant_",names(vec),".pdf")
  name_euler <- gsub("%","_percent",name_euler)
}
  
  combinations <- list(venn_list[c(1,4,7)],venn_list[c(2,5,8)],venn_list[c(3,6,9)],
                    venn_list[c(13,16)],venn_list[c(14,17)],venn_list[c(15,18)],
                    venn_list[c(19,22)],venn_list[c(20,23)],venn_list[c(21,24)],
                    venn_list[c(1,10)],venn_list[c(2,11)],venn_list[c(3,12)],
                    venn_list[c(13,10)],venn_list[c(14,11)],venn_list[c(15,12)],
                    venn_list[c(19,10)],venn_list[c(20,11)],venn_list[c(21,12)],
                    venn_list[c(1,2,3)],venn_list[c(4,5,6)],venn_list[c(7,8,9)],
                    venn_list[c(16,17,18)],venn_list[c(22,23,24)],
                    venn_list[c(13,14,15)],venn_list[c(19,20,21)]
                    )
  names_combinations <- c("Au_naked_vs_Au_Evs_5_min","Au_naked_vs_Au_Evs_4_h","Au_naked_vs_Au_Evs_24_h",
                          "Pd_naked_vs_Pd_Evs_5_min", "Pd_naked_vs_Pd_Evs_4_h", "Pd_naked_vs_Pd_Evs_24_h",
                          "Pt_naked_vs_Pt_Evs_5_min", "Pt_naked_vs_Pt_Evs_4_h", "Pt_naked_vs_Pt_Evs_24_h",
                          "Au_Evs_vs_No_particle_Evs_5_min", "Au_Evs_vs_No_particle_Evs_4_h", "Au_Evs_vs_No_particle_Evs_24_h",
                          "Pd_Evs_vs_No_particle_Evs_5_min", "Pd_Evs_vs_No_particle_Evs_4_h", "Pd_Evs_vs_No_particle_Evs_24_h",
                          "Pt_Evs_vs_No_particle_Evs_5_min", "Pt_Evs_vs_No_particle_Evs_4_h", "Pt_Evs_vs_No_particle_Evs_24_h",
                          "Au_Evs_over_time", "Au_Cit_over_time", "Au_Lys_over_time",
                          "Pd_over_time", "Pt_over_time", 
                          "Pd_Evs_over_time", "Pt_Evs_over_time")
  names(combinations) <- names_combinations
  
  fill <- c("chartreuse1","coral","yellow")
  
  for (comb in seq_along(combinations)) {
    
    name_labels <- names(combinations[[comb]])
    venn_list_plt <- combinations[[comb]]
    name <- names_combinations[comb]
    venn_list_eulerr <- venn(venn_list_plt)
    
    
    venn1 <- plot(venn_list_eulerr,quantities=list(type = c("counts", "percent"),
                                                  col="black", font=4, round=2, cex=0.5),
                  fills = list(fill = fill[1:length(combinations[[1]])], alpha = 0.6),
                  labels = list(col = "black", fontsize = 12),
                  col="black",
                  lty = 4:1,
                  lwd=3,
                  legend = list(labels = name_labels))
    venn1
    venn2 <- plot(venn_list_eulerr)
    
    name_euler_percent <- paste0("Euler_diagram_quant_",names(vec))
    name_euler_percent <- gsub("%","_percent",name_euler_percent)

    create_dir(file.path(output_dir,"Figures",venn_dir,name_euler_percent))
    region_dir <- file.path(output_dir,"Figures",venn_dir,name_euler_percent,name)
    create_dir(region_dir)
    
    pdf(file.path(output_dir,"Figures",venn_dir,name_euler_percent,paste0(name,".pdf")),width=9.5,height=5)
    print(venn1)
    print(venn2)
    dev.off()
    
    venn3 <- gplots::venn(venn_list_plt)
    regions <- attributes(venn3)$intersections
    
    for (nombre in names(regions)) {
      write.csv(regions[[nombre]], file = file.path(region_dir, paste0(nombre, ".csv")), 
                row.names = FALSE, quote = FALSE)
      writeLines(regions[[nombre]], file.path(region_dir, paste0(nombre, ".txt")))
    }
  }
  
  
  
  
  


  
  
  
  
  
  
  
  
  
  