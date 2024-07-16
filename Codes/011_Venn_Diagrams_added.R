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

Au_4h_exp_levels   <- data.frame(row.names=rownames(betas),betas$Au_4h)
PEG_4h_exp_levels  <- data.frame(row.names=rownames(betas),betas$PEG_4h)
Au_24h_exp_levels  <- data.frame(row.names=rownames(betas),betas$Au_24h)
PEG_24h_exp_levels <- data.frame(row.names=rownames(betas),betas$PEG_24h)

# Set the threshold 
umbrales=quantile(exp_norm,c(0.1,0.25,0.5,0.75,0.9))
print(umbrales)
vec <- umbrales[3:5]

for (umbral in seq_along(vec)) {
  
  exp_genes_Au_4h=rownames(Au_4h_exp_levels)[which(Au_4h_exp_levels>vec[umbral])]
  exp_genes_PEG_4h=rownames(PEG_4h_exp_levels)[which(PEG_4h_exp_levels>vec[umbral])]
  exp_genes_Au_24h=rownames(Au_24h_exp_levels)[which(Au_24h_exp_levels>vec[umbral])]
  exp_genes_PEG_24h=rownames(PEG_24h_exp_levels)[which(PEG_24h_exp_levels>vec[umbral])]
  
  venn_dir <- "005_Venn_diagram"
  create_dir(file.path(input_dir,venn_dir))
  quantil <- 
  write.table(exp_genes_Au_4h, file.path(input_dir,venn_dir,paste0(my_name(exp_genes_Au_4h),
                                                                   "_quant_",names(vec[umbral]),".txt")))
  write.table(exp_genes_PEG_4h, file.path(input_dir,venn_dir,paste0(my_name(exp_genes_PEG_4h),
                                                                    "_quant_",names(vec)[umbral],".txt")))
  write.table(exp_genes_Au_24h, file.path(input_dir,venn_dir,paste0(my_name(exp_genes_Au_24h),
                                                                    "_quant_",names(vec)[umbral],".txt")))
  write.table(exp_genes_PEG_24h, file.path(input_dir,venn_dir,paste0(my_name(exp_genes_PEG_24h),
                                                                     "_quant_",names(vec)[umbral],".txt")))
  
  create_dir(file.path(output_dir,"Figures",venn_dir))
  
  venn_list <- list(Au_4h=exp_genes_Au_4h,PEG_4h=exp_genes_PEG_4h,
                    Au_24h=exp_genes_Au_24h,PEG_24h=exp_genes_PEG_24h)
  
  name_euler <- paste0("Euler_diagram_quant_",names(vec)[umbral],".pdf")
  name_euler <- gsub("%","_percent",name_euler)
  
  fill <- c("dodgerblue1","chartreuse1","coral","yellow")
  name_labels <- c("Au_4h","PEG_4","Au_24h","PEG_24")
  
  venn1 <- plot(venn(venn_list),quantities=list(type = c("counts", "percent"),
                                                col="black", font=4, round=2, cex=0.5),
                fills = list(fill = fill, alpha = 0.6),
                labels = list(col = "black", fontsize = 12),
                col="black",
                lty = 4:1,
                lwd=3,
                legend = list(labels = name_labels))
  venn1
  venn2 <- plot(venn(venn_list))
  
  pdf(file.path(output_dir,"Figures",venn_dir,name_euler),width=9.5,height=5)
  print(venn1)
  print(venn2)
  dev.off()
  
  # Access specific regions
  only_Au_4h <- setdiff(setdiff(setdiff(venn_list$Au_4h, venn_list$PEG_4h), venn_list$Au_24h), venn_list$PEG_24h)
  only_PEG_4h <- setdiff(setdiff(setdiff(venn_list$PEG_4h, venn_list$Au_4h), venn_list$Au_24h), venn_list$PEG_24h)
  only_AU_24h <- setdiff(setdiff(setdiff(venn_list$Au_24h, venn_list$Au_4h), venn_list$PEG_4h), venn_list$PEG_24h)
  only_PEG_24h <- setdiff(setdiff(setdiff(venn_list$PEG_24h, venn_list$Au_4h), venn_list$PEG_4h), venn_list$Au_24h)
  
  Au_4h_and_PEG_4h_not_AU_24h_and_PEG_24h <- setdiff(setdiff(intersect(venn_list$Au_4h, venn_list$PEG_4h), venn_list$Au_24h), venn_list$PEG_24h)
  Au_4h_and_AU_24h_not_PEG_4h_and_PEG_24h <- setdiff(setdiff(intersect(venn_list$Au_4h, venn_list$Au_24h), venn_list$PEG_4h), venn_list$PEG_24h)
  Au_4h_and_PEG_24h_not_PEG_4h_and_AU_24h <- setdiff(setdiff(intersect(venn_list$Au_4h, venn_list$PEG_24h), venn_list$PEG_4h), venn_list$Au_24h)
  PEG_4h_and_AU_24h_not_Au_4h_and_PEG_24h <- setdiff(setdiff(intersect(venn_list$PEG_4h, venn_list$Au_24h), venn_list$Au_4h), venn_list$PEG_24h)
  PEG_4h_and_PEG_24h_not_Au_4h_and_AU_24h <- setdiff(setdiff(intersect(venn_list$PEG_4h, venn_list$PEG_24h), venn_list$Au_4h), venn_list$Au_24h)
  AU_24h_and_PEG_24h_not_Au_4h_and_PEG_4h <- setdiff(setdiff(intersect(venn_list$Au_24h, venn_list$PEG_24h), venn_list$Au_4h), venn_list$PEG_4h)
  
  Au_4h_and_PEG_4h_and_AU_24h_not_PEG_24h <- setdiff(intersect(intersect(venn_list$Au_4h, venn_list$PEG_4h), venn_list$Au_24h), venn_list$PEG_24h)
  Au_4h_and_PEG_4h_and_PEG_24h_not_AU_24h <- setdiff(intersect(intersect(venn_list$Au_4h, venn_list$PEG_4h), venn_list$PEG_24h), venn_list$Au_24h)
  Au_4h_and_AU_24h_and_PEG_24h_not_PEG_4h <- setdiff(intersect(intersect(venn_list$Au_4h, venn_list$Au_24h), venn_list$PEG_24h), venn_list$PEG_4h)
  PEG_4h_and_AU_24h_and_PEG_24h_not_Au_4h <- setdiff(intersect(intersect(venn_list$PEG_4h, venn_list$Au_24h), venn_list$PEG_24h), venn_list$Au_4h)
  
  all_four <- intersect(intersect(intersect(venn_list$Au_4h, venn_list$PEG_4h), venn_list$Au_24h), venn_list$PEG_24h)
  
  # List of regions and their names
  regions <- list(
    only_Au_4h = only_Au_4h,
    only_PEG_4h = only_PEG_4h,
    only_AU_24h = only_AU_24h,
    only_PEG_24h = only_PEG_24h,
    Au_4h_and_PEG_4h_not_AU_24h_and_PEG_24h = Au_4h_and_PEG_4h_not_AU_24h_and_PEG_24h,
    Au_4h_and_AU_24h_not_PEG_4h_and_PEG_24h = Au_4h_and_AU_24h_not_PEG_4h_and_PEG_24h,
    Au_4h_and_PEG_24h_not_PEG_4h_and_AU_24h = Au_4h_and_PEG_24h_not_PEG_4h_and_AU_24h,
    PEG_4h_and_AU_24h_not_Au_4h_and_PEG_24h = PEG_4h_and_AU_24h_not_Au_4h_and_PEG_24h,
    PEG_4h_and_PEG_24h_not_Au_4h_and_AU_24h = PEG_4h_and_PEG_24h_not_Au_4h_and_AU_24h,
    AU_24h_and_PEG_24h_not_Au_4h_and_PEG_4h = AU_24h_and_PEG_24h_not_Au_4h_and_PEG_4h,
    Au_4h_and_PEG_4h_and_AU_24h_not_PEG_24h = Au_4h_and_PEG_4h_and_AU_24h_not_PEG_24h,
    Au_4h_and_PEG_4h_and_PEG_24h_not_AU_24h = Au_4h_and_PEG_4h_and_PEG_24h_not_AU_24h,
    Au_4h_and_AU_24h_and_PEG_24h_not_PEG_4h = Au_4h_and_AU_24h_and_PEG_24h_not_PEG_4h,
    PEG_4h_and_AU_24h_and_PEG_24h_not_Au_4h = PEG_4h_and_AU_24h_and_PEG_24h_not_Au_4h,
    all_four = all_four
  ) 
  
  #define file name
  name_regions <- paste0("Regions_Venn_diagram_",names(vec)[umbral],".txt")
  name_regions <- gsub("%","_percent",name_regions)
  
  sink(file.path(input_dir,"txt",name_regions))
  print(regions)
  sink()
  
}






  
  
  
  
  
  
  
  
  
  