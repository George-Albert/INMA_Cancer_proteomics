
############################################################################
############################################################################
###                                                                      ###
###                   DIFFERENTIAL EXPRESSION ANALYSIS                   ###
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
}


#################
##  Functions  ##
################# 
{
create_dir=function(x){suppressWarnings(dir.create(x,recursive=TRUE))}
drcols=function(x){data.frame(colnames(x))}
ul=function(x,n=5){x[1:min(nrow(x),n),1:min(ncol(x),n)]}
options(width=1000)
my_name <- function(v1) {
  deparse(substitute(v1))
}
}

###################################
##  Setting wd and loading data  ##
###################################
getwd()
main_wd <-  "/Users/George_A/Documents/Portafolio/INMA_cancer/Analyses"
setwd(main_wd)
input_dir <-"Inputs" 
output_dir <- "Outputs"

feature_data=read.table(file.path(input_dir,"txt","feature_data.txt"))
metadata=read.table(file.path(input_dir,"txt","metadata_filtered.txt"))
reads_spec=read.table(file.path(input_dir,"txt","reads_spec_filtered.txt"))
reads_empai=read.table(file.path(input_dir,"txt","reads_empai.txt"))
reads_empai <- reads_empai[,colnames(reads_spec)]


length(which(colnames(reads_spec) != colnames(reads_empai)))
write.table(reads_empai, file.path(input_dir,"txt","reads_empai_filtered.txt"))

##################################################################
##                     SET INPUT PARAMETERS                     ##
##################################################################
reads_vec <- list(reads_empai,reads_spec)
reads <- reads_vec[[2]]

remove_extreme_vec <- c("variance","mean","both")
rm_ext <- remove_extreme_vec[3]

th_mean_vec <- c(0.2, 0.5, 1)
# th_mean <- th_mean_vec[1]
  
########################################
##  Differential Expression Analysis  ##
########################################

# Log Transformation
exp=log2(reads+1)
plotDensities(exp,legend=F)

exp_norm=normalize.quantiles.robust(as.matrix(exp),copy=FALSE, 
                                    remove.extreme=rm_ext,
                                    n.remove=1,use.median=FALSE,
                                    use.log2=FALSE)
exp_norm <- data.frame(exp_norm)

# ?normalizeQuantiles()
plotDensities(exp_norm,legend=F)

### number of samples and genes
n_samples <- ncol(reads)
n_genes <- nrow (reads)
n_Au <- nrow(metadata[which(metadata$Particle=="Au"),])
n_PEG <- n_samples-n_Au
th <- 0.05

### Design matrix 
metadata$short_setup <- factor(metadata$short_setup)
metadata$Time=factor(metadata$Time,levels=c("4h","24h"))
metadata$Particle=factor(metadata$Particle,levels=c("Au","PEG"))

group <- factor(metadata$short_setup)

design=model.matrix(~0+group,data=metadata$short_setup)
colnames(design) <- levels(group)
contrast_vec <- colnames(design)

fit=lmFit(exp_norm,design)
fit2=eBayes(fit,trend=T, robust=T)
plotSA(fit2)

### filter low expression data based on the mean
means=apply(exp_norm,1,mean)

for (th_mean in th_mean_vec ) {
  
  reads_filtered <- reads[which(means>th_mean),]
  
  exp=log2(reads_filtered+1)
  exp_norm_1=normalize.quantiles.robust(as.matrix(exp),copy=FALSE, 
                                        remove.extreme=rm_ext,
                                        n.remove=1,use.median=FALSE,
                                        use.log2=FALSE)
  
  fit=lmFit(exp_norm_1,design)
  fit3=eBayes(fit,trend=T, robust=T)
  
  create_dir(file.path(output_dir,"plot_SA"))
  SA_plt_dir <- paste0("Plot_SA_mean>",th_mean,".pdf")
  
  plotSA(fit3)
  
  pdf(file.path(output_dir,"plot_SA",SA_plt_dir), width = 8, height = 6)
  # 
  plotSA(fit3, xlab = "Average log-expression", ylab = "sqrt(sigma)", zero.weights = FALSE,
         pch = 19, cex = 0.6, col = c("black","red"),cex.axis = 1.5, cex.lab = 1.5)
  dev.off()
  
  ### make the Contrasts
  Au_24hvsAu_4h   <- makeContrasts(Au_24h-Au_4h, levels = design)
  PEG_24hvsPEG_4h <- makeContrasts(PEG_24h-PEG_4h, levels = design)
  PEG_4hvsAu_4h   <- makeContrasts(PEG_4h-Au_4h, levels = design)
  PEG_24hvsAu_24h <- makeContrasts(PEG_24h-Au_24h, levels = design)
  PEGvsAu         <- makeContrasts(PEG_24hvsPEG_4h-Au_24hvsAu_4h, levels = design)
  
  Contrast_matrix <- cbind(Au_24hvsAu_4h,PEG_24hvsPEG_4h,PEG_4hvsAu_4h,PEG_24hvsAu_24h,
                           PEGvsAu)
  colnames(Contrast_matrix) <- c("Au_24hvsAu_4h","PEG_24hvsPEG_4h",
                                 "PEG_4hvsAu_4h","PEG_24hvsAu_24h",
                                 "PEGvsAu")
  total_fit <- contrasts.fit(fit,Contrast_matrix)
  total_fit <- eBayes(total_fit)
  
  n_genes <- nrow(reads_filtered)
  
  deg_dir <- "003_DEG"
  dir.create(file.path(input_dir,deg_dir),recursive = T,showWarnings = F)
  
  deg_Au_24hvsAu_4h   <- topTable(total_fit,coef = "Au_24hvsAu_4h",number = n_genes,adjust.method = "BH")
  deg_PEG_24hvsPEG_4h <- topTable(total_fit,coef = "PEG_24hvsPEG_4h",number = n_genes,adjust.method = "BH")
  deg_PEG_4hvsAu_4h   <- topTable(total_fit,coef = "PEG_4hvsAu_4h",number = n_genes,adjust.method = "BH")
  deg_PEG_24hvsAu_24h <- topTable(total_fit,coef = "PEG_24hvsAu_24h",number = n_genes,adjust.method = "BH")
  deg_PEGvsAu         <- topTable(total_fit,coef = "PEGvsAu",number = n_genes,adjust.method = "BH")
  
  write.table(deg_Au_24hvsAu_4h,file = file.path(input_dir,deg_dir, paste0(my_name(deg_Au_24hvsAu_4h),"mean>",th_mean,".txt")))
  write.table(deg_PEG_24hvsPEG_4h,file = file.path(input_dir,deg_dir, paste0(my_name(deg_PEG_24hvsPEG_4h),"mean>",th_mean,".txt")))
  write.table(deg_PEG_4hvsAu_4h,file = file.path(input_dir,deg_dir, paste0(my_name(deg_PEG_4hvsAu_4h),"mean>",th_mean,".txt")))
  write.table(deg_PEG_24hvsAu_24h,file = file.path(input_dir,deg_dir, paste0(my_name(deg_PEG_24hvsAu_24h),"mean>",th_mean,".txt")))
  write.table(deg_PEGvsAu,file = file.path(input_dir,deg_dir, paste0(my_name(deg_PEGvsAu),"mean>",th_mean,".txt")))
  
  sig_deg_dir <- "004_Significative_DEG"
  dir.create(file.path(input_dir,sig_deg_dir),recursive = T,showWarnings = F)
  
  deg_Au_24hvsAu_4h   <- topTable(total_fit,coef = "Au_24hvsAu_4h",number = n_genes,adjust.method = "BH",p.value = 0.05)
  deg_PEG_24hvsPEG_4h <- topTable(total_fit,coef = "PEG_24hvsPEG_4h",number = n_genes,adjust.method = "BH",p.value = 0.05)
  deg_PEG_4hvsAu_4h   <- topTable(total_fit,coef = "PEG_4hvsAu_4h",number = n_genes,adjust.method = "BH",p.value = 0.05)
  deg_PEG_24hvsAu_24h <- topTable(total_fit,coef = "PEG_24hvsAu_24h",number = n_genes,adjust.method = "BH",p.value = 0.05)
  deg_PEGvsAu         <- topTable(total_fit,coef = "PEGvsAu",number = n_genes,adjust.method = "BH",p.value = 0.05)
  
  write.table(deg_Au_24hvsAu_4h,file = file.path(input_dir,sig_deg_dir, paste0(my_name(deg_Au_24hvsAu_4h),"mean>",th_mean,".txt")))
  write.table(deg_PEG_24hvsPEG_4h,file = file.path(input_dir,sig_deg_dir, paste0(my_name(deg_PEG_24hvsPEG_4h),"mean>",th_mean,".txt")))
  write.table(deg_PEG_4hvsAu_4h,file = file.path(input_dir,sig_deg_dir, paste0(my_name(deg_PEG_4hvsAu_4h),"mean>",th_mean,".txt")))
  write.table(deg_PEG_24hvsAu_24h,file = file.path(input_dir,sig_deg_dir, paste0(my_name(deg_PEG_24hvsAu_24h),"mean>",th_mean,".txt")))
  write.table(deg_PEGvsAu,file = file.path(input_dir,sig_deg_dir, paste0(my_name(deg_PEGvsAu),"mean>",th_mean,".txt")))
  
}





