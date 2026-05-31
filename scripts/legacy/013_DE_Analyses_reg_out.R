
############################################################################
############################################################################
###                                                                      ###
###     DIFFERENTIAL EXPRESSION ANALYSIS REGRESSION OUT                  ###
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

main_wd <-  getwd()
setwd(main_wd)
input_dir <-"Inputs" 
output_dir <- "Outputs"

feature_data=read.table(file.path(input_dir,"txt","feature_data_added.txt"))
metadata=read.table(file.path(input_dir,"txt","metadata_filtered_added.txt"))
reads_spec=read.table(file.path(input_dir,"txt","reads_spec_added.txt"))
reads_spec <- reads_spec[,rownames(metadata)]
reads_spec[is.na(reads_spec)] <- 0
length(which(colnames(reads_spec) != rownames(metadata)))
write.table(reads_spec, file.path(input_dir,"txt","reads_spec_added_filtered.txt"))

##################################################################
##                     SET INPUT PARAMETERS                     ##
##################################################################

reads <- reads_spec

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
th <- 0.05

### Design matrix 
metadata$short_setup <- factor(metadata$short_setup,levels = unique(metadata$short_setup))
metadata$Time=factor(metadata$Time,levels=c("5min","4h","24h"))
metadata$Sample.Order <- factor(metadata$Sample.Order)
group <- factor(metadata$short_setup)

design=model.matrix(~0+group+Sample.Order,data=metadata)
colnames(design)[1:24]<- levels(group)
contrast_vec <- colnames(design)

fit=lmFit(exp_norm,design)
fit2=eBayes(fit,trend=T, robust=T)
plotSA(fit2)

### filter low expression data based on the mean
means=apply(exp_norm,1,mean)

####################################
##  Make Contrasts (21 in total)  ##
####################################
{
##############################
##  Pd vs Pd_Evs x 3 times  ##
##############################
Pd_5minvsPd_Evs_5min   <- makeContrasts(Pd_Evs_Yes_5min - Pd_Evs_NO_5min, levels = design)
Pd_4hvsPd_Evs_4h   <- makeContrasts(Pd_Evs_Yes_4h - Pd_Evs_NO_4h, levels = design)
Pd_24hvsPd_Evs_24h   <- makeContrasts(Pd_Evs_Yes_24h - Pd_Evs_NO_24h, levels = design)

##############################
##  Pt vs Pt_Evs x 3 times  ##
##############################
Pt_5minvsPt_Evs_5min   <- makeContrasts(Pt_Evs_Yes_5min - Pt_Evs_NO_5min, levels = design)
Pt_4hvsPt_Evs_4h   <- makeContrasts(Pt_Evs_Yes_4h - Pt_Evs_NO_4h, levels = design)
Pt_24hvsPt_Evs_24h   <- makeContrasts(Pt_Evs_Yes_24h - Pt_Evs_NO_24h, levels = design)

##################################
##  Au_Cit vs Au_Evs x 3 times  ##
##################################
Au_Cit_5minvsAu_Evs_5min   <- makeContrasts(Au_Evs_Yes_5min - Au_Cit_Evs_NO_5min, levels = design)
Au_Cit_4hvsAu_Evs_4h   <- makeContrasts(Au_Evs_Yes_4h - Au_Cit_Evs_NO_4h, levels = design)
Au_Cit_24hvsAu_Evs_24h   <- makeContrasts(Au_Evs_Yes_24h - Au_Cit_Evs_NO_24h, levels = design)

##################################
##  Au_Lys vs Au_Evs x 3 times  ##
##################################
Au_Lys_5minvsAu_Evs_5min   <- makeContrasts(Au_Evs_Yes_5min - Au_Lys_Evs_NO_5min, levels = design)
Au_Lys_4hvsAu_Evs_4h   <- makeContrasts(Au_Evs_Yes_4h - Au_Lys_Evs_NO_4h, levels = design)
Au_Lys_24hvsAu_Evs_24h   <- makeContrasts(Au_Evs_Yes_24h - Au_Lys_Evs_NO_24h, levels = design)

#################################################
##  Particle_Evs vs No_Particle_Evs x 3 times  ##
#################################################
Au_5minvsNo_particle_5min   <- makeContrasts(Au_Evs_Yes_5min - No_particle_Evs_Yes_5min, levels = design)
Pd_5minvsNo_particle_5min    <- makeContrasts(Pd_Evs_Yes_5min - No_particle_Evs_Yes_5min, levels = design)
Pt_5minvsNo_particle_5min    <- makeContrasts(Pt_Evs_Yes_5min - No_particle_Evs_Yes_5min, levels = design)

Au_4hvsNo_particle_4h   <- makeContrasts(Au_Evs_Yes_4h - No_particle_Evs_Yes_4h, levels = design)
Pd_4hvsNo_particle_4h    <- makeContrasts(Pd_Evs_Yes_4h - No_particle_Evs_Yes_4h, levels = design)
Pt_4hvsNo_particle_4h    <- makeContrasts(Pt_Evs_Yes_4h - No_particle_Evs_Yes_4h, levels = design)

Au_24hvsNo_particle_24h   <- makeContrasts(Au_Evs_Yes_24h - No_particle_Evs_Yes_24h, levels = design)
Pd_24hvsNo_particle_24h    <- makeContrasts(Pd_Evs_Yes_24h - No_particle_Evs_Yes_24h, levels = design)
Pt_24hvsNo_particle_24h    <- makeContrasts(Pt_Evs_Yes_24h - No_particle_Evs_Yes_24h, levels = design)

Contrast_matrix <- cbind(Pd_5minvsPd_Evs_5min,Pd_4hvsPd_Evs_4h,Pd_24hvsPd_Evs_24h,
                         Pt_5minvsPt_Evs_5min,Pt_4hvsPt_Evs_4h,Pt_24hvsPt_Evs_24h,
                         Au_Cit_5minvsAu_Evs_5min,Au_Cit_4hvsAu_Evs_4h,Au_Cit_24hvsAu_Evs_24h,
                         Au_Lys_5minvsAu_Evs_5min,Au_Lys_4hvsAu_Evs_4h,Au_Lys_24hvsAu_Evs_24h,
                         Au_5minvsNo_particle_5min,Pd_5minvsNo_particle_5min,Pt_5minvsNo_particle_5min,
                         Au_4hvsNo_particle_4h,Pd_4hvsNo_particle_4h,Pt_4hvsNo_particle_4h,
                         Au_24hvsNo_particle_24h,Pd_24hvsNo_particle_24h,Pt_24hvsNo_particle_24h
)
colnames(Contrast_matrix) <- gsub(" - ", "vs",colnames(Contrast_matrix))
}

for (th_mean in th_mean_vec ) {
  
  reads_filtered <- reads[which(means>th_mean),]
  
  exp=log2(reads_filtered+1)
  exp_norm_1=normalize.quantiles.robust(as.matrix(exp),copy=FALSE, 
                                        remove.extreme=rm_ext,
                                        n.remove=1,use.median=FALSE,
                                        use.log2=FALSE)
  
  fit=lmFit(exp_norm_1,design)
  fit3=eBayes(fit,trend=T, robust=T)
  
  ### Create the beta matrix
  betas <- fit3$coefficients
  ### Select the batch effect to delete
  lote_col_to_delete <- grep(pattern = "Sample.Order", colnames(betas))
  ### get the new expression data
  exp_clean=exp_norm_1 - as.matrix(betas[,lote_col_to_delete]) %*% t(as.matrix(design[,lote_col_to_delete]))
  
  create_dir(file.path(output_dir,"plot_SA"))
  SA_plt_dir <- paste0("Plot_SA_mean_reg_out>",th_mean,".pdf")
  
  plotSA(fit3)
  
  pdf(file.path(output_dir,"plot_SA",SA_plt_dir), width = 8, height = 6)
  # 
  plotSA(fit3, xlab = "Average log-expression", ylab = "sqrt(sigma)", zero.weights = FALSE,
         pch = 19, cex = 0.6, col = c("black","red"),cex.axis = 1.5, cex.lab = 1.5)
  dev.off()
  
  total_fit <- contrasts.fit(fit,Contrast_matrix)
  total_fit <- eBayes(total_fit)
  
  n_genes <- nrow(reads_filtered)
  
  deg_dir <- "011_DEG_reg_out"
  sig_deg_dir <- "012_Significative_DEG_reg_out"
  
  dir.create(file.path(input_dir,deg_dir),recursive = T,showWarnings = F)
  dir.create(file.path(input_dir,sig_deg_dir),recursive = T,showWarnings = F)
  
  contrasts <- colnames(Contrast_matrix)
  
  create_dir(file.path(input_dir,deg_dir,paste0("mean_",th_mean),"txt"))
  create_dir(file.path(input_dir,deg_dir,paste0("mean_",th_mean),"csv"))
  
  create_dir(file.path(input_dir,sig_deg_dir,paste0("mean_",th_mean),"txt"))
  create_dir(file.path(input_dir,sig_deg_dir,paste0("mean_",th_mean),"csv"))
  
  # Iterate over each contrast and save the results
  for (contrast in contrasts) {
    # Extract the topTable results for the current contrast
    print(contrast)
    result <- topTable(total_fit, coef = contrast, number = n_genes, adjust.method = "BH")
    result_1 <- topTable(total_fit, coef = contrast, number = n_genes, adjust.method = "BH",p.value = 0.05)
    
    # Create the file name
    file_name_csv <- file.path(input_dir,deg_dir, paste0("mean_",th_mean),"csv",paste0(contrast,"_mean>",th_mean, ".csv"))
    file_name_txt <- file.path(input_dir,deg_dir,paste0("mean_",th_mean),"txt", paste0(contrast,"_mean>",th_mean,".txt"))
    
    file_name_csv.1 <- file.path(input_dir,sig_deg_dir, paste0("mean_",th_mean),"csv",paste0(contrast,"_mean>",th_mean, ".csv"))
    file_name_txt.1 <- file.path(input_dir,sig_deg_dir,paste0("mean_",th_mean),"txt", paste0(contrast,"_mean>",th_mean,".txt"))
    
    
    # Save the results to a CSV file and to txt file
    write.csv(result, file = file_name_csv, row.names = T)
    write.table(result,file = file_name_txt)
    
    write.csv(result_1, file = file_name_csv.1, row.names = T)
    write.table(result_1,file = file_name_txt.1)
    
  }
}





