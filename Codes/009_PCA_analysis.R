

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
  library(openxlsx)
  library(preprocessCore)
  library(stringr)
  library(ggrepel)
}


############################
### 1. Declare functions ### 
############################
dcols=function(x){data.frame(colnames(x))}
create_dir=function(x){suppressWarnings(dir.create(x,recursive=TRUE))}

getwd()
main_wd <-  getwd()
setwd(main_wd)
input_dir <-"Inputs" 
output_dir <- "Outputs"

##########################################################
##                     Loading Data                     ##
##########################################################
feature_data=read.table(file.path(input_dir,"txt","feature_data_added.txt"))
metadata=read.table(file.path(input_dir,"txt","metadata_added.txt"))
reads_spec=read.table(file.path(input_dir,"txt","reads_spec_added.txt"))
reads_empai=read.table(file.path(input_dir,"txt","reads_empai_added.txt"))

reads_vec <- list(reads_empai,reads_spec)
reads <- reads_vec[[2]]

remove_extreme_vec <- c("variance","mean","both")
rm_ext <- remove_extreme_vec[3]

# Set all the NAn values to zero
reads[is.na(reads)] <- 0

exp=log2(reads+1)
exp_norm=normalize.quantiles.robust(as.matrix(exp),copy=FALSE, 
                                    remove.extreme=rm_ext,
                                    n.remove=1,use.median=FALSE,
                                    use.log2=FALSE)
th_mean <- 0.5
means=apply(exp_norm,1,mean)
reads_filtered <- reads[which(means>th_mean),]

exp=log2(reads_filtered+1)
exp_norm_1=normalize.quantiles.robust(as.matrix(exp),copy=FALSE, 
                                      remove.extreme=rm_ext,
                                      n.remove=1,use.median=FALSE,
                                      use.log2=FALSE)

### number of samples and genes
n_samples <- ncol(reads)
th <- 0.05

### Design matrix 
metadata$short_setup <- factor(metadata$short_setup,levels = unique(metadata$short_setup))
metadata$Time=factor(metadata$Time,levels=c("5min","4h","24h"))

group <- factor(metadata$short_setup)

design=model.matrix(~0+group,data=metadata$short_setup)
colnames(design) <- levels(group)

tab_name <- "PCA_table_added"

exp <- exp_norm_1
pca <-prcomp(t(exp))
sum_pca=data.frame(summary(pca)$importance[,c(1:5)])
# Create a df with the meta data and the data from PCA
pca_df <- as.data.frame(pca$x)

# Bind the metadata and PCA data
combined_data <- bind_cols(metadata,pca_df[,1:6])
combined_data$Shape <- 21
combined_data[combined_data$Time == "4h",]$Shape <- 8
combined_data[combined_data$Time == "24h",]$Shape <- 2

pca_dir <- "PCA"
dir.create(file.path(output_dir,pca_dir,tab_name),recursive=TRUE,showWarnings = F)
write.table(combined_data,file=file.path(output_dir,pca_dir,tab_name,paste0(tab_name,".txt")))

fill_base <- combined_data$Fill[which(!duplicated(combined_data$short_setup))]
color_base <- combined_data$Color[which(!duplicated(combined_data$short_setup))]
shape <- combined_data$Shape[which(!duplicated(combined_data$Shape))]
label <- factor(unique(combined_data$short_setup),levels = unique(combined_data$short_setup))

combined_data$short_setup <- factor(combined_data$short_setup,levels = unique(combined_data$short_setup))

# Generate the plot
pca_plt <- ggplot(combined_data, aes(x = PC1, y = PC2, fill = short_setup, 
                                     color = short_setup, shape = Time)) +
  geom_point(size = 2,stroke=1.5) +
  # geom_label_repel(aes(label = short_setup),color="black",fill="white",
  #                  max.overlaps = 24, size = 1, nudge_x = 0.1, nudge_y = 0.1,
  #                  show.legend = F) +
  theme_bw() +
  scale_fill_manual(values = fill_base, name = "Samples",labels=label) +
  scale_color_manual(values = color_base, name = "Samples",labels=label) +
  scale_shape_manual(values = shape, name = "Time",labels=unique(combined_data$Time)) +
  xlab(paste0("PC1:",100*round(sum_pca$PC1[2],digits=3),"% variance explained"))+
  ylab(paste0("PC2:",100*round(sum_pca$PC2[2],digits=3),"% variance explained"))+theme_bw()+
  theme(
    axis.text.y   = element_text(size=14),
    axis.text.x   = element_text(size=14),
    axis.title.y  = element_text(size=14),
    axis.title.x  = element_text(size=14),
    # panel.background = element_blank(),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=1,linetype="solid"),
    # legend.title=element_blank(),
    # legend.position="none",
    legend.text=element_text(size=14),
    legend.key.size = unit(0.3, 'lines'))

print(pca_plt)
pca_name <- "PCA_all_particles_added_wo_label"
dir.create(file.path(output_dir,pca_dir),recursive=TRUE,showWarnings = F)
pdf(file = file.path(output_dir,pca_dir,paste0(pca_name,".pdf")),width=10,height=5)
print(pca_plt)
dev.off()

#########################################
##  PCA Analysis without the outliers  ##
#########################################

### Remove outlier samples 9,4,14
metadata_wo_outlier <- metadata[!(metadata$Sample.Names %in% c("Evs_Au_Control_1", "Evs_Control_1" , "Evs_Pd_Control_1",
                                                               "Evs_Pt_Control_1")),]
reads_wo_outlier <- reads_filtered[,rownames(metadata_wo_outlier)]

exp_wo_outlier=log2(reads_wo_outlier+1)
exp_norm_wo_outlier=normalize.quantiles.robust(as.matrix(exp_wo_outlier),copy=FALSE, 
                                    remove.extreme=rm_ext,
                                    n.remove=1,use.median=FALSE,
                                    use.log2=FALSE)

### Design matrix 
metadata_wo_outlier$short_setup <- factor(metadata_wo_outlier$short_setup,levels = unique(metadata_wo_outlier$short_setup))
metadata_wo_outlier$Time=factor(metadata_wo_outlier$Time,levels=c("5min","4h","24h"))

group <- factor(metadata_wo_outlier$short_setup)

design=model.matrix(~0+group,data=metadata_wo_outlier$short_setup)
colnames(design) <- levels(group)

tab_name <- "PCA_table_wo_outlier"
pca_name <- "PCA_all_particles_wo_outlier"

exp <- exp_norm_wo_outlier
pca <-prcomp(t(exp))
sum_pca=data.frame(summary(pca)$importance[,c(1:5)])
# Create a df with the meta data and the data from PCA
pca_df <- as.data.frame(pca$x)

# Bind the metadata and PCA data
combined_data <- bind_cols(metadata_wo_outlier,pca_df[,1:6])
combined_data$Shape <- 21
combined_data[combined_data$Time == "4h",]$Shape <- 8
combined_data[combined_data$Time == "24h",]$Shape <- 2

dir.create(file.path(output_dir,pca_dir,"PCA_table_added"),recursive=TRUE,showWarnings = F)
write.table(combined_data,file=file.path(output_dir,pca_dir,"PCA_table_added",paste0(tab_name,".txt")))

fill_base <- combined_data$Fill[which(!duplicated(combined_data$short_setup))]
color_base <- combined_data$Color[which(!duplicated(combined_data$short_setup))]
shape <- combined_data$Shape[which(!duplicated(combined_data$Shape))]
label <- factor(unique(combined_data$short_setup),levels = unique(combined_data$short_setup))
combined_data$short_setup <- factor(combined_data$short_setup,levels = unique(combined_data$short_setup))

# Generate the plot
pca_plt <- ggplot(combined_data, aes(x = PC1, y = PC2, fill = short_setup, 
                                     color = short_setup,shape = Time)) +
  geom_point(size = 2,stroke=1.5) +
  # geom_label_repel(aes(label = short_setup),color="black",fill="white",
  #                  max.overlaps = 21, size = 2, nudge_x = 0.1, nudge_y = 0.1,
  #                  show.legend = F) +
  theme_bw() +
  scale_fill_manual(values = fill_base,name = "Samples",labels=label) +
  scale_color_manual(values = color_base,name = "Samples",labels=label) +
  scale_shape_manual(values = shape, name = "Time",labels=unique(combined_data$Time)) +
  xlab(paste0("PC1:",100*round(sum_pca$PC1[2],digits=3),"% variance explained"))+
  ylab(paste0("PC2:",100*round(sum_pca$PC2[2],digits=3),"% variance explained"))+theme_bw()+
  theme(
    axis.text.y   = element_text(size=14),
    axis.text.x   = element_text(size=14),
    axis.title.y  = element_text(size=14),
    axis.title.x  = element_text(size=14),
    # panel.background = element_blank(),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=1,linetype="solid"),
    # legend.title=element_blank(),
    # legend.position="none",
    legend.text=element_text(size=14),
    legend.key.size = unit(1, 'lines'))

print(pca_plt)

dir.create(file.path(output_dir,pca_dir),recursive=TRUE,showWarnings = F)
pdf(file = file.path(output_dir,pca_dir,paste0(pca_name,"_added.pdf")),width=12,height=8)
print(pca_plt)
dev.off()


write.table(metadata_wo_outlier, file.path(input_dir,"txt","metadata_filtered_added.txt"))
write.table(reads_wo_outlier, file.path(input_dir,"txt","reads_spec_filtered_added.txt"))













