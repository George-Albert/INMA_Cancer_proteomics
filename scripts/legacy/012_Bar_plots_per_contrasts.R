
############################################################################
############################################################################
###                                                                      ###
###                   BAR PLOT WITH DEG PER CONTRAST                     ###
###                                                                      ###
############################################################################
############################################################################

###########################
##  0. Load dependences  ##
###########################
{
  library(tidyverse)
  library(stringr)
  library(tools)
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
main_wd <-  getwd()
setwd(main_wd)
input_dir <-"Inputs" 
output_dir <- "Outputs"
DEG_dir <- "009_Significative_DEG_added"

deg_mean_0.2 <- list.files(path = file.path(input_dir,DEG_dir,"mean_0.2","txt"),pattern = "txt")
deg_mean_0.5 <- list.files(path = file.path(input_dir,DEG_dir,"mean_0.5","txt"),pattern = "txt")
deg_mean_1 <- list.files(path = file.path(input_dir,DEG_dir,"mean_1","txt"),pattern = "txt")

list_of_files_mean_0.2 <- lapply(deg_mean_0.2, function(x) read.table(file.path(input_dir,DEG_dir,"mean_0.2","txt",x)))
list_of_files_mean_0.5 <- lapply(deg_mean_0.5, function(x) read.table(file.path(input_dir,DEG_dir,"mean_0.5","txt",x)))
list_of_files_mean_1 <- lapply(deg_mean_1, function(x) read.table(file.path(input_dir,DEG_dir,"mean_1","txt",x),blank.lines.skip = FALSE))

constrast_names_0.2 <- file_path_sans_ext(deg_mean_0.2)
constrast_names_0.5 <- file_path_sans_ext(deg_mean_0.5)
constrast_names_1   <- file_path_sans_ext(deg_mean_1)

list_of_files_mean_0.2 <- setNames(list_of_files_mean_0.2,constrast_names_0.2)
list_of_files_mean_0.5 <- setNames(list_of_files_mean_0.5,constrast_names_0.5)
list_of_files_mean_1 <- setNames(list_of_files_mean_1,constrast_names_1)

# create a table to introduce to histogram
get_num_rows <- function(df) {
  num_rows <- nrow(df)
  if (num_rows == 1 && all(is.na(df))) {
    return(0)
  } else {
    return(num_rows)
  }
}

vec <- list(list_of_files_mean_0.2,list_of_files_mean_0.5,list_of_files_mean_1)

for (i in seq_along(vec)) {
  
  print(i)
  
  num_rows_list <- sapply(vec[[i]], get_num_rows)
  num_rows_list <- data.frame(Contrasts= names(num_rows_list),Number.Genes=num_rows_list)
  ### remove the mean word from the names of the contrasts
  rownames(num_rows_list) <- gsub("_mean>1","",rownames(num_rows_list))
  rownames(num_rows_list) <- gsub("_mean>0.5","",rownames(num_rows_list))
  rownames(num_rows_list) <- gsub("_mean>0.2","",rownames(num_rows_list))
  
  num_rows_list$Contrasts <- rownames(num_rows_list)
  
  # Create histogram of genes per sample before filtering
  num_rows_list$Contrasts <- factor(num_rows_list$Contrasts,levels = num_rows_list$Contrasts)
  
  hist_plt_mean <- ggplot(num_rows_list, aes(x = Contrasts, y = Number.Genes,)) +
    geom_col(color="black",fill="skyblue") +
    geom_text(aes(label = Number.Genes), 
              vjust = -0.5, 
              size = 3.5, 
              color = "black") +  # Add labels
    labs(title = "Number of Genes by Contrast",
         x = "Contrast",
         y = "Number of Genes") +
    # ylim(0,150)+
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  hist_plt_mean
  
  create_dir(file.path(output_dir,"Figures","Bar_plot_contrasts"))
  mean_vec <- c(0.2,0.5,1)
  
  hist <- file.path(output_dir,"Figures","Bar_plot_contrasts",paste0("Number_of_Genes_by_Contrast_mean_",mean_vec[i],"_.pdf"))
  
  pdf(file =hist )
  print(hist_plt_mean)
  dev.off()
}


























