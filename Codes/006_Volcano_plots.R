############################################################################
############################################################################
###                                                                      ###
###                       EXPRESSION LEVELS SCRIPT                       ###
###                                                                      ###
############################################################################
############################################################################

###########################
##  0. Load dependences  ##
###########################
{
    library(tidyverse)
    library(ggplot2)
    library(EnhancedVolcano)
    library(writexl)
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
  volcan_plot <- function(data,x,y,xintercept,th,ymax,title){
    
    datos <- data.frame(data)
    y <- -(log10(y))
    y[is.infinite(y)] <- NA
    
    sizes <- c("UP" = 2, "DOWN" = 2, "Not Sig." = 0.5) 
    alphas <- c("UP" = 1, "DOWN" = 1, "Not Sig." = 0.5)
    #FFA373
    ##50486D
    mycolors <- c("UP" ="green","Not Sig." ="grey","DOWN" ="magenta")
    
    ggplot(datos,aes(x=x,y=y,
                     size=Association,
                     alpha=Association,
                     fill=Association))+
      geom_point(shape = 21, colour = "black", size = 3, alpha = 0.8)+
      ylab("-log10(FDR)")+xlab("Log2FC")+theme_minimal()+
      scale_fill_manual(values=mycolors)+
      scale_size_manual(values = sizes,guide="none")+
      scale_alpha_manual(values = alphas, guide="none")+
      geom_vline(xintercept=c(-xintercept, xintercept),col="black",linetype = "dashed")+
      geom_hline(yintercept=-log10(th), col="black",linetype = "dashed")+
      scale_x_continuous(limits = c(-(max(x)), max(x)))+
      scale_y_continuous(limits = c(-1, ymax),expand = expansion(0))+
      ggtitle(title)+
      theme(plot.title = element_text(hjust = 0.5),
            text = element_text(size = 16),
            axis.title = element_text(size = 18),
            axis.text = element_text(size = 14))
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

deg_0.2_files <- list.files(file.path(input_dir,"003_DEG"),pattern = ">0.2")
deg_0.5_files <- list.files(file.path(input_dir,"003_DEG"),pattern = ">0.5")
deg_1_files <- list.files(file.path(input_dir,"003_DEG"),pattern = ">1")

deg_1_list <- lapply(deg_1_files, function(x) read.table(file.path(input_dir,"003_DEG",x)))
deg_1_list <- lapply(deg_1_list, function(x) x[,c("logFC", "adj.P.Val")])

deg_0.5_list<-lapply(deg_0.5_files,function(name){
  x <- try(read.table(file.path(input_dir,"003_DEG",name)))
  if(inherits(x, "try-error"))
    return(NULL)
  else
    return(x)
})

deg_0.5_list <- Filter(Negate(is.null), deg_0.5_list)
deg_0.5_list <- lapply(deg_0.5_list, function(x) x[,c("logFC", "adj.P.Val")])


deg_0.2_list<-lapply(deg_0.2_files,function(name){
  x <- try(read.table(file.path(input_dir,"003_DEG",name)))
  if(inherits(x, "try-error"))
    return(NULL)
  else
    return(x)
})

deg_0.2_list <- Filter(Negate(is.null), deg_0.2_list)
deg_0.2_list <- lapply(deg_0.2_list, function(x) x[,c("logFC", "adj.P.Val")])

### Set the threshold
th=0.05
create_dir(file.path(output_dir,"Figures","006_Volcano_plot"))

# deg_0.2_list
# deg_0.5_list
# deg_1_list

list_to_plt <- deg_1_list
nombre <- deg_1_files
nombre <- gsub(".txt","",nombre)

for (index in seq_along(list_to_plt)) {
  
  df <- list_to_plt[[index]]
  
  name <- nombre[index]
  df$Association <- "Background"
  df[which(df$logFC>0),]$Association <- "UP"
  df[which(df$logFC<0),]$Association <- "DOWN"
  df[which(df$adj.P.Val > th),"Association"] <- "Not Sig."
  deg <- dim(df[which(df$Association!="Not Sig."),])[1]
  
  print(paste0 (name, " there are ",deg))
  #######################
  ### 3. Volcano plot ###
  #######################
  
  genes <- rownames(df[which(df$Association=="UP" |
                                          df$Association=="DOWN"),])

  df_genes <- df[which(df$Association=="UP" | 
                                    df$Association=="DOWN"),]

  ymax <- max(-(log10(df$adj.P.Val)))
  
  plt <- volcan_plot(df,x=df$logFC ,y=df$adj.P.Val, xintercept = 0.0,th=th,ymax=10,title = name )
  
  # pdf(file = file.path(output_dir,"Figures_paper","volcano_plot_DE_lipids_Fig_3B.pdf"),width=10,height=12)
  print(plt)
  # dev.off()
  
  
  plt1 <-EnhancedVolcano(df,x="logFC",y="adj.P.Val",lab = rownames(df),
                         title = "Volcano plot",
                         subtitle = bquote(italic(LCFA ~interaction)),
                         ylab = bquote(~-log[10]~ FDR),
                         legendLabels = c("NS", expression(Log[2] ~ FC), "FDR", expression(~-log[10]~ FDR)),
                         pCutoff = 0.05, FCcutoff = 0,
                         col = c( "grey30","grey30","royalblue","red"),
                         pointSize = 2,
                         colAlpha = 0.9,
                         drawConnectors = T,
                         arrowheads = F,
                         max.overlaps=20)
  plt1 <- plt1 +
    # geom_point(size = 3,color="red")+
    theme_bw() +
    theme(plot.title = element_text(size = 24, hjust = 0.5),
          plot.subtitle = element_text(size = 18, hjust = 0.5),
          axis.title = element_text(size = 24),
          legend.title = element_blank(),
          legend.text = element_text(size = 14),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 2)
    )
  
  plt1
  
  pdf(file = file.path(output_dir,"Figures","006_Volcano_plot",paste0(name,".pdf")),width=10,height=12)
  print(plt1)
  dev.off()
  
  
}



















