#This script creates a feature dictionary pairing RV_numbers, Gene_names and protein GI_number accessions.
# 0. Load dependences
# 1. ID conversion
## following https://support.bioconductor.org/p/64484/
# 0. Load dependences

## PARAMETROS
{
    library(readxl)
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
    
    create_dir=function(x){suppressWarnings(dir.create(x,recursive=TRUE))}
    dcols=function(x){data.frame(colnames(x))}
    ul=function(x,n=5){x[1:min(nrow(x),n),1:min(ncol(x),n)]}
    options(width=1000)
}


## Load the data, transform and filter lowly detected proteins.

{
    
    reads=read.table("Inputs/3_Merged/spectral_reads.txt")
    feature_data=read.table("Inputs/3_Merged/feature_data.txt")
    meta_data=read.table("Inputs/3_Merged/meta_data.txt")
    
    exp=log2(reads+1)
    
    exp=normalize.quantiles.robust(as.matrix(exp),copy=FALSE, remove.extreme="variance",n.remove=1,use.median=FALSE,use.log2=FALSE)
    
    meta_data$Time=factor(meta_data$Time,levels=c("4h","24h"))
    design <- model.matrix(~Particle+Time:Particle,data=meta_data)

    
    fit=lmFit(exp,design)
    fit=eBayes(fit,trend=TRUE, robust=TRUE)

    plotSA(fit)
    
    ## Tenemos que filtrar low exp.
    
    means=apply(exp,1,mean)
    
    ## Exploro la opcion liberal.
    reads=reads[which(means>0.2),]
    exp=log2(reads+1)
    
    exp=normalize.quantiles.robust(as.matrix(exp),copy=FALSE, remove.extreme="variance",n.remove=1,use.median=FALSE,use.log2=FALSE)
    design <- model.matrix(~Particle+Time:Particle,data=meta_data)

    
    fit=lmFit(exp,design)
    fit=eBayes(fit,trend=TRUE, robust=TRUE)

    plotSA(fit)
    
    ## 0.1 es 0, ergo el primero es cualquier expresión no nula.
    umbrales=quantile(exp,c(0.1,0.25,0.5,0.75,0.9))
    
    ## Contraste q captura el nivel de expresión de Au 4h
    vec_1=c(1,0,0,0)
    ## Contraste q captura el nivel de expresión de Au 24h
    vec_2=c(1,0,1,0)
    ## Contraste q captura el nivel de expresión de PEG 4h
    vec_3=c(1,1,0,0)
    ## Contraste q captura el nivel de expresión de PEG 24h
    vec_4=c(1,1,0,1)
    
    fit_1 <- contrasts.fit(fit, vec_1)
    fit_2 <- contrasts.fit(fit, vec_2)
    fit_3 <- contrasts.fit(fit, vec_3)
    fit_4 <- contrasts.fit(fit, vec_4)
    
    expresiones_promedio_1=coefficients(fit_1)
    expresiones_promedio_2=coefficients(fit_2)
    expresiones_promedio_3=coefficients(fit_3)
    expresiones_promedio_4=coefficients(fit_4)

    
    genes_expresados_en_Au_4h=rownames(expresiones_promedio_1)[which(expresiones_promedio_1>umbrales[1])]
    genes_expresados_en_Au_24h=rownames(expresiones_promedio_2)[which(expresiones_promedio_2>umbrales[1])]
    genes_expresados_en_PEG_4h=rownames(expresiones_promedio_3)[which(expresiones_promedio_3>umbrales[1])]
    genes_expresados_en_PEG_24h=rownames(expresiones_promedio_4)[which(expresiones_promedio_4>umbrales[1])]

    ## Ahora haríamos un venn diagram y lo repetiríamos para cada nivel del umbral.
    
    
    ## Etc.--> con esto, en toda muestra, hacemos un Venn diagram por umbral.
    
    pca <-prcomp(t(exp))
    sum_pca=data.frame(summary(pca)$importance[,c(1:5)])

    datos <- data.frame(pca$x)
    colnames(datos)=paste0("PC",c(1:ncol(datos)))
    length(which(rownames(datos)!=rownames(meta_data)))
    datos=cbind(meta_data,datos[,1:6])
    datos$Setup=paste0(datos$Particle,"_",datos$Time)
    

    colors_base=c("slategray2", "orange2","darkorchid2","springgreen2")
    datos$Colorby=c(rep(colors_base[1],3),rep(colors_base[2],3),rep(colors_base[3],3),rep(colors_base[4],3))


    PC1PC2_plot=ggplot(datos)+geom_point(aes(x=PC1,y=PC2,color=Setup),size=2)+
    scale_colour_manual(values=colors_base)+xlab(paste0("PC1: ",round(100*sum_pca$PC1[2],digits=2),"% variance"))+ylab(paste0("PC2: ",round(100*sum_pca$PC2[2],digits=2),"% variance"))
    
    PC1PC3_plot=ggplot(datos)+geom_point(aes(x=PC1,y=PC3,color=Setup),size=2)+
    scale_colour_manual(values=colors_base)+xlab(paste0("PC1: ",round(100*sum_pca$PC1[2],digits=2),"% variance"))+ylab(paste0("PC3: ",round(100*sum_pca$PC3[2],digits=2),"% variance"))

    plot3d(datos$PC1,datos$PC2,datos$PC3,col=datos$Colorby,size=15,xlab=paste0("PC1: ",round(100*sum_pca$PC1[2],digits=2),"% variance"),ylab=paste0("PC2: ",round(100*sum_pca$PC2[2],digits=2),"% variance"),zlab=paste0("PC3: ",round(100*sum_pca$PC3[2],digits=2),"% variance"))

    
    
    
    ## Extract DE stats

    ParticlePEG=topTable(fit,"ParticlePEG",number=nrow(reads))
    Time_at_Au=topTable(fit,"ParticleAu:Time24h",number=nrow(reads))
    Time_at_PEG=topTable(fit,"ParticlePEG:Time24h",number=nrow(reads))
    ## completar.
    
    ## Sacamos:
    ## PCA.
    ## Venn diagrams de niveles de expresión para los 5 umbrales propuestos. Guardamos para cada Venn diagram: sets de txt files con los genes en cada región del Venn.
    ## 5 volcanos (uno por contraste) , y un barplot direccional con la cantidad de señal (a 5% y a 10% FDR, en dos paneles)+ 5 files con los estadísticos (en excel)

