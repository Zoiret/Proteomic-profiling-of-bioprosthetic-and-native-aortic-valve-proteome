###This is a code to reproduce statistical analysis of proteomics data from Proteolytic degradation is the culprit behind a bioprosthetic heart valve failure

##Openinig_the_data
#We recommend you to set the working directory for easy reproduction  of the code
#setwd("your directory")

rm(list = ls())

{
  library(BiocManager)
  library(readxl)
  library(devtools)
  library(ggvenn)
  library(impute)
  library(RColorBrewer)
  library(limma)
  library(vsn)
  library(NMF)
  library(ggplot2)
  library(mixOmics)
  library(EnhancedVolcano)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(pheatmap)
  library(openxlsx)
}

### Functions

Tiff <- function(name = '123name.tiff',
                 task,
                 width=8,
                 height=8,
                 res=900) {
  tiff(name, 
       units="in", 
       width = width,
       height = height,
       res = res, 
       compression = 'lzw')
  task
  dev.off()
}

Graph_CC <- function(dataset,legend = colnames(dataset)) {
  f <- 1:length(dataset)
  for (i in 1:length(dataset)) {f[i] <- mean(complete.cases(dataset[,i]))}
  plot.ts(f,x = seq(1, length(dataset), 1), 
          xy.labels = legend, 
          ylab = "Complete cases", 
          xlab = "Samples", axes = F) 
  axis(1, at = 0:length(dataset))
  axis(2, at = seq(0,1,0.05))
}
#Uploading and preparing data

dat <- data.frame(read.csv("protein_data.csv"))

  dat[dat==0] <- NA 
  dat$UNIPROT <- sub("\\|.*", "", dat$Accession)
  dat$Protein.name <- sub(".*\\|", "", dat$Accession)
  dat$Gene.name <- mapIds(org.Hs.eg.db, keys = dat$UNIPROT, column = "SYMBOL", keytype = "UNIPROT")
  
  dat1 <- dat[,c(21:35)]
  dat1 <- dat1[,c(1,3,9,10,11,
                  4,6,12,13,15,
                  2,5,7,8,14)]
  head(dat1)
  
  str(dat1)
  rownames(dat1) <- paste(dat$UNIPROT, dat$Gene.name, sep="--")
  colnames(dat1) <- c('Native_1','Native_2','Native_3','Native_4','Native_5',
                      'Bio_1_1','Bio_1_2','Bio_1_3','Bio_1_4','Bio_1_5',
                      'Bio_2_1','Bio_2_2','Bio_2_3','Bio_2_4','Bio_2_5')
  head(dat1)
  
fact <- data.frame(read_excel("sample_info.xlsx")) 
  
  fact$Name <- fact[c(1,3,9,10,11,
                      4,6,12,13,15,
                      2,5,7,8,14),1]

  fact$Type <- factor(fact$Type, levels = c("Native", "Bio_1", "Bio_2"))
  fact$Type <- fact[c(1,3,9,10,11,
                      4,6,12,13,15,
                      2,5,7,8,14),2]
  table(fact$Type)
  
  
  fact$Type2 <- fact[c(1,3,9,10,11,
                       4,6,12,13,15,
                       2,5,7,8,14),3]
  fact$Type2 <- factor(fact$Type2, levels = c("Native", "Bio"))
  table(fact$Type2)
  
  fact$Type3 <- c('Native_1','Native_2','Native_3','Native_4','Native_5',
                  'Bio_1_1','Bio_1_2','Bio_1_3','Bio_1_4','Bio_1_5',
                  'Bio_2_1','Bio_2_2','Bio_2_3','Bio_2_4','Bio_2_5')
  fact$Type3 <- factor(fact$Type3, levels = c('Native_1','Native_2','Native_3','Native_4','Native_5',
                                              'Bio_1_1','Bio_1_2','Bio_1_3','Bio_1_4','Bio_1_5',
                                              'Bio_2_1','Bio_2_2','Bio_2_3','Bio_2_4','Bio_2_5'))
  
  rownames(fact) <- colnames(dat1)
  fact <- fact[,-c(1,3)]
  
  str(fact)

Graph_CC(dat1)

## Qualititative analysis

#detection of optimal way is to filter out proteins with too much NA
{  
  df <- data.frame(part=0:1000)
  c <- matrix(,ncol = 1001, nrow=3)
  for (i in 0:1000) 
  {
    Nat <- dat1[which(rowMeans(!is.na(dat1[,1:5])) >= i/1000), ]
    Bio_1 <- dat1[which(rowMeans(!is.na(dat1[,6:10])) >= i/1000), ]
    Bio_2 <- dat1[which(rowMeans(!is.na(dat1[,11:15])) >= i/1000), ]
    c[1,i+1] <- nrow(Nat)
    c[2,i+1] <- nrow(Bio_1)
    c[3,i+1] <- nrow(Bio_2)
  }
  df['Nat'] <- c[1,]
  df['Bio_1'] <- c[2,]
  df['Bio_2'] <- c[3,]
  ggplot(data = df) +
    geom_line(mapping = aes(x = part, y = Nat, colour = "Nat")) +
    geom_line(mapping = aes(x = part, y = Bio_1, colour = "Bio_1")) +
    geom_line(mapping = aes(x = part, y = Bio_2, colour = "Bio_2")) +
    geom_vline(xintercept = 600) +
    xlab("Decline NA (‰)") + ylab("Detected proteins") +
    scale_x_continuous(breaks = seq(0, 1000, 100)) +
    scale_y_continuous(breaks = seq(0, 2000, 200)) +
    labs(subtitle="Dependence of the number of detectable proteins on decline NA")
}

{  
  df <- data.frame(part=0:1000)
  c <- matrix(,ncol = 1001, nrow=2)
  for (i in 0:1000) 
  {
    Nat <- dat1[which(rowMeans(!is.na(dat1[,1:5])) >= i/1000), ]
    Bio <- dat1[which(rowMeans(!is.na(dat1[,6:15])) >= i/1000), ]
    c[1,i+1] <- nrow(Nat)
    c[2,i+1] <- nrow(Bio)
  }
  df['Nat'] <- c[1,]
  df['Bio'] <- c[2,]
  ggplot(data = df) +
    geom_line(mapping = aes(x = part, y = Nat, colour = "Nat")) +
    geom_line(mapping = aes(x = part, y = Bio, colour = "Bio")) +
    geom_vline(xintercept = 600) +
    xlab("Decline NA (‰)") + ylab("Detected proteins") +
    scale_x_continuous(breaks = seq(0, 1000, 100)) +
    scale_y_continuous(breaks = seq(0, 2000, 200)) +
    labs(subtitle="Dependence of the number of detectable proteins on decline NA")
}


#Extraction group-specific proteins
  Nat <- dat1[which(rowMeans(!is.na(dat1[,c(1:5)])) >= 0.6), ]
  Bio_1 <- dat1[which(rowMeans(!is.na(dat1[,c(6:10)])) >= 0.6), ]
  Bio_2 <- dat1[which(rowMeans(!is.na(dat1[,c(11:15)])) >= 0.6), ]

#Venn diagram
  vennn <- list(Native = rownames(Nat), Bio_1 = rownames(Bio_1), Bio_2 = rownames(Bio_2))
ggvenn(vennn, 
       fill_color = c("#0073C2FF", "#CD534CFF", 'green'),
       stroke_size = 0.5, set_name_size = 8, text_size = 7,)

  # Multiple set version of intersect
  Intersect <- function (x) {
    # x is a list
    if (length(x) == 1) {
      unlist(x)
    } else if (length(x) == 2) {
      intersect(x[[1]], x[[2]])
    } else if (length(x) > 2){
      intersect(x[[1]], Intersect(x[-1]))
    }
  }
  
  # Multiple set version of union
  Union <- function (x) {  
    # x is a list
    if (length(x) == 1) {
      unlist(x)
    } else if (length(x) == 2) {
      union(x[[1]], x[[2]])
    } else if (length(x) > 2) {
      union(x[[1]], Union(x[-1]))
    }
  }
  
  # Remove the union of the y's from the common x's 
  Setdiff <- function (x, y) {
    # x and y are lists of characters
    xx <- Intersect(x)
    yy <- Union(y)
    setdiff(xx, yy)
  }
  
  
  Nat_spec <- Setdiff(vennn[c("Native")], vennn[c('Bio_1','Bio_2')])
  Bio_1_spec <- Setdiff(vennn[c("Bio_1")], vennn[c('Native','Bio_2')])
  Bio_2_spec <- Setdiff(vennn[c("Bio_2")], vennn[c('Native','Bio_1')])


Datafile_1 <- list('Native HV' = Nat_spec,
              'Bioprosthetic HV 1' = Bio_1_spec,
              'Bioprosthetic HV 2' = Bio_2_spec)
# write.xlsx(Datafile_1, 'Supplementary Datafile 1 (Specific proteins for NHV, BHV1, BHV2).xlsx', rownames = T)
  

## Quantitative analysis

    #detection of optimal way is to filter out proteins with too much NA
  {  
    df <- data.frame(part=0:1000)
    c <- matrix(,ncol = 1001, nrow=2)
    for (i in 0:1000) 
    {
      completcases <- dat1[which(rowMeans(!is.na(dat1)) >= i/1000),]
      c[1,i+1] <- mean(complete.cases(completcases))
      c[2,i+1] <- length(rownames(completcases))
    }
    df['completcases'] <- c[1,]
    df['proteins'] <- c[2,]
    ggplot(data = df, aes(x = part)) +
      geom_line(mapping = aes(y = completcases, colour = "completcases")) + 
      geom_smooth(mapping = aes(y = completcases, colour = "completcases")) +
      geom_line(mapping = aes(y = proteins/c[2,2], colour = "proteins")) + 
      geom_smooth(mapping = aes(y = proteins/c[2,2], colour = "proteins")) +
      geom_hline(yintercept = 0.5) +
      xlab("Decline NA (‰)") + ylab("Complete cases") +
      scale_x_continuous(breaks = seq(0, 1000, 100)) +
      scale_y_continuous(breaks = seq(0, 1, 0.1), sec.axis = sec_axis(~ . * c[2,2], name = "Number of proteins", breaks = seq(0, 2000, 200))) +
      labs(subtitle="Dependence of the number of detectable proteins on decline NA")
  }

#Removing rows with a lot of missing values
  dat2 <- dat1[which(rowMeans(!is.na(dat1)) >= 0.75), ]
  dat2 <- rbind(dat2,dat1[c('P39900--MMP12','P08253--MMP2','Q99542--MMP19','P14780--MMP9','P22894--MMP8',
                            "P01033--TIMP1",'P16035--TIMP2'),])
  mean(complete.cases(dat2))
  NAsums <- data.frame(colSums(is.na(dat2))) #detection NA in experimental sample
  NAsums 
  str(dat2)

#knn imputation of missng values
  tdat <- t(dat2)
  dat_knn1 <- impute.knn(tdat, k = 5)
  dat_knn <- t(dat_knn1$data)
  mean(complete.cases(dat_knn))

#Normalization and data QC
  pal <- brewer.pal(n = 9, name = "Set1")
  cols <- pal[fact$Type]
boxplot(dat_knn, outline = FALSE, col = cols, main = "Raw data")
  legend("topright", levels(fact$Type), fill = pal, bty = "n", xpd = T)
  data.frame(colSums(dat_knn))
  
#Logarithm of data
  dat_log <- log2(dat_knn+1)
  head(dat_log)
  
  mean(complete.cases(dat_log))
boxplot(dat_log, outline = FALSE, col = cols, main = "Log-transformed data")
  legend("topright", levels(fact$Type), fill = pal, bty = "n", xpd = T)

#Quantile normalization
  dat_norm <- normalizeQuantiles(dat_log)
  head(dat_norm)
  
boxplot(dat_norm, col = cols, main = "Normalized data")
  legend("topright", levels(fact$Type), fill = pal, bty = "n", xpd = T)
  
  mean(complete.cases(dat_norm))
  colSums(is.na(dat_norm))

#MAplot (Log-expression)
maplot <- function(X1, X2, pch = 21, main = "MA-plot", xlab = "Average log-expression", ylab = "Expression log-ratio", lpars = list(col = "blue", lwd = 2), ...){
    X <- (rowMeans(X2) + rowMeans(X1)) / 2
    Y <- rowMeans(X2) - rowMeans(X1)
    scatter.smooth(x = X, y = Y,
                   main = main, pch = pch,
                   xlab = xlab, ylab = ylab,
                   lpars = lpars, ...)
    abline(h = c(-1, 0, 1), lty = c(2, 1, 2))
  }
  
maplot(dat_log[, rownames(fact)[fact$Type == "Native"]], dat_log[, rownames(fact)[fact$Type == "Bio_1"]])
maplot(dat_log[, rownames(fact)[fact$Type == "Native"]], dat_log[, rownames(fact)[fact$Type == "Bio_2"]])
maplot(dat_log[, rownames(fact)[fact$Type == "Bio_1"]], dat_log[, rownames(fact)[fact$Type == "Bio_2"]])

#MeanSd
meanSdPlot(as.matrix(dat_log))
meanSdPlot(as.matrix(dat_norm))

#Heatmap
aheatmap(cor(dat_norm), 
         color = "-RdBu:256", 
         # annCol = colnames(dat_norm), 
         fontsize = 10)

#Principle components analysis
  dat_pca <- pca(t(dat_norm), ncomp = 6, center = TRUE)
plot(dat_pca)

plotIndiv(dat_pca, 
          comp = c(1,2),
          ind.names = F, 
          group = fact$Type, 
          legend = TRUE, 
          ellipse = TRUE,
          title = '', 
          ellipse.level = 0.95,
          pch = c(2,1,5),
          col.per.group = c( "green","#0073C2FF", '#cd00cd')) 

#sPLS-DA clusterization
  t_dat1 <- t(dat_norm)

  ordination.optimum.splsda <- splsda(t_dat1, fact$Type, ncomp = 3, keepX = c(15,15,15))
plotIndiv(ordination.optimum.splsda, 
          ind.names = F, 
          ellipse = T, 
          title = "PLS-DA ordination of different donors", 
          legend=TRUE,
          pch = c(15,16, 17),
          col.per.group = c('#ff8800',"#0073C2FF",'#cd00cd'))

# Tiff('PLS-DA.tiff', plotIndiv(ordination.optimum.splsda,ind.names = F,ellipse = T,title = "PLS-DA ordination of different donors",legend=TRUE,pch = c(15,16, 17),col.per.group = c('#ff8800',"#0073C2FF",'#cd00cd')))
dev.off()



#Limma - differentially expressed proteins

  X <- model.matrix(~ fact$Type2)
  X
  
  fit <- lmFit(dat_norm, design = X, method = "robust", maxit = 10000)
  
  efit <- eBayes(fit)
  
  topTable(efit, coef = 2)
  numGenes <- length(dat_norm)
  full_list_efit <- topTable(efit, number = length(dat_norm[,1]))

  full_list_efit$UNIPROT <- sub("\\--.*", "", rownames(full_list_efit))
  full_list_efit<- merge(full_list_efit, dat[,c(56:58)], by = c('UNIPROT'), all.x = T, all.y=F)
  full_list_efit <- full_list_efit[!duplicated(full_list_efit$Gene.name),]
  
  rownames(full_list_efit) <- paste(full_list_efit$UNIPROT, full_list_efit$Gene.name, sep = '--')
  full_list_efit <- full_list_efit[-1,]
  full_list_efit <- full_list_efit[,c(1,9,8,2,3,6)]
  
  
  Datafile_2 <- list("NHV vs BHV" = full_list_efit)
# write.xlsx(Datafile_2, 'Supplementary Datafile 2 (Different proteins for NHV and BHV).xlsx', rownames = T)
  
EnhancedVolcano(full_list_efit,
                lab = full_list_efit$Gene.name,
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                FCcutoff = 2,
                xlim = c(-6, 10),
                ylim = c(0, 10),
                title ="Native versus Bio",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)

EnhancedVolcano(full_list_efit,
                lab = NA,
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                FCcutoff = 2,
                xlim = c(-6, 10),
                ylim = c(0, 10),
                title ="Native versus Bio",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)

EnhancedVolcano(full_list_efit,
                lab = full_list_efit$Gene.name,
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                FCcutoff = 1,
                xlim = c(-6, 10),
                ylim = c(0, 10),
                title ="Native versus Bio",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)

EnhancedVolcano(full_list_efit,
                lab = NA,
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                FCcutoff = 1,
                xlim = c(-6, 10),
                ylim = c(0, 10),
                title ="Native versus Bio",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)
  
#diff expression Bio_1_bio2
  fact_bio <- fact[6:15,]
  dat_bio <- dat_norm[,6:15]
  
  fact_bio$Type <- as.factor(as.character(fact_bio$Type))
  fact_bio$Type
  X_N2 <- model.matrix(~ fact_bio$Type)
  X_N2
  
  fit_N2 <- lmFit(dat_bio, design = X_N2, method = "robust", maxit = 10000)
  
  # Empirical Bayes statistics
  efit_N2 <- eBayes(fit_N2)
  
  # Dif_expr_table
  topTable(efit_N2, coef = 2)
  full_list_efit_N2 <- topTable(efit_N2, number = length(dat_bio))
  #write.csv(full_list_efit_N2,'Dif_expr_Native_vs_bio2.csv')
  head(full_list_efit_N2)
  str(full_list_efit_N2)

EnhancedVolcano(full_list_efit_N2,
                lab = rownames(full_list_efit_N2),
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                xlim = c(-5,5),
                ylim = c(0, 3),
                FCcutoff = 2,
                title ="Bio1 versus Bio2",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)

EnhancedVolcano(full_list_efit_N2,
                lab = NA,
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                xlim = c(-5,5),
                ylim = c(0, 3),
                FCcutoff = 2,
                title ="Bio1 versus Bio2",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)

EnhancedVolcano(full_list_efit_N2,
                lab = rownames(full_list_efit_N2),
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                xlim = c(-5,5),
                ylim = c(0, 3),
                FCcutoff = 1,
                title ="Bio1 versus Bio2",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)

EnhancedVolcano(full_list_efit_N2,
                lab = NA,
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                xlim = c(-5,5),
                ylim = c(0, 3),
                FCcutoff = 1,
                title ="Bio1 versus Bio2",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)
  
#Heatmap

  prot_for_heatmap_up <- rownames(full_list_efit[c('P01024','P0C0L5','P0C0L4','P10643','P07357','P07358','P07360','P02748','Q9BXR6',
                                                'P27918','P10909','P09871','P13671','Q03591','P36980','P62805','P04003',
                                                'P01903','P16401','Q71DI3','P00742','Q96IY4','P00740','P08709','Q9UK55','P04275',
                                                'P80511','P05109','P06702','P02788','P31151','P37840','P48061','P59665','Q92954',
                                                'P04004','O43866','P02649','P02042','P04040','Q06033','P48740','P23142','P33908','P08493','P07237','P08195'),])
  length(prot_for_heatmap_up)
  
  prot_for_heatmap_down <- rownames(full_list_efit[c('P18669','P14618','P36871','P12109','P00558','P60174','P09972','O95336',
                                                     'P15121','O43488','Q04760','Q14974','P14174','O75874','P09211',
                                                     'P06396','P20810','O60664','Q16658','Q99497','Q92597',
                                                     'Q09666','P08758','P04632','P19827','Q07954','P35555','P13611','P26447'),])
  length(prot_for_heatmap_down)
  
  prot_for_heatmap <- rownames(full_list_efit[c('P01024','P0C0L5','P0C0L4','P10643','P07357','P07358','P07360','P02748','Q9BXR6',
                                                   'P27918','P10909','P09871','P13671','Q03591','P36980','P62805','P04003',
                                                   'P01903','P16401','Q71DI3','P00742','Q96IY4','P00740','P08709','Q9UK55','P04275',
                                                   'P80511','P05109','P06702','P02788','P31151','P37840','P48061','P59665','Q92954',
                                                   'P04004','O43866','P02649','P02042','P04040','Q06033','P48740','P23142','P33908','P08493','P07237','P08195',
                                                'P18669','P14618','P36871','P12109','P00558','P60174','P09972','O95336',
                                                'P15121','O43488','Q04760','Q14974','P14174','O75874','P09211',
                                                'P06396','P20810','O60664','Q16658','Q99497','Q92597',
                                                'Q09666','P08758','P04632','P19827','Q07954','P35555','P13611','P26447'),])
  length(prot_for_heatmap)
  
pheatmap(dat_norm[prot_for_heatmap_up,], #main heatmap
         main = 'Heatmap upregulated proteins',
         annotation_col = fact,
         cluster_cols = F,
         cluster_rows = T,
         cellwidth = 12,
         cellheight = 10,
         border_color = "black",
         color = colorRampPalette(c("#00bfff", '#005aeb', '#240935', '#800080','#cd00cd',"#ff19ff"))(100))

pheatmap(dat_norm[prot_for_heatmap_down,], #main heatmap
         main = 'Heatmap downregulated proteins',
         annotation_col = fact,
         cluster_cols = F,
         cluster_rows = T,
         cellwidth = 12,
         cellheight = 10,
         border_color = "black",
         color = colorRampPalette(c("#00bfff", '#005aeb', '#240935', '#800080','#cd00cd',"#ff19ff"))(100))

pheatmap(dat_norm[prot_for_heatmap,], #main heatmap
         main = 'Heatmap',
         annotation_col = fact,
         cluster_cols = F,
         cluster_rows = T,
         cellwidth = 12,
         cellheight = 10,
         border_color = "black",
         color = colorRampPalette(c("#00bfff", '#005aeb', '#240935', '#800080','#cd00cd',"#ff19ff"))(100))
  

##Unique proteins and preparing table
 
#Areas and Unique proteins

  Uniq_prot <- dat[,c(21:35,56,58,57,36)]
  Uniq_prot <- Uniq_prot[,c(16,17,18,19,1,3,9,10,11,
                            4,6,12,13,15,
                            2,5,7,8,14)]
  rownames(Uniq_prot) <- paste(Uniq_prot$UNIPROT, Uniq_prot$'Gene name', sep="--")
  
  Uniq_prot_filtered <- Uniq_prot[which((rowMeans(!is.na(Uniq_prot[,5:19])) >0.75)),] 

  colnames(Uniq_prot) <- c('UNIPROT',
                           'Gene name',
                           'Protein',
                           'Peptides',
                           'Native_1','Native_2','Native_3','Native_4','Native_5',
                           'Bio_1_1','Bio_1_2','Bio_1_3','Bio_1_4','Bio_1_5',
                           'Bio_2_1','Bio_2_2','Bio_2_3','Bio_2_4','Bio_2_5')
  
  colnames(Uniq_prot_filtered) <- c('UNIPROT',
                                    'Gene name',
                                    'Protein',
                                    'Peptides',
                                    'Native_1','Native_2','Native_3','Native_4','Native_5',
                                    'Bio_1_1','Bio_1_2','Bio_1_3','Bio_1_4','Bio_1_5',
                                    'Bio_2_1','Bio_2_2','Bio_2_3','Bio_2_4','Bio_2_5')
  
  Uniq_prot_NvsB_maxCons <- Uniq_prot[(rowMeans(!is.na(Uniq_prot[,5:9])) >= 5/5) & 
                                        (rowMeans(!is.na(Uniq_prot[,10:19])) <= 0/10),]
  
  Uniq_prot_NvsB_medCons <- Uniq_prot[(rowMeans(!is.na(Uniq_prot[,5:9])) >= 4/5) & 
                                        (rowMeans(!is.na(Uniq_prot[,10:19])) <= 2/10),]
  
  Uniq_prot_NvsB_minCons <- Uniq_prot[(rowMeans(!is.na(Uniq_prot[,5:9])) >= 3/5) & 
                                        (rowMeans(!is.na(Uniq_prot[,10:19])) <= 4/10),]
  
  Uniq_prot_BvsN_maxCons <- Uniq_prot[(rowMeans(!is.na(Uniq_prot[,5:9])) <= 0/5) &
                                        (rowMeans(!is.na(Uniq_prot[,10:19])) >= 10/10),]
  
  Uniq_prot_BvsN_medCons <- Uniq_prot[(rowMeans(!is.na(Uniq_prot[,5:9])) <= 1/5) &
                                        (rowMeans(!is.na(Uniq_prot[,10:19])) >= 8/10),]
  
  Uniq_prot_BvsN_minCons <- Uniq_prot[(rowMeans(!is.na(Uniq_prot[,5:9])) <= 2/5) &
                                        (rowMeans(!is.na(Uniq_prot[,10:19])) >= 6/10),]

  
  Norm_data <- data.frame(dat_norm)
  Norm_data$UNIPROT <- sub("\\--.*", "", rownames(Norm_data))
  Norm_data <- merge(Norm_data, Uniq_prot[,1:3], by = c('UNIPROT'), all.x = T, all.y=F)
  Norm_data <- Norm_data[,c(1,17,18,2:16)]

  Datafile_3 <- list("Raw data" = dat,
                     "1614 Raw areas" = Uniq_pro,
                     "546 Filtered areas" = Uniq_prot_filtere,
                     "551 Norm areas (546 + 5 MMPs)" = Norm_data,
                     "2 Uniq prot NHV_5-5 BHV_0-10" = Uniq_prot_NvsB_maxCons,
                     "90 Uniq prot NHV_4-5 BHV_2-10" = Uniq_prot_NvsB_medCons,
                     "298 Uniq prot NHV_3-5 BHV_4-10" = Uniq_prot_NvsB_minCons,
                     "7 Uniq prot BHV_10-10 NHV_0-5" = Uniq_prot_BvsN_maxCons,
                     "53 Uniq prot BHV_8-10 NHV_1-5" = Uniq_prot_BvsN_medCons,
                     "129 Uniq prot BHV_6-10 NHV_2-5" = Uniq_prot_BvsN_minCons)
  
# write.xlsx(Datafile_3, 'Supplementary Datafile 3 (Proteins NHV and BHV).xlsx', rownames = T)
  