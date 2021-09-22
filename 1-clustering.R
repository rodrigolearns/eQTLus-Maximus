###############################################
##### seed transcript clustering analysis #####
###############################################


#### prepare the script working environment ####
  remove(list = ls())
  gc()
  set.seed(1000)  
  
  # set working directory ####
  work.dir <- "C:/Users/harta005/Projects/seed-germination-qtl"
  setwd(work.dir)

  # dependencies ####
  library(dplyr)
  library(ggplot2)
  library(limma)
  library(Biobase)
  library(stats)
  library(gplots)
  library(RColorBrewer)
  library(Rfast)
  library(amap)
  library(topGO)
  library(org.At.tair.db)
  library(colorspace)
  # unused libraries ####
  # library(Mfuzz)
  # library(factoextra)
  # library(doParallel)
  # library(forcats)
  # library(tidyr)
  # library(gridExtra)
  # library(Hmisc)


  
# load the data ####
  trait.matrix <- read.csv(file = 'files/trait-matrix.csv', row.names = 1)
  sample.list <- read.csv(file = 'files/sample-list.csv')
  sample.stage <- sample.list$stage
  stages <- c('pd', 'ar', 'im', 'rp')
  
# data presentation ####
  presentation <- theme(axis.text.x = element_text(size=6, face="bold", color="black"),
                        axis.text.y = element_text(size=6, face="bold", color="black"),
                        axis.title.x = element_text(size=7, face="bold", color="black"),
                        axis.title.y = element_text(size=7, face="bold", color="black"),
                        strip.text.x = element_text(size=7, face="bold", color="black"),
                        strip.text.y = element_text(size=7, face="bold", color="black"),
                        strip.text = element_text(size =7, face="bold", color="black"),
                        #plot.title = element_text(size=15, face="bold"),
                        panel.background = element_rect(fill = "white",color="black"),
                        panel.grid.major = element_line(colour = "grey80"),
                        panel.grid.minor = element_blank())
  
# create the expression object ####
  x <- trait.matrix # the expression matrix
  genes <- data.frame(trait = rownames(x)) # the gene list
  f <- data.frame(trait = rownames(x)) # the gene feature
  rownames(f) <- rownames(x)
  p <- data.frame(id = colnames(trait.matrix), stage = sample.stage)
  rownames(p) <- colnames(x)
  eset <- ExpressionSet(assayData = as.matrix(x),
                        phenoData = AnnotatedDataFrame(p),
                        featureData = AnnotatedDataFrame(f))
  
# PCA ####
  pr.out <- prcomp(x = (t(trait.matrix)),  center = T, scale. = F)
  pc.df <- data.frame(pc1 = pr.out$x[, 1], 
                          pc2 = pr.out$x[, 2], 
                          stage = as.character(sample.stage),
                          population = substring(colnames(trait.matrix), 1, 3))
  pc.sum <- summary(pr.out)
  pc1.var <- round(pc.sum$importance[2, 'PC1'] * 100, 2)
  pc2.var <- round(pc.sum$importance[2, 'PC2'] * 100, 2)
  
  pca.all <- ggplot(pc.df, aes(x = pc1, y = pc2, color = factor(stage, level = c('parent', 'pd', 'ar', 'im', 'rp')))) + 
                            geom_point(aes(shape = population)) +
                            scale_colour_manual(values = c('#ccbb44', '#228833', '#4477aa', '#cc3311')) +
                            scale_shape_manual(values = c(2, 19, 6)) +
                            labs(x = paste0("PC1 (", pc1.var, "%)"), y = paste0("PC2 (", pc2.var, "%)")) +
                            labs(colour = 'stage') +
                            theme(text = element_text(size = 10),
                                  panel.background = element_blank(),
                                  panel.border=element_rect(fill=NA),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(),
                                  strip.background=element_blank(),
                                  axis.text.x=element_text(colour="black"),
                                  axis.text.y=element_text(colour="black"),
                                  axis.ticks=element_line(colour="black"),
                                  plot.margin=unit(c(1,1,1,1),"line")) +
                            #presentation +
                            geom_hline(aes(yintercept = 0), linetype = 'dashed', size = .1) +
                            geom_vline(aes(xintercept = 0), linetype = 'dashed', size = .1) 
  tiff(file = paste0("figures/pca-all.tiff"), 
       width = 2250, 
       height = 1200, 
       units = 'px',
       res = 300,
       compression = 'lzw')
  print(pca.all)
  dev.off()  
  
# PCA per stage ####
  for (i in stages) {

    
    trait.stage <- trait.matrix[, which(sample.stage == i)]
    pr.stage <- prcomp(x = (t(trait.stage)),  center = T, scale. = F)
    pc.stage <- data.frame(pc1 = pr.stage$x[, 1], 
                        pc2 = pr.stage$x[, 2], 
                        population = c(rep('parent', 4), rep('RIL', ncol(trait.stage) - 4)))
    pc.sum <- summary(pr.stage)
    pc1.var <- round(pc.sum$importance[2, 'PC1'] * 100, 2)
    pc2.var <- round(pc.sum$importance[2, 'PC2'] * 100, 2)
    pca.stage.plot <- ggplot(pc.stage, aes(x = pc1, y = pc2)) + 
      geom_point(aes(shape = population)) +
      #scale_colour_manual(values = c('#ccbb44', '#228833', '#4477aa', '#cc3311')) +
      scale_shape_manual(values = c(17, 19)) +
      labs(x = paste0("PC1 (", pc1.var, "%)"), y = paste0("PC2 (", pc2.var, "%)")) +
      labs(colour = 'stage') +
      theme(text = element_text(size = 10),
            panel.background = element_blank(),
            panel.border=element_rect(fill=NA),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background=element_blank(),
            axis.text.x=element_text(colour="black"),
            axis.text.y=element_text(colour="black"),
            axis.ticks=element_line(colour="black"),
            plot.margin=unit(c(1,1,1,1),"line")) +
      #presentation +
      geom_hline(aes(yintercept = 0), linetype = 'dashed', size = .1) +
      geom_vline(aes(xintercept = 0), linetype = 'dashed', size = .1) 

    tiff(file = paste0("figures/pca-", i, ".tiff"), 
         width = 2250, 
         height = 1200, 
         units = 'px',
         res = 300,
         compression = 'lzw')
    print(pca.stage.plot)
    dev.off()  
  }
  
  

  # differential expresed gene analysis and hierarchiecal clustering ####
  group <- with (pData(eset), stage)
  group <- factor(group)
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  colSums(design) 
  
  # create pairwise contrast matrix among four different stages
  cm <- makeContrasts(pd_ar = pd - ar,
                      ar_im = ar - im,
                      im_rp = im - rp,
                      levels = design)
  coef <- colnames(cm)
  
  # fit the linear model to determine significant differentially expressed genes
  # for each pairwise contrast 
  fit <- lmFit(eset, design)
  fit2 <- contrasts.fit(fit, contrasts = cm) %>%
    eBayes()
  result <- decideTests(fit2)
  summary(result) # 1 = upregulate
  write.csv(result, 'files/differentally-expressed-genes.csv')

  # determine the differentially expressed genes for each pairwise contrast  
  # based on p value (>0.05) and  log fold change more than 1 sorted by fold change
  de.genes <- NA
  coef2 <- 
  for (i in 1:length(coef)) {
    print(i)
    deg.tmp <- topTable(fit2, lfc = 1, coef = coef[i], 
                        p.value = 0.05, 
                        sort.by = 'logFC', 
                        number = 100000)
    de.genes <- c(de.genes, as.character(deg.tmp$trait))
  }
  
  de.genes <- unique(de.genes) # as much as 990 genes are differentially expressed between one of the contrasts
  cluster.data <- trait.matrix[which(rownames(trait.matrix) %in% de.genes), ]
  colnames(cluster.data) <- sample.stage
  
# hierarchiecal tree clustering done for stage and genes ####
  
  # hc for stage
  dist.matrix.stage <- amap::Dist(t(cluster.data), 
                            method = 'pearson', 
                            upper = T, 
                            diag = T) 
  hc.stage <- hclust(dist.matrix.stage, method = 'ward.D2')
  #hc.stage$order <- c(hc.stage$order[46:180], hc.stage$order[1:45])
  reorder(x = as.dendrogram(hc.stage), c(46:180, 1:45))
  plot(reorder(x = as.dendrogram(hc.stage), 180:1, agglo.FUN = mean))
  
  # hc for gene
  dist.matrix.gene <- amap::Dist(cluster.data, 
                           method = 'pearson', # to group the genes based on patterns similarity across stages
                           upper = T, 
                           diag = T)  
  hc.gene <- hclust(dist.matrix.gene, method = 'ward.D2') # similar to average linkage-
  # it calculates the distance that is minimizing the variance within cluster
  # and maximazing the variance between clusters
  plot(hc.gene)
  hc.gene <- cutree(hc.gene, k =6) # if k > 5, the cluster will be diproporsional i.e. a cluster only have few member
  table(hc.gene)
  hc.table <- as.data.frame(table(true = rownames(cluster.data), cluster = hc.gene))
  hc.table <- hc.table[which(hc.table$Freq != 0), ]
  hc.table$cluster <- as.numeric(hc.table$cluster)
  hc.table$true <- as.character(hc.table$true)
  hc.table$Freq <- NULL
  table(hc.table$cluster)
  write.csv(hc.table, 'files/gene-cluster.csv')
  
  # heatmap and hierarchiecal clustering ####
  sample.color <- ifelse(grepl('pd', colnames(trait.matrix)), '#4daf4a', 
                         ifelse(grepl('ar', colnames(trait.matrix)), '#377eb8', 
                                ifelse(grepl('im', colnames(trait.matrix)), '#984ea3',
                                       ifelse(grepl('rp', colnames(trait.matrix)), '#e41a1c', NA))))
  
  tiff(file = paste0("figures/heatmap-genes.tiff"), 
       width = 2250, 
       height = 2000, 
       units = 'px',
       res = 300,
       compression = 'lzw')
  heatmap.2(x = as.matrix(cluster.data), cexRow = 0.8, cexCol = 1.2,
                  distfun = function(x) amap::Dist(x, method = 'pearson'),
                  hclust = function(x) hclust(x, method = 'ward.D2'),
                  scale = 'row', 
                  density.info = "density", 
                  trace = 'none', 
                  col = diverging_hcl(100, palette = 'Blue-Red3'),
                  keysize = 1,
                  key.title = 'Z-score of\nlog-intensities',
                  key.xlab = NA, 
                  Colv = reorder(x = as.dendrogram(hc.stage), 180:1, agglo.FUN = mean))
  dev.off() 
 
# GOE analysis for genes in each cluster ####
  x <- org.At.tairCHR
  all.genes <- as.list(rownames(trait.matrix))
  hc.go.list <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 8))
  colnames(hc.go.list) <- c('GO.ID', 'Term', 'Annotated', 'Significant', 'Expected', 'Fisher', 
                            'FDR', 'cluster')
  ontology <- 'BP'
  
  for (i in unique(hc.table$cluster)) {
    gene.set <- hc.table$true[which(hc.table$cluster == i)]
    gene.set <- factor(as.integer(all.genes %in% gene.set))
    names(gene.set) <- all.genes
    GOdata <- new("topGOdata",
                  description = "Analyzing clustering results", ontology = ontology,
                  allGenes = gene.set, 
                  annot = annFUN.org,mapping= "org.At.tair.db")
    resultFisher <- runTest(GOdata, algorithm = 'weight', statistic = "fisher")
    # GOE from topgo consider the general terms on the hierarchy
    # the weight algortihm maintain the balance between type I and II errors
    result.df <- GenTable(GOdata, Fisher = resultFisher,
                          orderBy = "Fisher", ranksOf = "Fisher", topNodes = length(resultFisher@score))
    result.df$FDR <- p.adjust(p = result.df$Fisher, method = 'fdr')
    result.df$cluster <- i
    result.df <- result.df[order(result.df$FDR), ]
    result.df <- filter(result.df, FDR <= 0.001)
    hc.go.list <- rbind(hc.go.list, result.df)
  }
  
  write.csv(hc.go.list, paste0('files/goe-tables', ontology, '.csv'))

  