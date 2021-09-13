##############################
##### seed eQTL analysis #####
##############################

#### prepare the script working environment ####
  remove(list = ls())
  gc()
  set.seed(1000)  
  
  # Set working directory ####
  getwd()
  work.dir <- "C:/Users/Rodrigo/Desktop/MiPRO/1_MASTER_DATA_MIPRO/2_Scripts/111 MyScripts/1A_Hartanto_Edited/"
  setwd(work.dir)

  # dependencies ####
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
   BiocManager::install("Mfuzz")
   BiocManager::install("topGO")
   BiocManager::install("org.At.tair.db")
   BiocManager::install("Rgraphviz")
   BiocManager::install("org.At.eg.db")
  alibrary(doParallel)
  library(dplyr)
  library(ggplot2)
  library(heritability)
  library(topGO)
  library(org.At.tair.db)

  
# load required function ####
  setwd('functions/')
  for(i in 1:length(dir())){
    source(dir()[i])
  } # read function from Mark
  setwd(work.dir)
  
  write.EleQTL <- function(map1.output,filename){
    
    selector <- cbind(trait = rownames(map1.output$LOD), pval = apply(map1.output$LOD,1,max,na.rm=T)) %>%
      data.frame()
    print('1')
    rownames(selector) <- NULL
    
    lod <- map1.output$LOD
    lod <- lod[rownames(lod) %in% selector[,1],]
    rownames(lod) <- selector$trait
    colnames(lod) <- map1.output$Marker[,1]
    print('2')
    eff <- map1.output$Effect
    eff <- eff[rownames(eff) %in% selector[,1],]
    rownames(eff) <- selector$trait
    colnames(eff) <- map1.output$Marker[,1]
    
    lod.eff <- lod*sign(eff)
    print('3')
    dat <- map1.output$Trait
    dat <- dat[rownames(dat) %in% selector[,1],]
    rownames(dat) <- selector$trait
    colnames(dat) <- colnames(map1.output$Map)
    print('4')
    map <- map1.output$Map
    rownames(map) <- map1.output$Marker[,1]
    
    marker <- map1.output$Marker
    print('5')
    help("connection")
    write.table(lod,file=paste(filename,"_lod.txt",sep=""),sep="\t",quote=F)
    write.table(eff,file=paste(filename,"_eff.txt",sep=""),sep="\t",quote=F)
    write.table(lod.eff,file=paste(filename,"_lodxeff.txt",sep=""),sep="\t",quote=F)
    write.table(marker,file=paste(filename,"_marker.txt",sep=""),sep="\t",quote=F)
    write.table(dat,file=paste(filename,"_data.txt",sep=""),sep="\t",quote=F)
    write.table(map,file=paste(filename,"_map.txt",sep=""),sep="\t",quote=F)    
    print('6')
  } # function to convert mapping result to tables
  map.per.marker <- function(trait, marker) {
    model <- lm(terms(trait ~ marker, keep.order = FALSE))
    summ <- summary(model)
    pval <- summ$coefficients[2, 4]
    lod <- -log10(pval)
    eff <- summ$coefficients[2, 1]
    output <- c(lod, eff)
    return(output)
  } # map the QTL at a marker location
  map.all.marker <- function(trait, markers) {
    eff.out <- rep(NA, nrow(markers))
    pval.out <- rep(NA, nrow(markers))
    for (i in 1:nrow(markers)) {
      if(i == 1) {
        out.tmp <- map.per.marker(trait, markers[i, ])
      }
      if( i != 1 & sum(abs(as.numeric(markers[i-1,])  - as.numeric(markers[i,])), na.rm = T) != 0 ) {
        out.tmp <- map.per.marker(trait, markers[i, ])
      }
      if( i != 1 & sum(abs(as.numeric(markers[i-1,]) - as.numeric(markers[i,])), na.rm = T) == 0 ) {
        out.tmp <- out.tmp
      }
      pval.out[i] <- out.tmp[1]
      eff.out[i] <- out.tmp[2]
      output.lod <- cbind(pval.out, eff.out)
      colnames(output.lod) <- c('LOD', 'Eff')
    }
    return(output.lod)
  } # map the QTL using genome wide markers
  threshold.determination <- function(trait, strain.map, n.perm){
    
    traits <- t(replicate(n.perm, trait))
    perm.trait <- permutate.traits(traits)
    
    ###Check for NAs
    pval.out <- matrix(NA,nrow(perm.trait),nrow(strain.map))
    
    for (i in 1:nrow(perm.trait)) {
      pval.out[i, ] <- fast.lod.all.marker(perm.trait[i, ], strain.map)
    }
    
    pval.distribution <- apply(pval.out, 1, max)
    threshold <- quantile(x = pval.distribution, probs = 0.95)
    return(threshold)
  } # determine threshold for a QTL


# Loop for sequential dataset analysis ------------------------------------
  #Load requiered datasets
    #Fill-in the names of all the studies
      studies <- c('Keurentjes','West','Cubillos','Snoek','Lowry','Imprialou','Rabanal') #'Hartanto',
    #Make a variable with the datasets inside each study
      Keurentjes <- c('keure') #alreadt run
      West <- c('wesw1','wesw2')
      Cubillos <- c('cub3b','cub3c','cub5b','cub5c')
      Snoek <- c('snoek')
      Lowry <- c('lowdr','lowwe')
      Hartanto <- c('harar','harim','harpd','harrp')
      Imprialou <- c('impri')
      Rabanal <- c('raban')
      
    #Placing the studies with their datasets within a list
      datasets <- list(Keurentjes,West,Cubillos,Snoek,Lowry,Hartanto,Imprialou,Rabanal)
      names(datasets) <- studies
      
  #for loop through studies
    for(study in studies){
      for(dtst in datasets[[study]]){
        #Load requiered datasets
        print(c('dtst',dtst))
          trait.matrix <- as.matrix(readRDS(file = paste0('files/',study,'/trait-matrix.',dtst,'.rds')))
          genetic.map <- as.matrix(readRDS(file = paste0('files/',study,'/genetic-map.',dtst,'.rds')))
          marker <- readRDS(paste0('files/',study,'/marker.',dtst,'.rds'))
          sample.list <- readRDS(paste0('files/',study,'/sample-list.',dtst,'.rds'))
          gene.info <- read.csv('files/trait.annotations.csv',row.names = 1)
          ril.dtst <- colnames(trait.matrix)
          trait <- trait.matrix
          map <- genetic.map

        #QTL mapping
          print(paste0('Starting processing of dataset ',dtst,' from ',study))
          n.cores <- detectCores() #- 1 
          print('running QTL.data.prep')
          qtl.data <- QTL.data.prep(trait.matrix = trait, 
                        strain.trait = colnames(trait), #Strain trait will still be the same one for all the datasets for now.
                        strain.map = map, 
                        strain.marker = marker)
          cluster <- makeCluster(n.cores, type = "PSOCK")
          registerDoParallel(cluster)
          output <- foreach(i = 1:nrow(trait), .combine = 'cbind') %dopar% {
            map.all.marker(trait = trait[i, ], markers = map)
          }
          stopCluster(cluster); print('cluster stopped')
          
          pval.out <- t(output[, which(colnames(output) == 'LOD')])
          eff.out <- t(output[, which(colnames(output) == 'Eff')])
          
          colnames(pval.out) <- rownames(marker); rownames(pval.out) <- rownames(trait)
          colnames(eff.out) <- rownames(marker); rownames(eff.out) <- rownames(trait)
          print(c('pval.pout',head(pval.out)))
          print(c('eff.out',head(eff.out)))
          print('colnames corrected')
          
          qtl.profile <- NULL; qtl.profile <- as.list(qtl.profile)
          qtl.profile[[1]] <- round(pval.out,digits=2)
          qtl.profile[[2]] <- round(eff.out,digits=3)
          qtl.profile[[3]] <- trait
          qtl.profile[[4]] <- map
          qtl.profile[[5]] <- marker
          names(qtl.profile) <- c("LOD","Effect","Trait","Map","Marker")
          print(c('qtl.profile'),str(qtl.profile))
          write.EleQTL(map1.output = qtl.profile, filename = paste0("qtl-tables/table_single-dtst-eqtl_", dtst))#The problem has been identified here.
          saveRDS(object = qtl.profile,file = paste0("qtl-profiles/profile_single-dtst-eqtl_", dtst,".rds")) 
          print(c('qtl.profile saved for',dtst))
              }
            }

# Nyan Cat ----------------------------------------------------------------
          browseURL('https://www.youtube.com/watch?v=QH2-TGUlwu4')
          str(qtl.profile)

          
# eQTL peak finder and table ####
  # threshold
  stud.list <- c('Keurentjes','West','West','Cubillos','Cubillos','Cubillos','Cubillos','Snoek','Lowry','Lowry','Hartanto','Hartanto','Hartanto','Hartanto','Imprialou','Rabanal')
  dtst.list <- c('keure','wesw1','wesw2','cub3b','cub3c','cub5b','cub5c','snoek','lowdr','lowwe','harar','harim','harpd','harrp','impri','raban')
  LOD.co.lst <- c(3.31,#int.qtl.keure
                  3.00,3.00,#int.qtl.wesw1,int.qtl.wesw2 
                  3.00,3.00,3.00,3.00,#int.qtl.cub3b,int.qtl.cub3c,int.qtl.cub5b,int.qtl.cub5c
                  3.09,#int.qtl.snoke
                  4.00,4.00,#int.qtl.lowdr,int.qtl.lowwe
                  4.2,4.1,4.2,4.3,#int.qtl.harar,int.qtl.harim,int.qtl.harpd,int.qtl.harrp
                  4.5, #Random 4.5 added for impar
                  4.5) #Random 4.5added to raban
  # I will have to re-run several datasets becasue they may have been run with the wrong threshold the first time around...
    #Revise this
  threshold.df <- data.frame('study'=stud.list,'dataset'=dtst.list,'LOD.co'=LOD.co.lst) #data frame with the lod cut off information.
  gene.info <- read.csv('files/trait.annotations.csv',row.names = 1) #This is needed for the loop below

  #For loop thorugh the studies and their datasets   
  for(study in studies){
    for(dtst in datasets[[study]]){
      print(paste0('Starting processing of dataset ',dtst,' from ',study))
      threshold <- threshold.df$LOD.co[threshold.df$dataset==dtst]
      #Loading its data
      print('1')
      trait.matrix <- as.matrix(readRDS(file = paste0('files/',study,'/trait-matrix.',dtst,'.rds')))
      genetic.map <- as.matrix(readRDS(file = paste0('files/',study,'/genetic-map.',dtst,'.rds')))
      marker <- readRDS(paste0('files/',study,'/marker.',dtst,'.rds'))
      sample.list <- readRDS(paste0('files/',study,'/sample-list.',dtst,'.rds'))
      gene.info <- read.csv('files/trait.annotations.csv',row.names = 1)
      ril.dtst <- colnames(trait.matrix)
      trait <- trait.matrix
      map <- genetic.map
      # the thresholds are based on multiple-testing correction using 100 permuted datasets
      print('2')
      qtl.profile <- readRDS(paste0('qtl-profiles/profile_single-dtst-eqtl_', dtst, '.rds'))
      qtl.peak <- mapping.to.list(map1.output = qtl.profile) %>% #The warning comes form this line
        peak.finder(threshold =  threshold) 
      print(c('threshold',threshold))
      qtl.peak <- na.omit(qtl.peak)
      saveRDS(object = qtl.peak,
              file = paste0("qtl-peaks/peak_single-dtst-eqtl_", dtst, ".rds"))
    
      ###eQTL table ###
      qtl.profile <- readRDS(paste0('qtl-profiles/profile_single-dtst-eqtl_', dtst, '.rds'))
      qtl.peak <- readRDS(paste0('qtl-peaks/peak_single-dtst-eqtl_', dtst, '.rds'))
      eqtl.table <- eQTL.table(peak.list.file = qtl.peak, trait.annotation = gene.info) %>%
        eQTL.table.addR2(QTL.prep.file = qtl.profile)
      eqtl.table$qtl_chromosome <- as.factor(eqtl.table$qtl_chromosome)
      eqtl.table$gene_chromosome <- as.factor(eqtl.table$gene_chromosome)
      tail(eqtl.table[,])
      
      dim(eqtl.table)
      ###heritability single stage removed###
      
      ###permutation ###section {removed}###
    
  # trans-bands identification ####
      print('3')
      window.nu <- 2e6
      maxsize <- 100e6
      chr.num <- 5
      print('mutate calculations')
      transband.id <- mutate(eqtl.table, interval = findInterval(qtl_bp, seq(1, maxsize, by = window.nu))) %>%
        group_by(qtl_chromosome, interval, qtl_type) %>%
        summarise(n.ct = length(unique(trait))) %>%
        data.frame() %>%
        group_by(qtl_type) %>%
        mutate(exp.ct = mean(as.numeric(unlist(n.ct)))) %>%
        data.frame() %>%
        mutate(transband_significance = ppois(n.ct, lambda = exp.ct, lower.tail = F)) %>%
        filter(transband_significance < 0.0001, qtl_type == "trans") 
      print('4')
      transband.id$transband_id <- with(transband.id, paste0("ch", qtl_chromosome, ":", 
                                                             (interval - 1) * 2, "-", 
                                                             interval * 2, "Mb"))
      transband.id$dtst <- dtst
      print('5')
      saveRDS(object = transband.id,
              file = paste0("trans-bands/trans.band_", dtst, ".rds"))
      
      #End the loop here for now. It will be a good spot to check out the results.
    }
  }
  
  test1 <- readRDS(paste0("qtl-peaks/peak_single-dtst-eqtl_", 'impri', ".rds"))
  test2 <- readRDS(paste0("trans-bands/trans.band_", 'impri', ".rds"))
  dim(test2)
  head(test1)
  head(eqtl.table)
  
head(eqtl.table)

# TransBand Results Analysis ----------------------------------------------
  for(study in studies){
    for(dtst in datasets[[study]]){
      #In the script, eqtl.table is not saved until later, so I have to run it here again
      ###eQTL table ###
        qtl.profile <- readRDS(paste0('qtl-profiles/profile_single-dtst-eqtl_', dtst, '.rds'))
        qtl.peak <- readRDS(paste0('qtl-peaks/peak_single-dtst-eqtl_', dtst, '.rds'))
        eqtl.table <- eQTL.table(peak.list.file = qtl.peak, trait.annotation = gene.info) %>%
          eQTL.table.addR2(QTL.prep.file = qtl.profile)
        eqtl.table$qtl_chromosome <- as.factor(eqtl.table$qtl_chromosome)
        eqtl.table$gene_chromosome <- as.factor(eqtl.table$gene_chromosome)
      #loading trans band file
      print(paste0('loading trans band file of ',dtst))
      tbs <- readRDS(paste0("trans-bands/trans.band_", dtst, ".rds"))

      #Interval table maker
        print('Filling interval table')
        window.nu <- 2e6; maxsize <- 100e6; chr.num <- 5 #This is the variables used earlier (make sure they match)
        chromosome <- tbs$qtl_chromosome
        interval <-   tbs$interval
        lower.limit <- seq(from = 0, to = 98e6,by =  2e6)
        upper.limit <- seq(from = 2e6, to = 100e6,by =  2e6)
        intervals <- data.frame(chromosome,interval,lower.limit[interval],upper.limit[interval])

      #trans.band annotation
      print('trans.band annotations on eqtl.table started')
       eqtl.table <- mutate(eqtl.table, trans_band = "none")
       
      for (tb in 1:nrow(tbs)) {
        trans.band <- tbs[tb,]
        print(paste0('transband ',tb,'being inserted'))   
        for (row in 1:nrow(eqtl.table)) {
            if(eqtl.table$qtl_type[row] == "trans" &
               eqtl.table$qtl_chromosome[row] == trans.band$qtl_chromosome &
               eqtl.table$qtl_bp[row] > intervals$lower.limit.interval.[tb] &
               eqtl.table$qtl_bp[row] <= intervals$upper.limit.interval.[tb])
                { eqtl.table$trans_band[row] <- trans.band$transband_id}
        }
                                       
      }
      print('Trans.band annotation ended,writing output')
      write.csv(x = eqtl.table,
                file = paste0("qtl-tables/table_single-dtst-eqtl_", dtst, ".csv"), 
                row.names = T)
      saveRDS(object = eqtl.table,
              file = paste0("qtl-tables/table_single-dtst-eqtl_", dtst, ".rds"))
      table(eqtl.table$qtl_type, eqtl.table$trans_band!="none")
    }
  }
  

    # GO for trans bands - single dtst ####
    ontology <- 'BP' #c('BP', 'CC', 'MF')
    
  trait.matrix <- readRDS('files/Hartanto/trait-matrix.harrp.rds')
  trait.matrix[1:10,]
    x <- org.At.tairCHR
    all.genes <- as.list(rownames(trait.matrix)) #loop through the triat matrices
    transband.go <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 9))
    colnames(transband.go) <- c('GO.ID', 'Term', 'Annotated', 'Significant', 'Expected', 'Fisher', 
                          'FDR', 'dtst', 'transband')
    transband.go.perdtst <- transband.go
    
    for (dtst in data.sets) {
      dtst <- 'keure'
      eqtl.table <- readRDS(paste0('qtl-tables/table_single-dtst-eqtl_', dtst, '.rds'))
      transband.id <- unique(eqtl.table$trans_band)
      transband.id <- transband.id[!transband.id %in% 'none']
      transband.go.perdtst <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 9))
      colnames(transband.go.perdtst) <- c('GO.ID', 'Term', 'Annotated', 'Significant', 'Expected', 'Fisher', 
                            'FDR', 'dtst', 'transband')
      
      for (i in 1:length(transband.id)) {
        gene.set <- eqtl.table[which(eqtl.table$trans_band == transband.id[i]), 'trait']
        gene.set <- factor(as.integer(all.genes %in% gene.set))
        names(gene.set) <- all.genes
        GOdata <- new("topGOdata",
                      description = "GOE for genes in regulated by trans bands", ontology = ontology,
                      allGenes = gene.set, 
                      annot = annFUN.org,mapping= "org.At.tair.db")
        resultFisher <- runTest(GOdata, algorithm = 'weight', statistic = "fisher")
        result.df <- GenTable(GOdata, Fisher = resultFisher,
                              orderBy = "Fisher", ranksOf = "Fisher", topNodes = length(resultFisher@score))
        result.df$dtst <- dtst
        result.df$transband <- transband.id[i]
        result.df$Fisher <- as.numeric(result.df$Fisher)
        result.df$FDR <- p.adjust(p = result.df$Fisher, method = 'fdr')
        result.df <- result.df[order(result.df$FDR), ]
        result.df <- dplyr::filter(result.df, Fisher <= 0.01)
        transband.go.perdtst <- rbind(transband.go.perdtst, result.df)
      }
      transband.go <- rbind(transband.go, transband.go.perdtst)
    }
    
    write.csv(transband.go, paste0('files/trans-bands-go-', ontology, '.csv'))
    
    help(topGO)
 