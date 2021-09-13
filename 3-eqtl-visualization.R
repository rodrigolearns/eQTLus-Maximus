##################################################
##### visualisation of eQTL analysis results #####
##################################################

#### prepare the script working environment ####
  remove(list = ls())
  gc()
  set.seed(1000)  
  
  # Set working directory ####
  work.dir <- "C:/Users/Rodrigo/Desktop/MiPRO/1_MASTER_DATA_MIPRO/2_Scripts/111 MyScripts/1A_Hartanto_Edited/"
  setwd(work.dir)

  # dependencies ####
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(doParallel)
  library(VennDiagram)
  library(UpSetR)
  library(forcats)
  library(grid)
  library(heritability)
  library(topGO)
  library(org.At.tair.db)
  # unused libraries ####
  # library(Mfuzz)
  # library(Biobase)
  # library(gplots)
  # library(RColorBrewer)
  # library(limma)
  # library(factoextra)
  # library(stats)
  # library(amap)
  # library(gridExtra)
  # library(RCy3)
  # library(corrplot)
  # library(reshape2)
  # library(Hmisc)
  # library(igraph)
  # library(threejs)
  # library(MASS)
  # 


  # load required function ####
  setwd('functions/')
  for(i in 1:length(dir())){
    source(dir()[i])
  } # read function from Mark
  setwd(work.dir)

# Overlap Function --------------------------------------------------------
  overlap <- function(table, dtst) {
    for (i in 1:length(dtst)) {
      print(dtst)
      print(c('nrow table',nrow(table)))
      if (nrow(table) != 0) {
        print(head((table[,dtst[i]])))
        table <- subset(table, table[,dtst[i]] == T)
      } else {
        break
      }
    }
    nrow(table)
    print(nrow(table))
  } # function to determine the overlapping qtl between 2 dtsts
  
  # load the data ####
    #seed.dtst changed to dtst.list
  dtst.list <- c('keure','wesw1','wesw2','cub3b','cub3c','cub5b','cub5c','snoek','lowdr','lowwe','harar','harim','harpd','harrp','impri','raban')
  

  # presentation
  presentation <- theme(axis.text.x = element_text(size=10, face="bold", color="black"),
                        axis.text.y = element_text(size=10, face="bold", color="black"),
                        axis.title.x = element_text(size=12, face="bold", color="black"),
                        axis.title.y = element_text(size=12, face="bold", color="black"),
                        strip.text.x = element_text(size=12, face="bold", color="black"),
                        strip.text.y = element_text(size=12, face="bold", color="black"),
                        strip.text = element_text(size =12, face="bold", color="black"),
                        #plot.title = element_text(size=15, face="bold"),
                        panel.background = element_rect(fill = "white",color="black"),
                        panel.grid.major = element_line(colour = "grey80"),
                        panel.grid.minor = element_blank())
  
# combine the eqtl table and plotting - single dtst ####
    
    all.eqtl.table <- data.frame(matrix(ncol = 20, nrow = 0))
    #colnames(all.eqtl.table) <- colnames(eqtl.table)

    for (dtst in dtst.list) {
      print(dtst)
      eqtl.table <- readRDS(paste0('qtl-tables/table_single-dtst-eqtl_', dtst, '.rds')) %>%
        mutate(dtst = dtst)
      all.eqtl.table <- rbind(all.eqtl.table, eqtl.table)
    }
    
    all.eqtl.table <- all.eqtl.table %>%
      mutate(qtl_type = ifelse(qtl_type == 'cis', 'local', 'distant'))
    
    saveRDS(all.eqtl.table, 'qtl-tables/table_single-dtst-eqtl_all.rds')
  
# comparing LOD, effect, and R2 of shared local and distant eQTL - single dtst ####
    eqtl.stat <- group_by(all.eqtl.table, dtst, qtl_type) %>%
                  summarise(median_sig = round(median(qtl_significance), 2),
                            #sd_sig = sd(qtl_significance),
                            median_eff = round(median(abs(qtl_effect)), 2),
                            #sd_eff = sd(qtl_effect),
                            median_r2 = round(median(qtl_R2_sm), 2)) %>%
                            #sd_r2 = sd(qtl_R2_sm),) %>%
                  mutate(pval_sig = NA, pval_eff = NA, pval_r2 = NA) %>%
                  as.data.frame()
    
    for (i in 1:nrow(eqtl.stat)) {
      dtst.row <- which(all.eqtl.table$dtst == eqtl.stat[i, 1])
      eqtl.stat[i, 6] <- wilcox.test(qtl_significance ~ qtl_type, data = all.eqtl.table[dtst.row, ])$p.value
      eqtl.stat[i, 7] <- wilcox.test(abs(qtl_effect) ~ qtl_type, data = all.eqtl.table[dtst.row, ])$p.value
      eqtl.stat[i, 8] <- wilcox.test(qtl_R2_sm ~ qtl_type, data = all.eqtl.table[dtst.row, ])$p.value
    }
    
# calculating shared local and distant eQTL ####
    
    qtl.table.summary <- data.frame(matrix(ncol = 18, nrow = 0))
    window.nu <- 4e6
    maxsize <- 100e6
    chr.num <- 5
  
    all.eqtl.table2 <- mutate(all.eqtl.table, interval = findInterval(qtl_bp, seq(1, maxsize, by = window.nu))) 
    all.eqtl.table2 <- subset(x = all.eqtl.table2, select = c(trait, qtl_chromosome, interval, qtl_type, 
                                                              dtst))
    
    local.table <- all.eqtl.table2[all.eqtl.table2$qtl_type == 'local', ] %>%
      count(trait, qtl_chromosome, interval, qtl_type, dtst) %>%
      spread(dtst, n) %>%
      dplyr::select(-qtl_type) %>%
      arrange(trait)
    duplicated.local <- duplicated(local.table$trait)
    
    for (i in 1:nrow(local.table)) {
      if (i == 1) {
        next
      }
      if (duplicated.local[i] && 
          local.table$qtl_chromosome[i] == local.table$qtl_chromosome[i-1] &&
          abs(local.table$interval[i] - local.table$interval[i-1]) == 1 &&
          local.table$trait[i] == local.table$trait[i-1]) {
        print(i)
        table.tmp <- rbind(local.table[i, ], local.table[i-1, ])
        table.tmp <- aggregate(table.tmp[4:ncol(local.table)], by = list(trait = table.tmp$trait), sum, na.rm = TRUE)
        table.tmp$qtl_chromosome <- local.table$qtl_chromosome[i]
        table.tmp$interval <- local.table$interval[i]
        table.tmp <- table.tmp[, c(colnames(local.table))]
        local.table[(i-1), 1] <- 'removed'
        local.table[i, ] <- table.tmp
      }
    }
    
    local.table <- local.table[which(local.table$trait != 'removed'), ]
    local.table[is.na(local.table)] <- 0
    write.csv(local.table, 'files/shared-local.csv')
    local.table <- local.table[, 4:ncol(local.table)]
    local.table <- apply(local.table, 2, as.integer)
    
    distant.table <- all.eqtl.table2[all.eqtl.table2$qtl_type == 'distant', ] %>%
      count(trait, qtl_chromosome, interval, qtl_type, dtst) %>%
      spread(dtst, n) %>%
      dplyr::select(-qtl_type) %>%
      arrange(trait)
    duplicated.distant <- duplicated(distant.table$trait)
    
    for (i in 1:nrow(distant.table)) {
      if (i == 1) {
        next
      }
      if (duplicated.local[i] && 
          distant.table$qtl_chromosome[i] == distant.table$qtl_chromosome[i-1] &&
          abs(distant.table$interval[i] - distant.table$interval[i-1]) == 1 &&
          distant.table$trait[i] == distant.table$trait[i-1]) {
        print(i)
        table.tmp <- rbind(distant.table[i, ], distant.table[i-1, ])
        table.tmp <- aggregate(table.tmp[4:ncol(local.table)], by = list(trait = table.tmp$trait), sum, na.rm = TRUE)
        table.tmp$qtl_chromosome <- distant.table$qtl_chromosome[i]
        table.tmp$interval <- distant.table$interval[i]
        table.tmp <- table.tmp[, c(colnames(distant.table))]
        distant.table[(i-1), 1] <- 'removed'
        distant.table[i, ] <- table.tmp
      }
    }
    
    distant.table <- distant.table[which(distant.table$trait != 'removed'), ]
    distant.table[is.na(distant.table)] <- 0
    write.csv(distant.table, 'files/shared-distant.csv')
    distant.table <- distant.table[, 4:ncol(local.table)]
    distant.table <- apply(distant.table, 2, as.integer)
    
    #subset(qtl.table.cis, qtl.table.cis[development.dtst[1]] == T)
    
    # venn diagram ####
        #Problem. All emty...

     draw.quad.venn(area1 = overlap(local.table, 'pd'), 
                   area3 = overlap(local.table, 'ar'), 
                   area4 = overlap(local.table, 'im'), 
                   area2 = overlap(local.table, 'rp'), 
                   n13 = overlap(local.table, c('pd', 'ar')), 
                   n14 = overlap(local.table, c('pd', 'im')), 
                   n12 = overlap(local.table, c('pd', 'rp')), 
                   n34 = overlap(local.table, c('ar', 'im')), 
                   n23 = overlap(local.table, c('ar', 'rp')), 
                   n24 = overlap(local.table, c('im', 'rp')), 
                   n134 = overlap(local.table, c('pd', 'ar', 'im')), 
                   n123 = overlap(local.table, c('pd', 'ar', 'rp')), 
                   n124 = overlap(local.table, c('pd', 'im', 'rp')), 
                   n234 = overlap(local.table, c('ar', 'im', 'rp')), 
                   n1234 = overlap(local.table, c('pd', 'ar', 'im', 'rp')), 
                   category = c('primary\ndormant', 'radicle\nprotrusion', 'after-\nripenned', '6 hours after\nimbibition'),
                   lty = 'blank', 
                   fill = c('#4daf4a', '#e41a1c', '#377eb8', '#984ea3'),
                   alpha = 0.5)
    
    draw.quad.venn(area1 = overlap(distant.table, 'pd'), 
                   area3 = overlap(distant.table, 'ar'), 
                   area4 = overlap(distant.table, 'im'), 
                   area2 = overlap(distant.table, 'rp'), 
                   n13 = overlap(distant.table, c('pd', 'ar')), 
                   n14 = overlap(distant.table, c('pd', 'im')), 
                   n12 = overlap(distant.table, c('pd', 'rp')), 
                   n34 = overlap(distant.table, c('ar', 'im')), 
                   n23 = overlap(distant.table, c('ar', 'rp')), 
                   n24 = overlap(distant.table, c('im', 'rp')), 
                   n134 = overlap(distant.table, c('pd', 'ar', 'im')), 
                   n123 = overlap(distant.table, c('pd', 'ar', 'rp')), 
                   n124 = overlap(distant.table, c('pd', 'im', 'rp')), 
                   n234 = overlap(distant.table, c('ar', 'im', 'rp')), 
                   n1234 = overlap(distant.table, c('pd', 'ar', 'im', 'rp')), 
                   category = c('primary\ndormant', 'radicle\nprotrusion', 'after-\nripenned', '6 hours after\nimbibition'),
                   lty = 'blank', 
                   fill = c('#4daf4a', '#e41a1c', '#377eb8', '#984ea3'),
                   alpha = 0.5, )
    
    pdf(file = 'figures/distant-local-venn-diagram.pdf', width = 10, height = 5)
    grid.arrange(local.plot, distant.plot, ncol = 2, heights = c(5, 5))
    dev.off() #Its saving wrong but its ok idgaf about this for now.
    head(local.table)
    # upset graph ####
    dim(local.table)
    upset.local <- upset(data = as.data.frame(local.table), 
                         mainbar.y.max = 750, 
                        sets =  c('cub3b' ,'cub3c' ,'cub5b', 'cub5c' ,'harar' ,'harim' ,'harpd' ,'harrp' ,'impri' ,'keure', 'lowdr' ,'lowwe' ,'raban' ,'snoek','wesw1','wesw2'),
                        empty.intersections = 'on',
                        mainbar.y.label = 'shared eQTL',
                        sets.x.label = 'eQTL per dtst', 
                        text.scale = 0.8, 
                        point.size = 1, 
                        line.size = 1,
                        keep.order = T,
                        set_size.scale_max = 1400)
    
    tiff(file = paste0("figures/upset-local.tiff"), 
         width = 1000, 
         height = 850, 
         units = 'px',
         res = 300,
         compression = 'lzw')
    upset.local
    dev.off()
    
    
    upset.distant <- upset(data = as.data.frame(distant.table), 
          mainbar.y.max = 750,
          sets = c('cub3b' ,'cub3c' ,'cub5b', 'cub5c' ,'harar' ,'harim' ,'harpd' ,'harrp' ,'impri' ,'keure', 'lowdr' ,'lowwe' ,'raban' ,'snoek','wesw1','wesw2'),
          empty.intersections = 'on',
          mainbar.y.label = 'shared eQTL',
          sets.x.label = 'eQTL per dtst', 
          text.scale = 0.8, 
          point.size = 1, 
          line.size = 1,
          keep.order = T, 
          set_size.scale_max = 1400)

    tiff(file = paste0("figures/upset-distant.tiff"), 
         width = 1000, 
         height = 850, 
         units = 'px',
         res = 300,
         compression = 'lzw')
    upset.distant
    dev.off()
    
# cis-trans plot - single dtst #####
    
    all.eqtl.table3 <- all.eqtl.table %>%
      mutate(dtst = replace(dtst, qtl_type == 'local', 'local\neQTL'))
    #all.eqtl.table3$dtst <- ordered(all.eqtl.table$dtst, dtst.list)
    
    all.eqtl.table3$qtl_type <- factor(all.eqtl.table3$qtl_type)
    #all.eqtl.table3 <- all.eqtl.table3 %>%
    #  mutate(allele = ifelse(sign(qtl_effect) <= 0, 'sha', #this is the funky part where it is determined where the allele came from.
    #                         ifelse(sign(qtl_effect) >=0, 'bay', 'null')))
    #head(all.eqtl.table3)
    unique(all.eqtl.table3$dtst)
    eqtl.plot <- ggplot(all.eqtl.table3, aes(x = qtl_bp, y = gene_bp)) + 
      geom_segment(aes(x = qtl_bp_left, y = gene_bp, xend = qtl_bp_right, yend = gene_bp),
                   alpha = 0.75, colour = "grey") +
      geom_point(aes(colour = ordered(dtst, levels = c("local\neQTL", 'cub3b' ,'cub3c' ,'cub5b', 'cub5c' ,'harar' ,'harim' ,'harpd' ,'harrp' ,'impri' ,'keure', 'lowdr' ,'lowwe' ,'raban' ,'snoek','wesw1','wesw2')) 
                     #The part above is all fucked up because it is study specific...
                     
                     #When looking at all.eqtl.table3 I realize that it indicates all eqtls allele as bay or sha... look where that was datermined and change it if it is used in somethinf important
                     
                     #, shape = allele,
                     #size = abs(qtl_effect),
                     #alpha = log10(qtl_significance))
                 )) +
      #scale_alpha(range = c(0.3, 1)) +
      facet_grid(fct_rev(gene_chromosome) ~ qtl_chromosome, space = "free", scales = "free") +
      presentation +
      theme(legend.title= element_text(size = 12, face = "bold"),
            legend.text = element_text(size = 12)) +
      scale_colour_manual('dtst', values = c( "#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", 
                                              "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
                                              "#D14285")) +
      labs(x = "eQTL peak position (Mb)", y = "Gene position (Mb)") +
      scale_x_continuous(breaks=c(5, 10, 15, 20, 25, 30, 35, 40)*10^6,labels=c(5, 10, 15, 20, 25, 30, 35, 40)) +
      scale_y_continuous(breaks=c(5, 10, 15, 20, 25, 30, 35, 40)*10^6,labels=c(5, 10, 15, 20, 25, 30, 35, 40))
    

    tiff(file = paste0("figures/eqtl-plot.tiff"), 
         width = 2250, 
         height = 1600, 
         units = 'px',
         res = 300,
         compression = 'lzw')
    
    eqtl.plot
    dev.off()
    browseURL('https://www.youtube.com/watch?v=QH2-TGUlwu4')
# histogram - single dtst ####
     
    hist.threshold <- data.frame(dtst =  dtst.list, exp = c(7.45098, 5.464286, 7, 6.470588), threshold = rep(NA, 4))
    for (i in 1:nrow(hist.threshold)) {
      hist.threshold$threshold[i] <- qpois(p = 0.0001, lambda = hist.threshold$exp[i], lower.tail = F)
    } # create horizontal line as a threshold
    
    
    eqtl.hist <- ggplot(all.eqtl.table3 %>% filter(qtl_type != 'local'), 
                        aes(x=qtl_bp, fill=ordered(dtst, levels = c("pd", "ar", "im", "rp")))) + # 750 x 400
      geom_histogram(binwidth = 2000000, right = T, origin = 0, alpha = 1) +
      facet_grid(factor(dtst, levels =  c("pd", "ar", "im", "rp")) ~ qtl_chromosome, space = "free",scales="free") +
      presentation +
      scale_fill_manual('dtst', values = c('#ccbb44', '#228833', '#4477aa', '#cc3311')) +
      theme(legend.position = "right", 
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12, face = 'bold')) +
      labs(x="eQTL peak position (Mb)",y="eQTL counts") +
      ylim(c(0, 100)) +
      geom_hline(data = hist.threshold, aes(yintercept = threshold), 
                 linetype = 'dashed', 
                 color = 'black', 
                 size = .5) +
      scale_x_continuous(breaks=c(5, 10, 15, 20, 25, 30, 35, 40)*10^6,labels=c(5, 10, 15, 20, 25, 30, 35, 40))
    
    tiff(file = paste0("figures/eqtl-hist.tiff"), 
         width = 2250, 
         height = 1200, 
         units = 'px',
         res = 300,
         compression = 'lzw')
    eqtl.hist
    dev.off()
    
# compare the R2, effect, andsignificance of local and distant eQTL ####
    r2.plot <- ggplot(all.eqtl.table, aes(x = qtl_R2_sm,
                                             fill = factor(qtl_type, levels = c('local', 'distant')),
                                             color = factor(qtl_type, levels = c('local', 'distant')))) +
      #geom_histogram(aes(y = ..density..), binwidth = 0.001, position = 'identity', alpha = 0.3) +
      geom_density(alpha = 0.5, size = 0) +
      facet_wrap( ~ factor(dtst, levels = c('pd', 'ar', 'im', 'rp')), ncol = 2) +
      geom_vline(data = group_by(all.eqtl.table, dtst, qtl_type) %>%
                   summarise(median = median(qtl_R2_sm)), aes(xintercept = median, 
                                                              color = factor(qtl_type, levels = c('local', 'distant'))),
                 linetype = 'dashed', size = .5) +
      presentation +
      scale_color_brewer(palette = 'Set1') +
      scale_fill_brewer(palette = 'Set1') +
      xlab('explained phenotypic variance (R2)') +
      ylab('density') +
      labs(fill = 'eQTL type', color = 'eQTL type') +
      #ylim(0, 1) + 
      xlim(0, 1) +
      theme(strip.text.x = element_text(size = 12))
    
    tiff(file = paste0("figures/r2-plot.tiff"), 
         width = 2250, 
         height = 1200, 
         units = 'px',
         res = 300,
         compression = 'lzw')
    r2.plot
    dev.off()
      
    effect.plot <- ggplot(all.eqtl.table, aes(x = log(abs(qtl_effect)),
                               fill = factor(qtl_type, levels = c('local', 'distant')),
                               color = factor(qtl_type, levels = c('local', 'distant')))) +
      #geom_histogram(aes(y = ..density..), binwidth = 0.001, position = 'identity', alpha = 0.3) +
      geom_density(alpha = 0.5, size = 0) +
      facet_wrap( ~ factor(dtst, levels = c('pd', 'ar', 'im', 'rp')), ncol = 2) +
      geom_vline(data = group_by(all.eqtl.table, dtst, qtl_type) %>%
                   summarise(median = median(log(abs(qtl_effect)))), aes(xintercept = median, 
                                                              color = factor(qtl_type, levels = c('local', 'distant'))),
                 linetype = 'dashed', size = .5) +
      presentation +
      scale_color_brewer(palette = 'Set1') +
      scale_fill_brewer(palette = 'Set1') +
      xlab('absolute eQTL effect') +
      ylab('density') +
      labs(fill = 'eQTL type', color = 'eQTL type') +
      #ylim(0, 1) + 
      #xlim(0, 1.2) +
      theme(strip.text.x = element_text(size = 12))
    
    tiff(file = paste0("figures/eff-plot.tiff"), 
         width = 2250, 
         height = 1200, 
         units = 'px',
         res = 300,
         compression = 'lzw')
    effect.plot
    dev.off()
    

    sig.plot <- ggplot(all.eqtl.table, aes(x = qtl_significance,
                               fill = factor(qtl_type, levels = c('local', 'distant')),
                               color = factor(qtl_type, levels = c('local', 'distant')))) +
      #geom_histogram(aes(y = ..density..), binwidth = 0.001, position = 'identity', alpha = 0.3) +
      geom_density(alpha = 0.5, size = 0) +
      facet_wrap( ~ factor(dtst, levels = c('pd', 'ar', 'im', 'rp')), ncol = 2) +
      geom_vline(data = group_by(all.eqtl.table, dtst, qtl_type) %>%
                   summarise(median = median(qtl_significance)), aes(xintercept = median, 
                                                               color = factor(qtl_type, levels = c('local', 'distant'))),
                 linetype = 'dashed', size = .5) +
      presentation +
      scale_color_brewer(palette = 'Set1') +
      scale_fill_brewer(palette = 'Set1') +
      xlab('-log10(p)') +
      ylab('density') +
      labs(fill = 'eQTL type', color = 'eQTL type') +
      #ylim(0, 1) + 
      xlim(0, 36) +
      theme(strip.text.x = element_text(size = 12))
    
    tiff(file = paste0("figures/sig-plot.tiff"), 
         width = 2250, 
         height = 1200, 
         units = 'px',
         res = 300,
         compression = 'lzw')
    sig.plot
    dev.off()
    
    # ggplot(all.eqtl.table, aes(x = h2_REML, y = qtl_R2_sm,
    #                                             fill = factor(qtl_type, levels = c('local', 'distant')),
    #                                             color = factor(qtl_type, levels = c('local', 'distant')))) +
    #   geom_point(size = 1, alpha = 0.5) +
    #   facet_wrap( ~ factor(dtst, levels = c('pd', 'ar', 'im', 'rp')), ncol = 2) +
    #   geom_vline(data = group_by(all.eqtl.table, dtst, qtl_type) %>%
    #                summarise(med = median(h2_REML, na.rm = T)), aes(xintercept = med, 
    #                                                                 color = factor(qtl_type, levels = c('local', 'distant'))),
    #              linetype = 'dashed', size = .5) +
    #   geom_hline(data = group_by(all.eqtl.table, dtst, qtl_type) %>%
    #                summarise(med = median(qtl_R2_sm, na.rm = T)), aes(yintercept = med, 
    #                                                                   color = factor(qtl_type, levels = c('local', 'distant'))),
    #              linetype = 'dashed', size = .5) +
    #   geom_smooth(model = lm) +
    #   geom_abline(aes(slope = 1, intercept = 0), size = .75, linetype = 'dashed') +
    #   presentation +
    #   scale_color_brewer(palette = 'Set1') +
    #   scale_fill_brewer(palette = 'Set1') +
    #   xlab('narrow-sense heritabllity (h2)') +
    #   ylab('explained phenotypic variance (R2)') +
    #   labs(fill = 'eQTL type', color = 'eQTL type') +
    #   ylim(0, 1) + 
    #   xlim(0, 1) + 
    #   theme(strip.text.x = element_text(size = 8))
 