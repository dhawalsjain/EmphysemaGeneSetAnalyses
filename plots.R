rm(list = ls())
setwd("~/external/AdelPeter")

# installs
if(F){
  pckgs <- c('ggplot2','plotly','reshape2','dplyr','knitr',
             'data.table','DT','gprofiler2','RColorBrewer',
             'ggpubr','Seurat')
  pckgs <- pckgs[!pckgs%in%installed.packages()]
  install.packages(pckgs,dependencies = T)
  rm(pckgs)
  
  pckgs <- c('gage','pathview','gageData','org.Hs.eg.db')
  pckgs <- pckgs[!pckgs%in%installed.packages()]
  BiocManager::install(pckgs)
  rm(pckgs)
}



## enrichments using gprofiler2
if(F){
  library(gprofiler2)
  d <- read.csv("LTRC_DE_dtabsGene_10_26_2022.csv")  
  
  dup <- d[d$adj.P.Val<0.1 & d$logFC>0, ]
  res <- gprofiler2::gost(query = dup$hgnc_symbol,organism = 'hsapiens',
                          ordered_query = F,exclude_iea = T,
                          domain_scope = 'annotated',significant = F,
                          user_threshold = 0.05,correction_method = 'g_SCS')
  table(res$result[res$result$significant==T,]$source)
  hterms = res$result[res$result$significant==T & res$result$source%in%c('KEGG','WP','REAC'),]$term_id
  p = gostplot(gostres = res,capped = F,interactive = F)
  p = publish_gostplot(p,highlight_terms = hterms)
  p
  
  ddn <- d[d$adj.P.Val<0.1 & d$logFC<0, ]
  resa <- gprofiler2::gost(query = ddn$hgnc_symbol,organism = 'hsapiens',
                          ordered_query = F,exclude_iea = T,
                          domain_scope = 'annotated',significant = F,
                          user_threshold = 0.05,correction_method = 'g_SCS')
  table(resa$result[resa$result$significant==T,]$source)
  htermsa = resa$result[resa$result$significant==T & resa$result$source%in%c('KEGG','WP','REAC'),]$term_id
  q = gostplot(gostres = resa,capped = F,interactive = F)
  q = publish_gostplot(q,highlight_terms = htermsa)
  q
  
  
}


## KEGG pathways, test plots
if(F){
  rm(list=ls())
  library('gage')
  library('pathview')
  library('gageData')
  library(org.Hs.eg.db)
  library(dplyr)
  data("kegg.sets.hs")
  data("sigmet.idx.hs")
  
  kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
  #names(kegg.sets.hs)
  kegg_human = kegg.gsets(species = 'hsa',id.type = 'kegg')
  kegg.sets.hs = kegg_human$kg.sets[unique(c(kegg_human$sigmet.idx,kegg_human$sig.idx,kegg_human$met.idx,kegg_human$dise.idx))]
  #names(kegg.sets.hs)
  
  deg <- read.csv("LTRC_DE_dtabsGene_10_26_2022.csv")  
  
  
  #--- upregulated genes
  res <- deg[deg$logFC>0 & deg$adj.P.Val < 0.1,]
  res$entrez = mapIds(org.Hs.eg.db,keys = res$hgnc_symbol,column = 'ENTREZID',keytype = 'SYMBOL',multiVals = 'first')
  res <- res[!is.na(res$entrez),]
  foldchages = res$logFC
  names(foldchages) = res$entrez
  #head(foldchages)
  keggres = gage(foldchages,gsets = kegg.sets.hs,same.dir = T)
  #lapply(keggres,head)
  keggrespathways = data.frame(id=rownames(keggres$greater),keggres$greater) %>%
    tbl_df() %>%
    filter(row_number()<=10) %>%
    .$id %>%
    as.character()
  #keggrespathways
  keggresids = substr(keggrespathways,start = 1,stop = 8)
  
  #--- downregulated genes
  res <- deg[deg$logFC>0 & deg$adj.P.Val < 0.1,]
  res$entrez = mapIds(org.Hs.eg.db,keys = res$hgnc_symbol,column = 'ENTREZID',keytype = 'SYMBOL',multiVals = 'first')
  res <- res[!is.na(res$entrez),]
  foldchages = res$logFC
  names(foldchages) = res$entrez
  #head(foldchages)
  keggres = gage(foldchages,gsets = kegg.sets.hs,same.dir = T)
  #lapply(keggres,head)
  keggrespathways = data.frame(id=rownames(keggres$greater),keggres$greater) %>%
    tbl_df() %>%
    filter(row_number()<=10) %>%
    .$id %>%
    as.character()
  #keggrespathways
  keggresids = substr(keggrespathways,start = 1,stop = 8)
  
  

  ##--- selected pathways, log2FC
  res <- deg
  res$entrez = mapIds(org.Hs.eg.db,keys = res$hgnc_symbol,column = 'ENTREZID',keytype = 'SYMBOL',multiVals = 'first')
  res <- res[!is.na(res$entrez),]
  foldchages = res$logFC
  names(foldchages) = res$entrez
  range(foldchages)
  tmp=pathview(gene.data = foldchages,pathway.id = 'hsa04068',species = 'hsa',
               limit = list(0.005,0.005))
  tmp=pathview(gene.data = foldchages,pathway.id = 'hsa04350',species = 'hsa',
               limit = list(0.005,0.005))
  
  ##--- selected pathways, FDR
  foldchages = -log10(res$adj.P.Val)
  names(foldchages) = res$entrez
  tmp=pathview(gene.data = foldchages,pathway.id = 'hsa04068',species = 'hsa',
               limit = list(5,5),mid = list('#ffffcc','#ffffcc'),
               low=list('white','white'),high = list('#800026','#800026'))
  tmp=pathview(gene.data = foldchages,pathway.id = 'hsa04350',species = 'hsa',
               limit = list(5,5),mid = list('#ffffcc','#ffffcc'),
               low=list('white','white'),high = list('#800026','#800026'))

  
}

## specificity of gene expression
if(F){
  library(data.table)
  library(Seurat)
  so <- readRDS("/home/rstudio/data/SCS/Kaminski_IPF2020/Kaminski_IPF_2020.rds")
  so@active.ident <- factor(so$cell_type)
  names(so@active.ident) <- rownames(so@meta.data)
  so$Disease <- so$Disease_Identity
  
  DimPlot(so, reduction = "umap", label = TRUE,repel=T, pt.size = 0.5,raster = T) + NoLegend()
  
  res <- deg
  res$entrez = mapIds(org.Hs.eg.db,keys = res$hgnc_symbol,column = 'ENTREZID',keytype = 'SYMBOL',multiVals = 'first')
  res <- res[!is.na(res$entrez),]
  res$entrez <- as.character(res$entrez)
  #resa <- res[res$entrez %in%   kegg.sets.hs[[grep('hsa04068',names(kegg.sets.hs))]], ]
  resa <- res[res$adj.P.Val<0.05 & res$logFC>0, ]
  resb <- res[res$adj.P.Val<0.05 & res$logFC < 0, ]
  
  qqa <- local({
    gc()
    e <- FetchData(so,vars = unique(resa$hgnc_symbol))
    cntr_thresh <- 0
    #plot(density(unlist(e),na.rm=T))
    dat <- so@active.ident
    dat <- data.frame(cell=names(dat),cluster=dat)
    tsnepl <- so@reductions$umap@cell.embeddings %>% as.data.frame
    tsnepl$cell <- rownames(tsnepl)
    tsnepl <- merge(tsnepl,dat,by="cell")
    e <- e %>% mutate(cell = rownames(e))
    dat <- merge(dat,e,by="cell",all.x=T)
    dat$cell <- NULL
    dat[is.na(dat)] <- 0
    dat <- as.data.table(dat)
    ## select expressed proteins in single cells
    dat <- dat[, lapply(.SD, function(x) sum(x>cntr_thresh) ), by=list(cluster)] %>%as.data.frame 
    rownames(dat) <- dat$cluster
    dat$cluster <- NULL
    dat <- dat %>% t %>% as.data.frame
    tsnepl <- merge(tsnepl,e,by="cell")
    
    y <- unname(table(tsnepl$cluster)[names(dat)]) 
    pvals <- t(apply(dat,1,function(x){
      out <- c()
      for(i in 1:length(x)){
        xx <- matrix(c( x[i], y[i]-x[i],  sum(x),sum(y)-sum(x)), nrow = 2,byrow = F)
        out <- c(out, fisher.test(xx,alternative = "greater")$p.val)
      }
      out
    }))
    pvals <- as.data.frame(pvals)
    padj <- matrix(p.adjust(unlist(pvals),method = "fdr"),nrow=nrow(pvals),byrow = F) %>% as.data.frame
    names(padj) <- names(pvals)
    rownames(padj) <- rownames(pvals)
    names(padj) <- names(pvals) <- paste0(names(dat))
    qq=list(pvals=pvals,padj=padj)
  })
  qqb <- local({
    gc()
    e <- FetchData(so,vars = unique(resb$hgnc_symbol))
    cntr_thresh <- 0
    #plot(density(unlist(e),na.rm=T))
    dat <- so@active.ident
    dat <- data.frame(cell=names(dat),cluster=dat)
    tsnepl <- so@reductions$umap@cell.embeddings %>% as.data.frame
    tsnepl$cell <- rownames(tsnepl)
    tsnepl <- merge(tsnepl,dat,by="cell")
    e <- e %>% mutate(cell = rownames(e))
    dat <- merge(dat,e,by="cell",all.x=T)
    dat$cell <- NULL
    dat[is.na(dat)] <- 0
    dat <- as.data.table(dat)
    ## select expressed proteins in single cells
    dat <- dat[, lapply(.SD, function(x) sum(x>cntr_thresh) ), by=list(cluster)] %>%as.data.frame 
    rownames(dat) <- dat$cluster
    dat$cluster <- NULL
    dat <- dat %>% t %>% as.data.frame
    tsnepl <- merge(tsnepl,e,by="cell")
    
    y <- unname(table(tsnepl$cluster)[names(dat)]) 
    pvals <- t(apply(dat,1,function(x){
      out <- c()
      for(i in 1:length(x)){
        xx <- matrix(c( x[i], y[i]-x[i],  sum(x),sum(y)-sum(x)), nrow = 2,byrow = F)
        out <- c(out, fisher.test(xx,alternative = "greater")$p.val)
      }
      out
    }))
    pvals <- as.data.frame(pvals)
    padj <- matrix(p.adjust(unlist(pvals),method = "fdr"),nrow=nrow(pvals),byrow = F) %>% as.data.frame
    names(padj) <- names(pvals)
    rownames(padj) <- rownames(pvals)
    names(padj) <- names(pvals) <- paste0(names(dat))
    qq=list(pvals=pvals,padj=padj)
  })
  
  
  m <- -log10(qqa[[2]])
  m <- as.matrix(m)
  m[is.infinite(m)] <- 100 
  m <- ifelse(m>100,100,m)
  ComplexHeatmap::Heatmap(t(m),
                          #column_order = names(n),
                          name = 'enrichment\n-log10(FDR)',
                          row_names_gp = grid::gpar(fontsize = 5),
                          column_names_gp = grid::gpar(fontsize = 5),
                          column_title_gp =  grid::gpar(fontsize = 6),
                          col=circlize::colorRamp2(range(m,na.rm=T),
                                                   c("#ffffcc", "#800026")))
  
  
  
}

## Projecting pathway score on Kaminski data, test cases
if(F){
  kegg_human = kegg.gsets(species = 'hsa',id.type = 'kegg')
  kegg.sets.hs = kegg_human$kg.sets[unique(c(kegg_human$sigmet.idx,kegg_human$sig.idx,kegg_human$met.idx,kegg_human$dise.idx))]
  
  res <- deg
  res$entrez = mapIds(org.Hs.eg.db,keys = res$hgnc_symbol,column = 'ENTREZID',keytype = 'SYMBOL',multiVals = 'first')
  res <- res[!is.na(res$entrez),]
  res$entrez <- as.character(res$entrez)
  
  resa <- res[res$entrez %in%   kegg.sets.hs[[grep('hsa04068',names(kegg.sets.hs))]], ]
  so <- AddModuleScore(so,features = list(resa$hgnc_symbol),
                         name="FoxO signalling")
  
  resa <- res[res$entrez %in%   kegg.sets.hs[[grep('hsa04350',names(kegg.sets.hs))]], ]
  so <- AddModuleScore(so,features = list(resa$hgnc_symbol),
                       name="TGF Beta signalling")
  
  FeaturePlot(so,features = c("FoxO.signalling1"), 
              label = T, repel = T,reduction = 'umap')+
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  
  VlnPlot(so,features = c("FoxO.signalling1"),split.by = 'Disease')  
  
  library(ggplot2)
  
  pl <- so@meta.data
  pl[!pl$cell_type %in% c('Aberrant_Basaloid','Ionocyte'),] %>%
    ggplot(aes(x=cell_type,y=FoxO.signalling1,fill=Disease))+
    geom_boxplot(outlier.shape = NA)+
    theme_bw()+
    ggtitle('FoxO Signalling')+ylab("")+xlab("Avg. normalized expression")+
    geom_hline(yintercept=median(pl$FoxO.signalling1),col='red')+
    #scale_color_manual(values=c('Control'='#91cf60','COPD'='#fc8d59'))+
    scale_fill_manual(values=c('Control'='#91cf60','COPD'='#fc8d59'))+
    stat_compare_means(label = "p.signif", method = "t.test",
                        hide.ns = T,label.y=seq(0.2,0.35,length.out=37)) + 
    theme(axis.text=element_text(color='black'),
          axis.text.x=element_text(angle=90,hjust=1))
  
  
  
  
  kegg.sets.hs[[grep('hsa04350',names(kegg.sets.hs))]]
  
}

## composite pathway expression changes
if(F){
  library(data.table)
  
  pl <- so@meta.data
  pl$foxo.score  <- ifelse(pl$FoxO.signalling1<quantile(pl$FoxO.signalling1,probs=c(0.66)),0,1)
  pl$tgf.score <- ifelse(pl$TGF.Beta.signalling1<quantile(pl$TGF.Beta.signalling1,probs=c(0.66)),0,1)
  pl <- pl[,c('cell_type','foxo.score','tgf.score')]
  pl <- data.table(pl)
  pl <- pl[,lapply(.SD,sum),by=list(cell_type)]
  pl <- as.data.frame(pl)
  pl$foxo.perc <- round(pl$foxo.score*100/sum(pl$foxo.score),2)
  pl$tgf.perc <- round(pl$tgf.score*100/sum(pl$tgf.score),2)
  
  pl <- pl[order(pl$tgf.perc),]
  pl$cell_type <- factor(pl$cell_type,levels=pl$cell_type)
  
  pl %>% 
    ggplot(aes(x=cell_type,y=tgf.perc,fill=tgf.perc))+
    geom_bar(stat='identity')+
    theme_bw()+
    xlab("")+ylab("specificity, in%")+ggtitle('TGF-Beta Signalling score')+
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Blues")))+
    theme(axis.text=element_text(color='black'),
          axis.text.x=element_text(angle=90,hjust=1))
  
  
  
  
  res <- deg
  res$entrez = mapIds(org.Hs.eg.db,keys = res$hgnc_symbol,column = 'ENTREZID',keytype = 'SYMBOL',multiVals = 'first')
  res <- res[!is.na(res$entrez),]
  res$entrez <- as.character(res$entrez)
  resa <- res[res$entrez %in%   kegg.sets.hs[[grep('hsa04068',names(kegg.sets.hs))]], ]
  resa
  
}


##-------------------------------------------
## additional pathways to look up
##-------------------------------------------
if(F){
  data("kegg.sets.hs")
  data("sigmet.idx.hs")
  kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
  #names(kegg.sets.hs)
  kegg_human = kegg.gsets(species = 'hsa',id.type = 'kegg')
  kegg.sets.hs = kegg_human$kg.sets[unique(c(kegg_human$sigmet.idx,kegg_human$sig.idx,kegg_human$met.idx,kegg_human$dise.idx))]
  #names(kegg.sets.hs)
  
  #mypaths <- read.csv("ltrc_egsea_camera_kegg.csv",header=T,sep=",")
  mypaths <- read.csv("LTRC_pathways_10_26_2022.csv",header=T,sep=",")
  mypaths$pathway <- mypaths$X
  mypaths$kegg <- NA
  for(i in 1:nrow(mypaths)){
    a <- names(kegg.sets.hs)[grep(mypaths$pathway[i],names(kegg.sets.hs))]
    if(length(a)>0){
      mypaths$kegg[i] <- paste0(a,collapse=";")
    }
    rm(a)
  }
  rm(i)
  mypaths$kegg[1] <- 'hsa03010 Ribosome'
  names(kegg.sets.hs)[grep('hsa04630',names(kegg.sets.hs))]
  mypaths$kegg[20] <-  names(kegg.sets.hs)[grep('hsa04630',names(kegg.sets.hs))]
  dput(mypaths[is.na(mypaths$kegg),]$pathway)
  #- pathways not found
  #c("Parkinson's disease", "Huntington's disease", "Alzheimer's disease", 
  #"Dorso-ventral axis formation", "Non-alcoholic fatty liver disease (NAFLD)", 
  #"Regulation of autophagy", "Inflammatory bowel disease (IBD)", 
  #"Chagas disease (American trypanosomiasis)", "Cell adhesion molecules (CAMs)", 
  #"Herpes simplex infection", "HTLV-I infection", "Hypertrophic cardiomyopathy (HCM)", 
  #"D-Glutamine and D-glutamate metabolism", "Prion diseases", "Arrhythmogenic right ventricular cardiomyopathy (ARVC)", 
  #"Mucin type O-Glycan biosynthesis", "Glycosphingolipid biosynthesis - globo series", 
  #"Butirosin and neomycin biosynthesis", "Synthesis and degradation of ketone bodies", 
  #"Citrate cycle (TCA cycle)", "Lysine biosynthesis", "Amyotrophic lateral sclerosis (ALS)", 
  #"Glycosylphosphatidylinositol(GPI)-anchor biosynthesis")
  mypaths <- mypaths[!is.na(mypaths$kegg),]
  
  #mypaths <- mypaths[-c(14,20),]
  mypaths$keggid = substr(mypaths$kegg,start = 1,stop = 8)
  mypaths$kegg_makenames <- paste0(make.names(mypaths$kegg),'1')
  mypaths <- mypaths[mypaths$p.adj<0.15,]
  
  
  
  res <- deg
  res$entrez = mapIds(org.Hs.eg.db,keys = res$hgnc_symbol,column = 'ENTREZID',keytype = 'SYMBOL',multiVals = 'first')
  res <- res[!is.na(res$entrez),]
  res$entrez <- as.character(res$entrez)
  for(i in 1:nrow(mypaths)){
    cat(i,"\t",mypaths$kegg[i],"\n")
    resa <- res[res$entrez %in% kegg.sets.hs[[mypaths$kegg[i]]], ]
    so <- AddModuleScore(so,features = list(resa$hgnc_symbol),
                         name=mypaths$kegg[i])
    rm(resa)
  }
  saveRDS(so@meta.data,file='Kaminski_AdelPathwayScores_metadata.rds')
  saveRDS(mypaths,file='/home/rstudio/external/AdelPeter/AdditionalPathways.rds')
  
  mypaths <- readRDS("AdditionalPathways.rds")
  foxo = res[res$entrez%in%unique(kegg.sets.hs[[grep('FoxO',names(kegg.sets.hs))]]),]
  tgfb = res[res$entrez%in%unique(kegg.sets.hs[[grep('TGF',names(kegg.sets.hs))]]),]
  datatable(foxo,rownames = F,extensions = 'Buttons',
            options = list(searching = TRUE,pageLength = 10,dom = 'Bfrtip', 
                           buttons = list(list(extend = 'csv',  
                                               filename = paste0("FoxO_signalling_genes"))),
                           scrollX = TRUE,scrollCollapse = TRUE))
  datatable(tgfb,rownames = F,extensions = 'Buttons',
            options = list(searching = TRUE,pageLength = 10,dom = 'Bfrtip', 
                           buttons = list(list(extend = 'csv',  
                                               filename = paste0("TGF-Beta_signalling_genes"))),
                           scrollX = TRUE,scrollCollapse = TRUE))
  
  
  ##--- log2FC
  res <- deg
  res$entrez = mapIds(org.Hs.eg.db,keys = res$hgnc_symbol,column = 'ENTREZID',keytype = 'SYMBOL',multiVals = 'first')
  res <- res[!is.na(res$entrez),]
  res$entrez <- as.character(res$entrez)
  #mypaths <- mypaths[-3,]
  ## save images for l2fc
  foldchages = res$logFC
  names(foldchages) = res$entrez
  for(i in 1:nrow(mypaths)){
    cat(i," ",mypaths$kegg[i],"\n")
    tmp=pathview(gene.data = foldchages,pathway.id = mypaths$keggid[i],
                 species = 'hsa',out.suffix='l2fc',
                 mid = list('#ffffcc','#ffffcc'),
                 low=list('#4575b4','#4575b4'),high = list('#e31a1c','#e31a1c'),
                 limit = list(0.5,0.5))
  }
  
  ##--- selected pathways, FDR
  foldchages = -log10(res$adj.P.Val)
  names(foldchages) = res$entrez
  for(i in 1:nrow(mypaths)){
    cat(i,"\n")
    tmp=pathview(gene.data = foldchages,pathway.id = mypaths$keggid[i],
                 species = 'hsa',out.suffix='fdr',
                 limit = list(5,5),mid = list('#ffffcc','#ffffcc'),
                 low=list('white','white'),high = list('#e31a1c','#e31a1c'))
  }
  
  ##--- selected pathways, LogFC, for genes with FDR<10%
  resa <- res[res$adj.P.Val<0.1,]
  foldchages = resa$logFC
  names(foldchages) = resa$entrez
  for(i in 1:nrow(mypaths)){
    cat(i,"\n")
    tmp=pathview(gene.data = foldchages,pathway.id = mypaths$keggid[i],
                 species = 'hsa',out.suffix='l2fc_statdiff',
                 mid = list('#ffffcc','#ffffcc'),
                 low=list('#4575b4','#4575b4'),high = list('#e31a1c','#e31a1c'),
                 limit = list(0.5,0.5))
  }
  
  
  
  

  
}

##-- automated plotting in the report
if(F){
  
  i=1
  
  
  length(unique(kegg.sets.hs[[mypaths$kegg[i] ]] ))
  nrow(res[res$entrez%in%unique(kegg.sets.hs[[mypaths$kegg[i] ]]),])
  
  
  myplotfun1 <- function(i=1,mypaths,so){
    P = FeaturePlot(so,features = mypaths$kegg_makenames[i], 
                label = T, repel = T,reduction = 'umap')+
      scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+
      ggtitle(gsub('hsa\\d*','',mypaths$kegg[i]))
    print(P)
  }
  
  myplotfun2 <- function(i=1,mypaths,so,plt=1){
    pl <- so@meta.data
    pl <- pl[!pl$cell_type %in% c('Aberrant_Basaloid','Ionocyte'),]
    mx = max(pl[,mypaths$kegg_makenames[i]])
    titl = gsub('hsa\\d*','',mypaths$kegg[i])
    subtitl = paste0("(",
      nrow(res[res$entrez%in%unique(kegg.sets.hs[[mypaths$kegg[i] ]]),]),
      "/",
      length(unique(kegg.sets.hs[[mypaths$kegg[i] ]] )),
      " genes annotated in the pathway are expressed)"
    )
    
    if(plt==1){
      P=ggplot(pl,aes(x=cell_type,y=pl[,mypaths$kegg_makenames[i]],fill=Disease))+
        geom_boxplot(outlier.shape = NA)+
        theme_bw()+
        stat_compare_means(label = "p.signif", method = "t.test",
                           hide.ns = T,label.y=seq(0.9*mx,1.5*mx,length.out=37)) + 
        ggtitle(label = titl,subtitle=subtitl)+
        ylab("")+xlab("Avg. normalized expression")+
        geom_hline(yintercept=median(pl[,mypaths$kegg_makenames[i]]),col='red')+
        scale_fill_manual(values=c('Control'='#91cf60','COPD'='#fc8d59'))+
        theme(axis.text=element_text(color='black'),
              axis.text.x=element_text(angle=90,hjust=1))
      print(P)
    }else{
      Q=ggplot(pl,aes(x=cell_type,y=pl[,mypaths$kegg_makenames[i]],fill=Disease))+
        geom_boxplot(outlier.shape = NA)+
        theme_bw()+
        ggtitle(label = titl,subtitle=subtitl)+
        ylab("")+xlab("Avg. normalized expression")+
        geom_hline(yintercept=median(pl[,mypaths$kegg_makenames[i]]),col='red')+
        scale_fill_manual(values=c('Control'='#91cf60','COPD'='#fc8d59'))+
        theme(axis.text=element_text(color='black'),
              axis.text.x=element_text(angle=90,hjust=1))
      print(Q)
    }
  }
  
  myplotfun2(i=2,mypaths,so,plt=2)

}


if(F){
}


