library(DESeq2)
library(sqldf)
library(stringr)
library(gplots)
library(RColorBrewer)
library(ggVennDiagram)
library(ggplot2)
library(reshape2)
library(cowplot)
library(patchwork)
library(enrichR)
library(gridExtra)


################################################################################
############################# Read rnaseq data #################################
################################################################################

########### Count table
cnt <- read.table("rnaseqdata/cnt_clean.csv")


########### Condition data
cond <- read.csv("rnaseqdata/cond.csv")
cond$raw_sample_name <- str_replace_all(cond$raw_sample_name, "-","_")
cond$raw_sample_name[cond$mod_name!=""] <- cond$mod_name[cond$mod_name!=""]
rownames(cond) <- cond$raw_sample_name
cond <- cond[colnames(cnt),]

cond$conc_time <- sprintf("%s_%s",cond$h2s_conc,cond$time)
cond$conc_time_rep <- sprintf("%s_%s_%s",cond$h2s_conc,cond$time,cond$replicate)

write.csv(cond[,c("treatment","h2s_conc","time","replicate","conc_time","conc_time_rep")],
          "rnaseqdata/cond_clean.csv")

#### Only keep h2s samples; remove anoxia samples
cond$h2s_conc[cond$h2s_conc=="con"] <- "0"  
keepsample <- cond$treatment!="anoxia"
cond <- cond[keepsample,]
cnt <- cnt[,keepsample]     

# Problematic samples
keepsample <- !(cond$sample_name %in% c("S_50ppm_2h_5","S_50ppm_2h_6"))
cond <- cond[keepsample,]
cnt <- cnt[,keepsample]     

################################################################################
######################### Rename to gene symbol ################################
################################################################################


gid <- read.csv("c_elegans.PRJNA13758.current.geneIDs.txt", header = FALSE)
gid <- gid[gid[,5]=="Live",c(2,3,4,6)]
colnames(gid) <- c("Gene.stable.ID","genesym","transcript","genetype")
gid_entrez <- read.csv("mart_export.txt.gz",sep="\t")[,c(2,3)] #to also get NCBI gene id
gid <- unique(merge(gid_entrez,gid))
gid

# use genesym, otherwise gene stable id, otherwise wormbaseid, for counts
gid <- gid[gid$Gene.stable.ID %in% rownames(cnt),]
gid <- gid[!duplicated(gid$transcript),]
rownames(gid) <- gid$Gene.stable.ID
gid$genesym[gid$genesym==""] <- gid$transcript[gid$genesym==""]

gid <- gid[!duplicated(gid$genesym),]

newnames <- gid[rownames(cnt),]$genesym
newnames[is.na(newnames)] <- rownames(cnt)[is.na(newnames)]
df_name_mapping <- data.frame(wbid=rownames(cnt), genesym=newnames)
rownames(cnt) <- newnames




################################################################################
########################## Clustering ##########################################
################################################################################


if(FALSE){
  dds <- DESeqDataSetFromMatrix(countData = cnt, colData = cond, design = ~1)
  dds <- DESeq(dds)
  
  ncnt <- counts(dds, normalize=TRUE)
  logcnt <- log(1+ncnt)
  logcnt2 <- logcnt
  colnames(logcnt2) <- cond$sample_name
  plot(hclust(dist(t(logcnt2)))) # could keep the genes with most variation too
}



################################################################################
################### Volcano plots and GO analysis ##############################
################################################################################


setEnrichrSite("WormEnrichr")
listEnrichrDbs()
dbs <- c("GO_Biological_Process_2018")


dds <- DESeqDataSetFromMatrix(countData = cnt, colData = cond, design = ~sample_group)
dds <- DESeq(dds)
allgo <- NULL
for(curtime in c(1,2,12)){
  for(curconc in c(50,150)){
    g1 <- sprintf("S_%sppm_%sh", curconc, curtime)
    g2 <- sprintf("S_con_%sh", curtime)
    
    
    print(paste(g1,"vs",g2))  

    ### DE testing
    res <- as.data.frame(results(dds, contrast = c("sample_group", g1, g2)))
    res <- res[!is.na(res$pvalue),]
    res <- res[order(res$pvalue, decreasing = FALSE),]

    write.csv(res,sprintf("newout/de/%s vs %s.csv",g1,g2))
    
    ### Volcano plot
    res$log_padj <- -log10(res$padj)

    ################### First volcano, long list of genes ######################
    if(TRUE){
      highlight_genes <-   c("sqrd-1","ethe-1","gst-19","gst-4","gst-30","ftn-1","semo-1","cysl-1","cysl-2","clec-66","clec-4","hsp-70","rhy-1",
                             "hsp-16.41","hsp-16.2","nhr-57","mpst-3","mpst-5 and sysm-1","irg-1","irg-2")
      
      res$txt <- ""
      res$txt[rownames(res) %in% highlight_genes] <- rownames(res)[rownames(res) %in% highlight_genes]
      
      res$color <- "gray"
      res$color[!is.na(res$padj) & res$padj<1e-10 & res$log2FoldChange<0] <- "lightblue"
      res$color[!is.na(res$padj) & res$padj<1e-10 & res$log2FoldChange>0] <- "lightgreen"
      res$color[rownames(res) %in% highlight_genes] <- "red"
      res <- res[order(res$color),]
      
      print(paste("#down",sum(res$color=="lightblue"),"#up",sum(res$color=="lightgreen")))
      
      ggplot(res, aes(x=log2FoldChange, y=log_padj, label=txt)) + 
        geom_point(color=res$color) + 
        geom_text() + 
        xlab(sprintf("Log2FC %s vs %s",g1, g2)) +
        ylab("-Log10 p.adj") + 
        xlim(c(-10,10)) 
      ggsave(sprintf("newout/firstvolcano/%s vs %s.pdf",g1,g2))      
    }

    ############################# Second volcano, 3 genes ######################
    if(TRUE){
      highlight_genes <- c("ftn-1", "ftn-2","smf-3")
      
      res$txt <- ""
      res$txt[rownames(res) %in% highlight_genes] <- rownames(res)[rownames(res) %in% highlight_genes]
      
      res$color <- "gray"
      res$color[!is.na(res$padj) & res$padj<1e-10 & res$log2FoldChange<0] <- "lightblue"
      res$color[!is.na(res$padj) & res$padj<1e-10 & res$log2FoldChange>0] <- "lightgreen"
      res$color[rownames(res) %in% highlight_genes] <- "red"
      res <- res[order(res$color),]
      
      ggplot(res, aes(x=log2FoldChange, y=log_padj, label=txt)) + 
        geom_point(color=res$color) + 
        geom_text() + 
        xlab(sprintf("Log2FC %s vs %s",g1, g2)) +
        ylab("-Log10 p.adj") + 
        xlim(c(-10,10)) 
      ggsave(sprintf("newout/secondvolcano/%s vs %s.pdf",g1,g2))      
    }
    
    ### GO analysis
    glist <- as.character(na.omit(rownames(res)[res$padj<1e-10]))  
    enriched <- enrichr(glist, dbs)
    
    for(cur_db in dbs){
      df <- enriched[[cur_db]]
      df <- df[df$P.value<0.01,]
      df <- df[order(df$P.value,decreasing = TRUE),]
      df$logp <- -log10(df$P.value)
      
      df$Term <- str_trim(df$Term)
      df <- df[order(df$logp,decreasing = FALSE),]
      df$Term<-factor(df$Term, levels = unique(df$Term))
      
      p <- ggplot(df, aes(x = Term, y = logp)) + 
        geom_bar(stat = "identity") + coord_flip() + 
        xlab("")+
        ylab("-Log10(p-value)")+
        theme_bw()
      p
      ggsave(plot=p, sprintf("newout/go.all/go %s vs %s - %s.pdf",g1,g2, cur_db))
      
      df$g1 <- g1
      df$g2 <- g2
      allgo <- rbind(allgo, df[,c("g1","g2","Term","logp")])
    }
  }
}
write.csv(allgo,"newout/go.all/golist.csv")


################################################################################
################### GO plotting ################################################         
################################################################################

### Use curated terms
allgo <- read.csv("newout/go.curated/curatedgo.csv")         
for(g1 in unique(allgo$g1)){
  for(g2 in unique(allgo$g2)){
    p <- ggplot(df, aes(x = Term, y = logp)) + 
      geom_bar(stat = "identity") + coord_flip() + 
      xlab("")+
      ylab("-Log10(p-value)")+
      theme_bw()
    p
    ggsave(plot=p, sprintf("newout/go.curated/go %s vs %s.pdf",g1,g2))
  }
}


################################################################################
################### Heatmaps of some genes #####################################  
################################################################################

allde <- NULL
for(curtime in c(1,2,12)){
  for(curconc in c(50,150)){
    g1 <- sprintf("S_%sppm_%sh", curconc, curtime)
    g2 <- sprintf("S_con_%sh", curtime)
    
    
    print(paste(g1,"vs",g2))  
    
    ### DE testing
    res <- as.data.frame(results(dds, contrast = c("sample_group", g1, g2)))
    res <- res[!is.na(res$pvalue),]
    res <- res[order(res$pvalue, decreasing = FALSE),]
    res$geneid <- rownames(res)
    res$time <- curtime
    res$conc <- curconc
    allde <- rbind(allde,res)
  }
}




plot_heatmap_highlight_genes <- function(highlight_genes){
  allde_sub <- allde[allde$geneid %in% highlight_genes,c("geneid","log2FoldChange","time","conc")]
  allde_sub$h_time <- sprintf("%sh",allde_sub$time)
  allde_sub$h_time <- factor(allde_sub$h_time, levels=c("1h","2h","12h"))
  
  forconc <- allde_sub[allde_sub$conc==50,c("geneid","log2FoldChange","h_time")]
  p1 <- ggplot(forconc, aes(h_time, geneid, fill=log2FoldChange)) + geom_tile(color = "black") +
    scale_fill_gradient2(low = "#075AFF",mid = "#FFFFFF",high = "#FF0000")+ xlab("Time, 50ppm")+ylab("")
  
  forconc <- allde_sub[allde_sub$conc==150,c("geneid","log2FoldChange","h_time")]
  p2 <- ggplot(forconc, aes(h_time, geneid, fill=log2FoldChange)) + geom_tile(color = "black") +
    scale_fill_gradient2(low = "#075AFF",mid = "#FFFFFF",high = "#FF0000")+ xlab("Time, 150ppm")+ylab("")
  
  ptot <- p1|p2  
}



highlight_genes <- unique(allde$geneid[allde$padj<1e-40])
ptot <- plot_heatmap_highlight_genes(highlight_genes)
ptot
ggsave("newout/heatmap/heatmap_tc_top.pdf",ptot, height = 10)




################################################################################
################### Venn diagram ############################################### 
################################################################################

read_one_de_list <- function(fname){
  delist <- read.csv(fname)
  delist <- delist[!is.na(delist$padj),]
  delist <- delist[order(delist$padj),]
  delist$X[!is.na(delist$padj) & delist$padj<1e-10]
}

list4venn <- list(
  c50vsC1h  = read_one_de_list("out.de.3/S_50ppm_1h vs S_con_1h.csv"),
  c50vsC2h  = read_one_de_list("out.de.3/S_50ppm_2h vs S_con_2h.csv"),
  c50vsC12h = read_one_de_list("out.de.3/S_50ppm_12h vs S_con_12h.csv")
)
ggVennDiagram(list4venn, label_alpha = 0)
ggsave(sprintf("newout/venn/123.pdf"))

list4venn <- list(
  c150vsC1h  = read_one_de_list("out.de.3/S_150ppm_1h vs S_con_1h.csv"),
  c150vsC2h  = read_one_de_list("out.de.3/S_150ppm_2h vs S_con_2h.csv"),
  c150vsC12h = read_one_de_list("out.de.3/S_150ppm_12h vs S_con_12h.csv")
)
ggVennDiagram(list4venn, label_alpha = 0)
ggsave(sprintf("newout/venn/456.pdf"))


list4venn <- list(
  c50vsC1h  = read_one_de_list("out.de.3/S_50ppm_1h vs S_con_1h.csv"),
  c150vsC1h  = read_one_de_list("out.de.3/S_150ppm_1h vs S_con_1h.csv")
)
ggVennDiagram(list4venn, label_alpha = 0)
ggsave(sprintf("newout/venn/12.pdf"))


list4venn <- list(
  c50vsC2h  = read_one_de_list("out.de.3/S_50ppm_2h vs S_con_2h.csv"),
  c150vsC2h  = read_one_de_list("out.de.3/S_150ppm_2h vs S_con_2h.csv")
)
ggVennDiagram(list4venn, label_alpha = 0)
ggsave(sprintf("newout/venn/34.pdf"))


list4venn <- list(
  c50vsC12h  = read_one_de_list("out.de.3/S_50ppm_12h vs S_con_12h.csv"),
  c150vsC12h  = read_one_de_list("out.de.3/S_150ppm_12h vs S_con_12h.csv")
)
ggVennDiagram(list4venn, label_alpha = 0)
ggsave(sprintf("newout/venn/56.pdf"))



###########################################################################################
################### TF analysis ###########################################################
###########################################################################################


targets_hif1 <- read.csv("tftarget/hif1.csv")$target
targets_skn1 <- read.csv("tftarget/skn1.csv")$target


#### DE analysis
dds <- DESeqDataSetFromMatrix(countData = cnt, colData = cond, design = ~sample_group)
dds <- DESeq(dds)
delist <- NULL
for(curconc in c(50,150)){
  for(curtime in c(1,2,12)){
    g1 <- sprintf("S_%sppm_%sh", curconc, curtime)
    g2 <- sprintf("S_con_%sh", curtime)
    print(paste(g1,"vs",g2))  
    
    ### DE testing
    res <- as.data.frame(results(dds, contrast = c("sample_group", g1, g2)))
    res <- res[!is.na(res$pvalue),]
    res <- res[order(res$pvalue, decreasing = FALSE),]
    res$log_padj <- -log10(res$padj)
    
    res$gene <- rownames(res)
    res$conc <- curconc
    res$time <- curtime
    delist <- rbind(delist, res)
  }
}


### 
show_target_confidence <- function(targets_forgene,tfname){
  outstat <- NULL
  for(curconc in c(50,150)){
    for(curtime in c(1,2,12)){
      de <- delist[delist$conc==curconc & delist$time==curtime,]
      de$istarget <- de$gene %in% targets_forgene
      de$log_pvalue <- log10(de$pvalue)
      onestat <- t.test(
        de$log_pvalue[de$istarget],
        de$log_pvalue[!de$istarget]
      )
      genecnt <- sum(!is.na(de$padj) & de$padj<1e-3 & de$istarget)
      thelab <- paste(
        round(log10(onestat$p.value), digits = 2), 
        " (",genecnt,")", sep="")
      outstat <- rbind(outstat, 
                       data.frame(
                         p=onestat$p.value, 
                         time=curtime, 
                         lab=thelab,
                         genecnt=genecnt,
                         conc=curconc))
    }
  }
  outstat$log_p <- log10(outstat$p)
  outstat$time_h <- factor(sprintf("%sh", outstat$time), levels=c("1h","2h","12h"))

  ggplot(outstat,  aes(time_h, genecnt, fill=paste(conc))) + 
    geom_bar(stat="identity", position=position_dodge()) + 
    labs(title=tfname, x="Time", y="Gene count") + labs(fill = "Conc (ppm)")
}

show_target_confidence(targets_hif1, "Hif-1")/
  show_target_confidence(targets_skn1, "Skn-1") 
ggsave("newout/tftarget/barplot_tf_time.pdf", width = 3)

