#===============================================================================
#         FILE: iclip_seq step6
#
#  DESCRIPTION: count distribution of clip reads in autr vs. cutr which are two 3'utr regions defined by different pA sites
#
#
#       AUTHOR: Wei Wang  wwei320@gmail.com
#       
#      VERSION: 2.0
#      CREATED: 2019-Nov-01
#===============================================================================

########to count distribution of clip reads in autr vs. cutr 
library(tidyr)
library(ggplot2)
library(gridExtra)
require(VennDiagram)
library(dplyr)
require(GenomicRanges)
library(AnnotationDbi)


####################functions for cUTR ##########################################
####################functions for non-coding (optional), 3utr, 5utr, intron, intergenic (others), cds, cutr, autr ##########################################
create_noncoding_from_pAs <- function(pA.df,geno){
  if(geno == "mm9"){
  ncTx = read.csv("/scratch/ww346/analyze/Structure/dinghai/scriptfromDH/annotation/mm9_ncTx.csv", as.is = T)
  }else if (geno == "hg19"){
  ncTx = read.csv("/scratch/ww346/analyze/Structure/dinghai/scriptfromDH/annotation/hg19_ncTx.csv", as.is = T)
  }
  
  ncTx = GRanges(seqnames = ncTx$chrom, strand = ncTx$strand,
               ranges = IRanges(start = ncTx$txStart, end = ncTx$txEnd),
               gene_symbol = ncTx$name) #ncTx$name?

# remove duplicates
  ncTx = unique(ncTx)
  names(ncTx)=ncTx$gene_symbol
  ncTx               
}


create_5utr_from_pAs <- function(pA.df,geno){
  txdb = loadDb(paste0("/scratch/ww346/analyze/annotation/ucsc.txdb/",geno,".refGene.txdb.sqlite")) 
  if (geno=='mm9'){
   require(org.Mm.eg.db)
  IDDB <- org.Mm.eg.db 
  } else if (geno=='hg19'){
    require(org.Hs.eg.db)
    IDDB <- org.Hs.eg.db 
  }
  
  fiveUTRs = fiveUTRsByTranscript(txdb, use.names=T) 
  x = unlist(fiveUTRs)
  names(x) = mapIds(IDDB, keys = names(x), keytype =  "ACCNUM", column =  "SYMBOL")

  fiveUTRs = split(x, names(x))
  fiveUTRs = reduce(fiveUTRs)

  # refdir="/scratch/ww346/analyze/APAanalyzer/REF/"
  # refdf= read.table(paste0(refdir,"SRS.",geno,".conserved.FL.txt"),header=TRUE)

  # df = refdf%>%dplyr::select(gene_symbol)
  # dfsymbol = unique(df$gene_symbol)
  # refsymbol = names(fiveUTRs)
  # mergedsymbol<- intersect(dfsymbol, refsymbol) 
  # gr=unlist(fiveUTRs[mergedsymbol])
  gr=unlist(fiveUTRs)                    
  gr$gene_symbol = names(gr)
  gr
}

create_3utr_from_pAs <- function(pA.df,geno){
  txdb = loadDb(paste0("/scratch/ww346/analyze/annotation/ucsc.txdb/",geno,".refGene.txdb.sqlite")) 
  if (geno=='mm9'){
   require(org.Mm.eg.db)
  IDDB <- org.Mm.eg.db 
  } else if (geno=='hg19'){
    require(org.Hs.eg.db)
    IDDB <- org.Hs.eg.db 
  }
  
  threeUTRs = threeUTRsByTranscript(txdb, use.names=T) 
  x = unlist(threeUTRs)
  names(x) = mapIds(IDDB, keys = names(x), keytype =  "ACCNUM", column =  "SYMBOL")

  threeUTRs = split(x, names(x))
  threeUTRs = reduce(threeUTRs)

  
  gr=unlist(threeUTRs)                    
  gr$gene_symbol = names(gr)
  gr
}


create_intron_from_pAs <- function(pA.df,geno){
  txdb = loadDb(paste0("/scratch/ww346/analyze/annotation/ucsc.txdb/",geno,".refGene.txdb.sqlite")) 
  if (geno=='mm9'){
   require(org.Mm.eg.db)
  IDDB <- org.Mm.eg.db 
  } else if (geno=='hg19'){
    require(org.Hs.eg.db)
    IDDB <- org.Hs.eg.db 
  }
  
  introns = intronsByTranscript(txdb, use.names=T)
  x = unlist(introns)
  names(x) = mapIds(IDDB, keys = names(x), keytype =  "ACCNUM", column =  "SYMBOL")

  introns = split(x, names(x))
  introns = reduce(introns)

  
  gr=unlist(introns)                    
  gr$gene_symbol = names(gr)
  gr
}


create_cds_from_pAs <- function(pA.df,geno){
  txdb = loadDb(paste0("/scratch/ww346/analyze/annotation/ucsc.txdb/",geno,".refGene.txdb.sqlite")) 
  if (geno=='mm9'){
   require(org.Mm.eg.db)
  IDDB <- org.Mm.eg.db 
  } else if (geno=='hg19'){
    require(org.Hs.eg.db)
    IDDB <- org.Hs.eg.db 
  }
  
  
  CDSbygene = cdsBy(txdb, by=("gene"))
  x = unlist(CDSbygene)
  acc2sym = AnnotationDbi::select(IDDB, keys = names(x), keytype =  "ENTREZID", columns =  "SYMBOL")
  acc2sym[is.na(acc2sym$SYMBOL),]$SYMBOL = acc2sym[is.na(acc2sym$SYMBOL), ]$ENTREZID
  names(x) = acc2sym$SYMBOL
  CDSbygene = split(x, names(x))
  CDSbygene = reduce(CDSbygene)

  refdir="/scratch/ww346/analyze/APAanalyzer/REF/"
  refdf= read.table(paste0(refdir,"SRS.",geno,".conserved.FL.txt"),header=TRUE)

  df = refdf%>%dplyr::select(gene_symbol)
  dfsymbol = unique(df$gene_symbol)
  refsymbol = names(CDSbygene)
  mergedsymbol<- intersect(dfsymbol, refsymbol) 
  gr=unlist(CDSbygene[mergedsymbol])
  gr$gene_symbol = names(gr)
  gr
}

create_cUTR_from_pAs = function(pA.df,geno){
  refdir="/scratch/ww346/analyze/APAanalyzer/REF/"
  df= read.table(paste0(refdir,"SRS.",geno,".conserved.FL.txt"),header=TRUE)
  df$strand = df$Strand
  df$chr = df$Chrom
  df$start = as.numeric(ifelse(df$strand == "+", df$cdsend + 1, df$First))
  df$end = as.numeric(ifelse(df$strand == "+", df$First, df$cdsend - 1))   # DINGHAI'S CODE need to be checked by dinghai
  #df$end = as.numeric(ifelse(df$strand == "+", df$pos, df$cds_end - 1))   ###added by ww
  df = df[df$end >= df$start, ]
  df $ length <- df$end - df$start
  df = subset(df,select=-c(Strand, Chrom, First, Last))
  #df[, c("pos", "cds_end")] = NULL     ###added by ww
  
  makeGRangesFromDataFrame(df, keep.extra.columns=T)
}

####################functions for aUTR##########################################

create_aUTR_from_pAs = function(pA.df,geno){
  refdir="/scratch/ww346/analyze/APAanalyzer/REF/"
  df= read.table(paste0(refdir,"SRS.",geno,".conserved.FL.txt"),header=TRUE)
  df$strand = df$Strand
  df$chr = df$Chrom
  df$start = as.numeric(ifelse(df$strand == "+", df$First + 1, df$Last))
  df$end = as.numeric(ifelse(df$strand == "+", df$Last, df$First - 1))   # DINGHAI'S CODE need to be checked by dinghai
  #df$end = as.numeric(ifelse(df$strand == "+", df$pos, df$cds_end - 1))   ###added by ww
  df = df[df$end >= df$start, ]
  df $ length <- df$end - df$start
  df = subset(df,select=-c(Strand, Chrom, First, Last))
  
  makeGRangesFromDataFrame(df, keep.extra.columns=T)
}


#######distribution of #reads with no clustering
add_hg19_iCLIP_Stau1 = function(pA.df){
  library(rtracklayer)
  # list.files("~/projects/fud/RIPiT/")
  Stau1_p = import("/scratch/ww346/analyze/Str_data/Stau1_iCLIP/rawout/HEK293_count_sorted.p.bw", format = "BigWig") # hg19
  # summary(stau1_p$score)
  # hist(log10(stau1_p$score), breaks=50)
  strand(Stau1_p)='+'

  Stau1_n = import("/scratch/ww346/analyze/Str_data/Stau1_iCLIP/rawout/HEK293_count_sorted.n.bw", format = "BigWig") # hg19
  strand(Stau1_n)='-'

  Stau1.hg19 = c(Stau1_p, Stau1_n)

  #####calculate total count for rpm calculation

  total=sum(Stau1.hg19$score)
  noncoding = create_noncoding_from_pAs(pA.df,geno='hg19')
  utr5 = create_5utr_from_pAs(pA.df,geno='hg19')
  utr3 = create_3utr_from_pAs(pA.df,geno='hg19')

  intron=create_intron_from_pAs(pA.df,geno='hg19')

  CDS =create_cds_from_pAs(pA.df,geno='hg19')

  cUTR = create_cUTR_from_pAs(pA.df,geno='hg19')
  aUTR = create_aUTR_from_pAs(pA.df,geno='hg19')



  #################################3utr region
  ncread= as.data.frame(mergeByOverlaps(noncoding, Stau1.hg19)[, c("gene_symbol","score")]) %>%
    summarise(Stau1_iCLIP_count_noncoding = sum(score))

  utr5read= as.data.frame(mergeByOverlaps(utr5, Stau1.hg19)[, c("gene_symbol","score")]) %>%
    summarise(Stau1_iCLIP_count_5utr = sum(score))

  utr3read= as.data.frame(mergeByOverlaps(utr3, Stau1.hg19)[, c("gene_symbol","score")]) %>%
    summarise(Stau1_iCLIP_count_3utr = sum(score))

  intronread= as.data.frame(mergeByOverlaps(intron, Stau1.hg19)[, c("gene_symbol","score")]) %>%
    summarise(Stau1_iCLIP_count_intron = sum(score))

  CDSread= as.data.frame(mergeByOverlaps(CDS, Stau1.hg19)[, c("gene_symbol","score")]) %>%
    summarise(Stau1_iCLIP_count_cds = sum(score))


  cUTRread= as.data.frame(mergeByOverlaps(cUTR, Stau1.hg19)[, c("gene_symbol","score")]) %>%
    summarise(Stau1_iCLIP_count_cutr = sum(score))
  
  aUTRread= as.data.frame(mergeByOverlaps(aUTR, Stau1.hg19)[, c("gene_symbol","score")]) %>%
    summarise(Stau1_iCLIP_count_autr = sum(score))


  interread=total-utr5read[1,]-utr3read[1,]-intronread[1,]-CDSread[1,]-ncread[1,]

disdf = data.frame(region=c('5utr','CDS','3utr','intron','nc','intergenic','cUTR','aUTR'), 
  order=c(1,2,3,4,5,6,7,8),
  read_dis=c(utr5read[1,],CDSread[1,],utr3read[1,],intronread[1,],ncread[1,],interread,cUTRread[1,], aUTRread[1,]))


disdf
}



#######distribution of #reads with no clustering
cluster_hg19_iCLIP_Stau1 = function(pA.df){
  library(rtracklayer)
  # list.files("~/projects/fud/RIPiT/")
  setwd('/scratch/ww346/analyze/Str_data/Stau1_iCLIP/rawsam')
  stau1= read.csv('clusters30.unique.reads.csv', as.is=T)

  Stau1= stau1%>% 
      dplyr::rename(chr=chromosome, 
        start=position,
        score=ERR605258)%>%
      mutate(end=start)

  Stau1.hg19 = makeGRangesFromDataFrame(Stau1, keep.extra.columns=T)

  #####calculate total count for rpm calculation

  total=sum(Stau1.hg19$score)

  noncoding = create_noncoding_from_pAs(pA.df,geno='mm9')
  utr5 = create_5utr_from_pAs(pA.df,geno='mm9')
  utr3 = create_3utr_from_pAs(pA.df,geno='mm9')

  intron=create_intron_from_pAs(pA.df,geno='mm9')
  CDS =create_cds_from_pAs(pA.df,geno='hg19')

  cUTR = create_cUTR_from_pAs(pA.df,geno='hg19')
  aUTR = create_aUTR_from_pAs(pA.df,geno='hg19')


  #################################3utr region
  ncread= as.data.frame(mergeByOverlaps(noncoding, Stau1.hg19)[, c("gene_symbol","score")]) %>%
    summarise(Stau1_iCLIP_cluster_noncoding = nrow(.),
      Stau1_iCLIP_count_noncoding = sum(score))

  utr5read= as.data.frame(mergeByOverlaps(utr5, Stau1.hg19)[, c("gene_symbol","score")]) %>%
    summarise(Stau1_iCLIP_cluster_5utr = nrow(.),
      Stau1_iCLIP_count_5utr = sum(score))

  utr3read= as.data.frame(mergeByOverlaps(utr3, Stau1.hg19)[, c("gene_symbol","score")]) %>%
    summarise(Stau1_iCLIP_cluster_3utr = nrow(.),
      Stau1_iCLIP_count_3utr = sum(score))

  intronread= as.data.frame(mergeByOverlaps(intron,Stau1.hg19)[, c("gene_symbol","score")]) %>%
    summarise(Stau1_iCLIP_cluster_intron = nrow(.),
      Stau1_iCLIP_count_intron = sum(score))

  CDSread= as.data.frame(mergeByOverlaps(CDS, Stau1.hg19)[, c("gene_symbol","score")]) %>%
    summarise(Stau1_iCLIP_cluster_cds = nrow(.),
      Stau1_iCLIP_count_cds = sum(score))

  
  cUTRread= as.data.frame(mergeByOverlaps(cUTR, Stau1.hg19)[, c("gene_symbol","score")]) %>%
    summarise(Stau1_iCLIP_cluster_cUTR = nrow(.),
      Stau1_iCLIP_count_cutr = sum(score))
  
  

  aUTRread= as.data.frame(mergeByOverlaps(aUTR, Stau1.hg19)[, c("gene_symbol","score")]) %>%
    summarise(Stau1_iCLIP_cluster_aUTR = nrow(.),
      Stau1_iCLIP_count_autr = sum(score))


  interread=total-utr5read[,2]-utr3read[,2]-intronread[,2]-CDSread[,2]-ncread[,2]
  intercluster=length(Stau1.hg19)-utr5read[,1]-utr3read[,1]-intronread[,1]-CDSread[,1]-ncread[,1]

disdf = data.frame(region=c('5utr','CDS','3utr','intron','nc','intergenic','cUTR','aUTR'), 
    order=c(1,2,3,4,5,6,7,8),
    cluster_dis=c(utr5read[,1],CDSread[,1],utr3read[,1],intronread[,1],ncread[,1],intercluster,cUTRread[,1], aUTRread[,1]),
    read_dis=c(utr5read[,2],CDSread[,1],utr3read[,2],intronread[,2],ncread[,2],interread,cUTRread[,2], aUTRread[,2]))


disdf
}
###12/18/19 add stau1_clip data in stau1 rip data set


i=2

######use ref seq table in APAlyzer  12/29/19
pA.df=1


########plot distribution of unclustered #reads 

max_RPM_df=add_hg19_iCLIP_Stau1(pA.df)
print (max_RPM_df)

setwd('/scratch/ww346/analyze/Str_data/Stau1_iCLIP/plot')

p1 <- ggplot(data=max_RPM_df, aes(x=reorder(region,order),y=read_dis)) + 
  geom_bar(stat="identity")+
  #scale_y_continuous(limits = c(0, 1))+
  ggtitle(paste0("Stau1_iCLIP_distribution")) + labs(x="iCLIP reads distribution", y="Sum of iCLIP reads")+
  # geom_hline(yintercept=0, color='black')+
  # geom_vline(xintercept=0, color='black')+
  # theme(axis.title.x = element_blank(),
  #       panel.background = element_blank())
  theme_classic()

  jpeg(paste0('Stau1_iCLIP_reads_distribution.tiff'), res=300, width=2000, height=2000, type="cairo")

  p1
dev.off()




########plot distribution of cluster (30nt window) 

max_RPM_df=cluster_hg19_iCLIP_Stau1(pA.df)
print (max_RPM_df)

setwd('/scratch/ww346/analyze/Str_data/Stau1_iCLIP/plot')

p1 <- ggplot(data=max_RPM_df, aes(x=reorder(region,order),y=cluster_dis)) + 
  geom_bar(stat="identity")+
  #scale_y_continuous(limits = c(0, 1))+
  ggtitle(paste0("Stau1_iCLIP_distribution_30ntclusters")) + labs(x="Stau1 iCLIP site distribution", y="Sum of Stau1 binding site")+
  # geom_hline(yintercept=0, color='black')+
  # geom_vline(xintercept=0, color='black')+
  # theme(axis.title.x = element_blank(),
  #       panel.background = element_blank())
  theme_classic()

  jpeg(paste0('Stau1_iCLIP_clusters30_distribution.tiff'), res=300, width=2000, height=2000, type="cairo")

  p1
dev.off()



