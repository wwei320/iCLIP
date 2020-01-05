rm(list = ls())
suppressMessages(library("GenomicAlignments"))
suppressMessages(library("GenomicRanges"))
suppressMessages(library("Rsamtools"))
suppressMessages(library("rtracklayer"))
suppressMessages(library("BiocParallel"))
suppressMessages(library("matrixStats"))
suppressMessages(library("magrittr"))
suppressMessages(library("readr"))
suppressMessages(library("dplyr"))
register(MulticoreParam(workers = 24))
args = commandArgs(trailingOnly=TRUE)

setwd(args[1])    ###  '/scratch/ww346/analyze/PPP/GSE100007_try/'
geno = args[2]  ## hg, rn, dm    ###hg19
bamdir =args[3]    ##  '/scratch/ww346/analyze/PPP/GSE100007_try/rawsam/'
REFdir =args[4]     ##  '/scratch/ww346/analyze/APAanalyzer/REF'
CPUS =args[5]      ##  '24'
samplename=args[6]


CPUS =as.integer(CPUS)
print (paste0('----------------',CPUS, '---------------------'))
tails="conserved"

###read strand###
mystring <- read_file(paste0(getwd(),"/test/strand.txt"))
mystring=gsub("\n","",mystring)

### bam files from star results
#fls = args[6]
setwd(bamdir)
gr <- import(paste0(samplename, ".bed"),format="BED") ##added 'asMates' for pair-end

grplus <- gr[strand(gr)=='+']
grminus <- gr[strand(gr)=='-']

end(grplus) <- start(grplus)
grplusid<-unique(grplus)
grplusid$count=countOverlaps(grplusid,grplus,type='equal')


start(grminus) <- end(grminus)
grminusid <- unique(grminus)
grminusid$count=countOverlaps(grminusid,grminus,type='equal')

grid <- c(grplusid,grminusid)

grid$rpm=grid$count/sum(grid$count)*1000000

grdf <- as.data.frame(grid)

#write.csv(grdf,paste0(samplename,'_count_rpm.csv'),sep='\t',row.names=FALSE)


# grdf1<-grdf%>%filter(strand=='+')%>%select(seqnames,start,end,rpm)%>%mutate(end=end+1)

# options(scipen=999)
# write.table(grdf1,paste0(samplename,'_rpm.p.bedGraph'),sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)

# grdf2<-grdf%>%filter(strand=='-')%>%select(seqnames,start,end,rpm)%>%mutate(start=start-1)

# write.table(grdf2,paste0(samplename,'_rpm.m.bedGraph'),sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)


#########use count instead of rpm to make count-based bedgraphs 

grdf1<-grdf%>%filter(strand=='+')%>%select(seqnames,start,end,count)%>%mutate(end=end+1)

options(scipen=999)
write.table(grdf1,paste0(samplename,'_count.p.bedGraph'),sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)

grdf2<-grdf%>%filter(strand=='-')%>%select(seqnames,start,end,count)%>%mutate(start=start-1)

write.table(grdf2,paste0(samplename,'_count.m.bedGraph'),sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)






