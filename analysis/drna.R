require(ggplot2)
require(Biostrings)
require(GenomicAlignments)
require(tidyverse)
require(gridExtra)

datadir='~/data/ncov/'
dbxdir='~/Dropbox/'

##devel from plot for timp ubiome
bam2df <- function(bam) {
    dat.raw=readGAlignments(bam, use.names=T)
    dat.gr=granges(dat.raw)
    chroms=seqlevels(dat.gr)
    perchrom.len=data.frame(chr=factor(as.character(seqnames(dat.gr)), levels=chroms), rlength=width(dat.gr))

    if (dim(perchrom.len)[1] < 1) {
        perchrom.len=rbind(perchrom.len, data.frame(chr='none', rlength=0))
    }
    return(perchrom.len)
}

samps=c('HP4', 'HP20')

alignlens=data.frame(chr=as.character(), rlength=as.numeric(), samp=as.character())
numreads=data.frame(samp=as.character(), numreads=as.numeric())
fullcovs=data.frame(samp=as.character(), cov=as.numeric(), pos=as.numeric())

for (i in samps) {
    bamfile=paste0(datadir, i, '/align/', i, '.sorted.bam')
    alignlen=bam2df(bamfile)
    alignlen$samp=i

    numalign=dim(alignlen)[1]
    numreads=rbind(numreads, data.frame(samp=i, numreads=numalign))
    
    alignlens=rbind(alignlen, alignlens)

    align=readGAlignments(bamfile, use.names=T)
    cov=as.numeric(coverage(align, drop.D.ranges=TRUE)[[1]])
    pos=seq(1, length(cov), 1)
    fullcov=data.frame(samp=i, cov, pos)
}


