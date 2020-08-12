library(tidyverse)

##dbxdir='~/Dropbox/timplab_data/ncov/nextera/'
dbxdir='~/Dropbox/timplab_data/ncov/nextera_v2/'

##datadir='/uru/Data/Nanopore/projects/ncov/nextera/'
datadir='/uru/Data/Nanopore/projects/ncov/nextera_v2/'

##covdir=paste0(datadir, 'cov/trimmed/')
covdir=paste0(datadir, 'cov/untrimmed/')
covfiles=list.files(covdir)

allcover=tibble(
    chr=as.character(),
    pos=as.numeric(),
    cov=as.numeric(),
    sample=as.character(),
    set=as.character())


##pdf(paste0(dbxdir, 'genome_cov_trimmed.pdf'), h=6, w=16)
for (i in covfiles) {
    covfile=paste0(covdir, i)
    
    sampname=strsplit(i, '_', fixed=TRUE)[[1]][1]
    sampset=paste0(strsplit(sampname, '-', fixed=TRUE)[[1]][1])

    if (endsWith(sampname, 'norm')) {
        cond='n'
    }else{
        cond=''
    }
    
    cover=read_tsv(covfile, col_names=c('chr', 'pos', 'cov')) %>%
        mutate(sample=sampname) %>%
        mutate(set=sampset) %>%
        mutate(cond=cond)
    
    allcover=rbind(allcover, cover)
}
        
sets=unique(allcover$set)
pdf(paste0(dbxdir, 'genome_cov.pdf'), h=6, w=16)
for (i in sets) {
    setcover=allcover %>%
        filter(set==i)
    plot=ggplot(setcover, aes(x=pos, y=cov, colour=sample)) +
        geom_line() +
        ggtitle(paste0('Sample ', i)) +
        xlab('Genome Position') +
        ylab('Coverage') +
        theme_bw()
    print(plot)
}
dev.off()



##pdf(paste0(dbxdir, 'genome_cov_all.pdf'), h=6, w=16)
pdf(paste0(dbxdir, 'genome_cov_all_trimmed.pdf'), h=6, w=16)
plot=ggplot(allcover, aes(x=pos, y=cov, colour=sample)) +
    geom_line() +
    ggtitle('All Samples') +
    xlab('Genome Position') +
    ylab('Coverage') +
    theme_bw()
print(plot)
plot=ggplot(allcover, aes(x=pos, y=cov, colour=sample)) +
    geom_line() +
    scale_y_continuous(trans='pseudo_log') +
    ggtitle('All Samples') +
    xlab('Genome Position') +
    ylab('Coverage') +
    theme_bw()
print(plot)
dev.off()


pdf(paste0(dbxdir, 'genome_cov_hist.pdf'), h=6, w=16)
plot=ggplot(allcover, aes(x=sample, y=cov, colour=sample, fill=sample, alpha=.3)) +
    geom_violin(scale='width') +
    ggtitle('All Samples') +
    xlab('Sample') +
    ylab('Coverage') +
    theme_bw()
print(plot)

##pseudocount to do log scale
plot=ggplot(allcover, aes(x=sample, y=cov, colour=sample, fill=sample, alpha=.3)) +
    geom_violin(scale='width') +
    scale_y_continuous(trans='pseudo_log') +
    ggtitle('All Samples') +
    xlab('Sample') +
    ylab('Coverage') +
    theme_bw()
print(plot)
dev.off()
    

##coverage by pool
bedfile='~/software/artic-ncov2019/primer_schemes/nCoV-2019/V3/nCoV-2019.bed'
colnames=c('chr', 'start', 'end', 'primer', 'pool', 'strand')
pools=read_tsv(bedfile, col_names=colnames) %>%
    rowwise() %>%
    mutate(amplicon=paste0(strsplit(primer, '_', fixed=TRUE)[[1]][1], '_' ,strsplit(primer, '_', fixed=TRUE)[[1]][2])) %>%
    group_by(amplicon) %>%
    summarise(start=min(start), end=max(end), pool=paste0('pool', strsplit(pool[1], '_', fixed=TRUE)[[1]][2]))

find_pool  <- function(position, pools) {
    pool= position>=pools$start & position<=pools$end
    if (sum(pool)==0) {
        status='none'
    } else if (sum(pool)==2) {
        status='primer'
    } else if (sum(pool)==1) {
        status=pools$pool[pool]
    } else {
        status=sum(pool)
    }
    return(status)
}    


poolcover=allcover %>%
    rowwise() %>%
    mutate(pool=find_pool(pos, pools)) %>%
    filter(pool!='none' && pool !='primer') %>%
    mutate(finalpool=paste0(cond, pool))


pdf(paste0(dbxdir, 'genome_cov_histpools.pdf'), h=10, w=25)
plot=ggplot(poolcover, aes(x=finalpool, y=cov, colour=pool, fill=pool, alpha=cond)) +
    geom_violin(scale='width') +
    scale_alpha_discrete(range=c(.2,.5)) +
    facet_wrap(~set, nrow=2) +
    ggtitle('Coverage Distribution') +
    xlab('Sample') +
    ylab('Coverage') +
    theme_bw()
print(plot)
dev.off()
