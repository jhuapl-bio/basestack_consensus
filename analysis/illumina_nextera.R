library(tidyverse)

dbxdir='~/Dropbox/timplab_data/ncov/nextera/'

datadir='/uru/Data/Nanopore/projects/ncov/nextera/'
covdir=paste0(datadir, 'cov/')

##covfiles=list.files(covdir)
covfiles=list.files(covdir, pattern="\\.trimmed.cov$")

allcover=tibble(
    chr=as.character(),
    pos=as.numeric(),
    cov=as.numeric(),
    sample=as.character())

##pdf(paste0(dbxdir, 'genome_cov.pdf'), h=6, w=16)
pdf(paste0(dbxdir, 'genome_cov_trimmed.pdf'), h=6, w=16)
for (i in covfiles) {
    covfile=paste0(covdir, i)
    cover=read_tsv(covfile, col_names=c('chr', 'pos', 'cov')) %>%
        mutate(sample=i)

    allcover=rbind(allcover, cover)
    
    plot=ggplot(cover, aes(x=pos, y=cov)) +
        geom_line() +
        ggtitle(i) +
        xlab('Genome Position') +
        ylab('Coverage') +
        theme_bw()
    plot=ggplot(cover, aes(x=pos, y=cov)) +
        geom_line() +
        scale_y_continuous(trans='pseudo_log') +
        ggtitle(i) +
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


##pdf(paste0(dbxdir, 'genome_cov_hist.pdf'), h=6, w=16)
pdf(paste0(dbxdir, 'genome_cov_hist_trimmed.pdf'), h=6, w=16)
plot=ggplot(allcover, aes(x=sample, y=cov, colour=sample, fill=sample, alpha=.3)) +
    geom_violin(scale='width') +
    ggtitle('All Samples') +
    xlab('Sample') +
    ylab('Coverage') +
    theme_bw()
print(plot)

##pseudocount to do log scale
allcover=allcover %>%
    mutate(cov1=cov+1)
plot=ggplot(allcover, aes(x=sample, y=cov, colour=sample, fill=sample, alpha=.3)) +
    geom_violin(scale='width') +
    scale_y_continuous(trans='pseudo_log') +
    ggtitle('All Samples') +
    xlab('Sample') +
    ylab('Coverage') +
    theme_bw()
print(plot)
dev.off()
    
