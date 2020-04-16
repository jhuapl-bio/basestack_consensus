library(ggplot2)
library(tidyverse)
library(gridExtra)
library(ShortRead)
library(ggpubr)

##calldir='~/data/ncov/initial-clinical/v3_run1/20200405_0422_GA30000_FAN30842_00fb1614/guppy_demux/'
dbxdir='~/Dropbox/timplab_data/ncov/barcoding/'

runs=c('20200405_0422_GA30000_FAN30842_00fb1614', '20200410_2018_X4_FAN32204_327837a0')
modes=c('guppy_demux', 'guppy_demux_both')
for (run in runs) {
    for (mode in modes) {
        calldir=paste0('~/data/ncov/initial-clinical/', run, '/', mode)

        
        bcdata=read_csv(paste0(calldir, '/counts.csv'), col_names=c('bc', 'num')) %>%
            mutate(prop=format(round((num/sum(num))*100, 2),nsmall=2)) %>%
            mutate(propnum=(num/sum(num))*100)
        barcodes=bcdata$bc
        
    
        plots=vector(mode='list')
        for (i in 1:length(barcodes)) {
            fqdir=paste0(calldir, '/', barcodes[i])
            fqfile=paste0(fqdir,'/',list.files(fqdir))
            reads=readFastq(fqfile)
            readlen=width(sread(reads))
            rlens=tibble(samp=barcodes[i], readlen=readlen, propnum=bcdata$propnum[bcdata$bc==barcodes[i]])
            
            if (rlens$propnum[1]<4) {
                col='#F8766D'
            }else{
                col='#00BFC4'
            }
            
            plot=ggplot(rlens, aes(x=readlen)) +
                geom_density(alpha=.5, colour=col, fill=col) +
                ggtitle(paste0('guppy ', barcodes[i], ' ', bcdata$prop[bcdata$bc==barcodes[i]], '%')) +
                xlab('length (bp)') +
                xlim(0,1500) +
                theme_bw() +
                theme(plot.title=element_text(face='bold', colour=col)) +
                theme(legend.position = "none")
            
        
            plots[[i]]=plot
        }
        
        plotfile=paste0(dbxdir, run, '_', mode, '.pdf')
        pdf(plotfile, width=25, height=10)
        print(ggarrange(plotlist=plots, ncol=5, nrow=3, align="hv"))
        dev.off()
    }
}
