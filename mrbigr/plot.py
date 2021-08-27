import pandas as pd
import numpy as np
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
import rpy2.robjects as robjects
from rpy2.rinterface_lib.embedded import RRuntimeError
import os, re


def forest_plot(MR_res, mTrait, snp, o):
    robjects.r('''
        p <- function(d, mTrait, snp, o)
        {
            suppressMessages(library(ggplot2))
            suppressMessages(library(scales))
            options(warn=-1)
            #d < - d[o,]
            d$pTrait <- factor(d$pTrait, levels = d$pTrait)
            d$std <- sqrt((d$effect) ^ 2 / d[, 5])
            d$min <- d$effect - d$std
            d$max <- d$effect + d$std

            # effect
            ggplot(data=d, aes(x=effect, y=pTrait, color=pTrait)) +
            geom_vline(xintercept=0, linetype="dashed", color="gray", size=1) +
            geom_point(size=1.5, shape=15) +
            geom_errorbarh(aes(xmin=min, xmax=max), height=.1) +
            theme_bw() +
            theme(legend.title = element_blank(),
                  legend.position = 'none',
                  plot.title = element_text(hjust=0.5),
                  panel.grid = element_blank(),
                  axis.line = element_line(colour='black', size=0.1),
                  axis.text.x = element_text(color='black'),
                  axis.text.y = element_text(color=hue_pal()(nrow(d))))+
                  xlab('') + ylab('') +
                  ggtitle(paste(mTrait, snp, sep='_'))
            ggsave(paste(o, mTrait, snp, 'forestplot.pdf', sep='_'), device='pdf', height=6, width=4)
        }
    ''')
    p = robjects.r('p')
    try:
        p(MR_res, mTrait, snp, o)
        return True
    except RRuntimeError:
        return False


def scatter_plot_mr(MR_res, g, mTrait, snp, o):
    robjects.r('''
        p <- function(d, a, mTrait, snp, o){
            suppressMessages(library(ggplot2))
            suppressMessages(library(ggrepel))
            if(nrow(a)==0){
                res <- d
                res$group <- 'group'
            }else{
                res <- merge(a,d,by.x ='id',by.y='pTrait')
                res$group <- factor(res$group)
                res <- res[order(res$group),]
            }
            res$log10p <- -log10(res$pvalue)
            res[res$effect<0,'log10p'] <- -res[res$effect<0,'log10p']
            #res$chrom <- factor(res$chrom, levels = paste('chr',sort(unique(res$chr)),sep=''))
            #res <- res[order(-res$log10p),]
            res$pos <- seq(1,nrow(res))
            scatter_plot <- ggplot(data=res,aes(x=pos,y=log10p,color=group))+
                  geom_hline(yintercept=log10(0.05/nrow(res)), linetype="dashed", color = "red", size=1)+
                  geom_hline(yintercept=-log10(0.05/nrow(res)), linetype="dashed", color = "red", size=1)+
                  geom_point(size=1.5,shape=19)+
                  theme_bw()+
                  theme(legend.title = element_blank(),
                        plot.title = element_text(hjust = 0.5),
                        panel.grid = element_blank(),
                        axis.line = element_line(colour = 'black',size=0.1),
                        axis.text = element_text(color = 'black'),
                        axis.text.x = element_blank(),
                        axis.ticks.length.x = unit(0,'cm'))+
                  xlab('')+ylab('')
            if(nrow(a)!=0){
                #scatter_plot <- scatter_plot + geom_text(aes(label=id),hjust=0, vjust=0)
                scatter_plot <- scatter_plot + geom_label_repel(aes(label=id))
            }
            ggsave(paste(o, mTrait, snp, 'scatterplot.pdf', sep='_'), plot=scatter_plot, device = 'pdf',height = 5,width = 9)
        }
    ''')
    p = robjects.r('p')
    try:
        p(MR_res, g, mTrait, snp, o)
        return True
    except RRuntimeError:
        return False


def manhattan_plot(d_fn, o):
    try:
        data_table = importr('data.table')
        base = importr('base')
        d = data_table.fread(d_fn, data_table=base.getOption("datatable.fread.datatable", False))
        #d.columns = ['SNP', 'Chromosome', 'Position', o]
        thresholdi = robjects.FloatVector([1.0 / d.shape[0], 1e-6, 1e-5])
        lim = -np.log10(d[d.columns[3:]].min().values) + 2
        base.sink('/dev/null')
        for path in os.environ.get('PATH').split(':'):
            if re.search(r'MRBIGR/utils', path):
                robjects.r('source("' + path + '/CMplot.r")')
        CMplot = robjects.r['CMplot']
        CMplot(d, plot_type='m', col=robjects.StrVector(["grey30", "grey60"]), ylim=robjects.FloatVector([2, lim]),
               threshold=thresholdi,
               cex=robjects.FloatVector([0.5, 0.5, 0.5]), signal_cex=robjects.FloatVector([0.5, 0.5, 0.5]),
               threshold_col=robjects.StrVector(['red', 'green', 'blue']), chr_den_col=robjects.rinterface.NULL,
               amplify=True,
               signal_pch=robjects.IntVector([19, 19, 19]), dpi=300,
               signal_col=robjects.StrVector(['red', 'green', 'blue']), multracks=False, LOG10=True, file='pdf', file_prefix=o)
        base.sink()
    except RRuntimeError:
        return False
    except ValueError:
        return False
    else:
        return True


def qq_plot(d_fn, o):
    try:
        data_table = importr('data.table')
        base = importr('base')
        d = data_table.fread(d_fn, data_table=base.getOption("datatable.fread.datatable", False))
        #d.columns = ['SNP', 'Chromosome', 'Position', o]
        base.sink('/dev/null')
        for path in os.environ.get('PATH').split(':'):
            if re.search(r'MRBIGR/utils', path):
                robjects.r('source("' + path + '/CMplot.r")')
        CMplot = robjects.r['CMplot']
        CMplot(d, plot_type='q', col='grey30', conf_int_col='gray', signal_col='red', multracks=False, LOG10=True,
               file='pdf', dpi=300, file_prefix=o)
        base.sink()
    except RRuntimeError:
        return False
    except ValueError:
        return False
    else:
        return True


def scatter_plot_ps(d, g, o):
    robjects.r('''
        p <- function(d, g, o){
            suppressMessages(library(ggplot2))
            d <- d[c(1,2,3)]
            if(nrow(g)!=0){
                d$group <- g[match(d[, 1], g[, 1]), 2]
            }else{
                d$group <- 'one'
            }
            names(d) <- c('sample','PC1','PC2','group')
            scatter_plot <- ggplot(data = d,aes(x=PC1,y=PC2,color=group))+
                geom_point(shape=19, size=1.5)+
                theme_bw() +
                theme(axis.text.x = element_text(colour = "black", size = 12, vjust =1 ),
                      axis.text.y = element_text(colour = "black", size = 12, hjust =1 ),
                      axis.title = element_text(margin=margin(10, 5, 0, 0), color = "black", size = 14),
                      axis.title.y = element_text(angle=90),
                      panel.border = element_rect(color = 'black', fill=NA, size = 1)
                     )+
                xlab('')+ylab('')
            if(nrow(g)==0){
                scatter_plot <- scatter_plot + theme(legend.position='none')
            }
            ggsave(paste(o,'scatter_ps.pdf', sep='_'), plot=scatter_plot, width = 4.2, height = 3)
        }
    ''')
    p = robjects.r('p')
    try:
        p(d, g, o)
        return True
    except RRuntimeError:
        return False


def hist_plot(d, o):
    robjects.r('''
        p <- function(d, o){
            suppressMessages(library(ggplot2))
            d <- na.omit(d)
            ggplot(data=d, aes_string(x=names(d)[2]))+
                geom_histogram(bins = 10, fill='gray80', color='black', alpha=0.9, size=0.5)+
                theme_bw() +
                theme(axis.text.x = element_text(colour = "black", size = 12, vjust =1 ),
                      axis.text.y = element_text(colour = "black", size = 12, hjust =1 ),
                      axis.title = element_text(margin=margin(10, 5, 0, 0), color = "black", size = 14),
                      axis.title.y = element_text(angle=90),
                      panel.border = element_rect(color = 'black', fill=NA, size = 1),
                      panel.grid = element_blank()
                )+
                ylab('Freq')
            ggsave(paste(o, 'hist.pdf', sep='_'),width = 3, height = 3)
        }
    ''')
    p = robjects.r('p')
    try:
        p(d, o)
        return True
    except RRuntimeError:
        return False


def box_plot(d, g, o):
    robjects.r('''
        p <- function(d, g, o){
            suppressMessages(library(ggplot2))
            if(nrow(g)!=0){
                d$group <- g[match(d[,1], g[,1]),2]
            }else{
                d$group <- names(d)[2]
            }
            group_num <- length(unique(d$group))
            if(group_num==1){
                width <- 2
            }else if(group_num==2){
                width <- 2.5
            }else{
                width <- 2.5 + 0.5 * (group_num-2)
            }
            if(width >= 8){
                width <- 8
            }
            d <- na.omit(d)
            ggplot(data=d,aes_string(y=names(d)[2],x='group', fill='group'))+
              stat_boxplot(geom = 'errorbar', width=0.2)+
              geom_boxplot()+
              theme_bw() +
              theme(axis.text.x = element_text(colour = "black", size = 12, vjust =1 ),
                    axis.text.y = element_text(colour = "black", size = 12, hjust =1 ),
                    axis.title = element_text(margin=margin(10, 5, 0, 0), color = "black", size = 14),
                    axis.title.y = element_text(angle=90),
                    panel.border = element_rect(color = 'black', fill=NA, size = 1),
                    legend.position = 'none'
              )+
              xlab('')
            ggsave(paste(o, 'boxplot.pdf', sep='_'), width = width, height = 3)
        }
    ''')
    p = robjects.r('p')
    try:
        p(d, g, o)
        return True
    except RRuntimeError:
        return False


def snp_density_plot(d_fn, o):
    robjects.r('''
        bin_range <- function(d){
            d <- d[order(as.numeric(d[, 2]), as.numeric(d[, 3])), ]
            chr <- as.numeric(d[, 2])
            pos <- as.numeric(d[, 3])
            chr.num <- unique(chr)
            chorm.maxlen <- max(pos)
            bin <- 1e6
            pos.x <- list()
            col.index <- list()
            maxbin.num <- NULL
            for(i in 1 : length(chr.num)){
              pos.x[[i]] <- pos[which(chr == chr.num[i])]
              cut.len <- ceiling((max(pos.x[[i]]) - min(pos.x[[i]])) / bin)
              if(cut.len <= 1){
                maxbin.num <- c(maxbin.num,length(pos.x[[i]]))
                col.index[[i]] <- rep(length(pos.x[[i]]), length(pos.x[[i]]))
              }else{
                cut.r <- cut(pos.x[[i]], cut.len, labels=FALSE)
                eachbin.num <- table(cut.r)
                maxbin.num <- c(maxbin.num, max(eachbin.num))
                col.index[[i]] <- rep(eachbin.num, eachbin.num)
              }
            }
            return(quantile(unlist(col.index), 0.95))
        }
    ''')
    bin_range = robjects.r('bin_range')
    try:
        data_table = importr('data.table')
        base = importr('base')
        d = data_table.fread(d_fn, data_table=base.getOption("datatable.fread.datatable", False))
        bin_max = bin_range(d)
        base.sink('/dev/null')
        for path in os.environ.get('PATH').split(':'):
            if re.search(r'MRBIGR/utils', path):
                robjects.r('source("' + path + '/CMplot.r")')
        CMplot = robjects.r['CMplot']
        CMplot(d, plot_type='d', bin_size=1e6, bin_max=bin_max[0], col=robjects.StrVector(['darkgreen', 'yellow', 'red']),
               file='pdf', memo='', dpi=300, file_output=True, verbose=True, file_prefix=o)
        base.sink()
    except RRuntimeError:
        return False
    except ValueError:
        return False
    else:
        return True
