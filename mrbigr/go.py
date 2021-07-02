import pandas as pd
import rpy2.robjects.packages as packages
from rpy2.robjects import pandas2ri
import rpy2.robjects as robjects
from rpy2.rinterface_lib.embedded import RRuntimeError
import warnings
import shutil
import os


warnings.filterwarnings("ignore")
pandas2ri.activate()
base = packages.importr('base')


def build_orgdb(go_anno, gene_info):
    if os.path.exists('./org.Custom.eg.db'):
        shutil.rmtree('./org.Custom.eg.db')
    go_anno.columns = ['GID', 'GO']
    go_anno.loc[:, 'EVIDENCE'] = 'IEA'
    gene_chr = gene_info[['gene_id', 'chr']]
    gene_chr.columns = ['GID', 'CHROMOSOME']
    gene_anno = gene_info[['gene_id', 'gene_id', 'annotation']]
    gene_anno.columns = ['GID', 'SYMBOL', 'GENENAME']
    try:
        robjects.r('''
            build_db <- function(gene_anno, gene_chr, go_anno){
                options(warn = - 1) 
                suppressMessages(library(AnnotationForge))
                sink('/dev/null')
                suppressMessages(makeOrgPackage(gene_info=gene_anno, chromosome=gene_chr, go=go_anno,
                     version='0.1',
                     maintainer='Custom <custom@someplace.org>',
                     author='Custom <custom@someplace.org>',
                     outputDir = '.',
                     tax_id='0000',
                     genus='C',
                     species='ustom',
                     goTable='go',
                     verbose=F))
                sink()
                install.packages('./org.Custom.eg.db',repos=NULL,type="source",unlink=TRUE,quiet = T)
            }
        ''')
        build_db = robjects.r('build_db')
        build_db(gene_anno, gene_chr, go_anno)
        shutil.rmtree('./org.Custom.eg.db')
        return True
    except RRuntimeError:
        return False


def cheak_orgdb(orgdb):
    try:
        packages.importr(orgdb)
        return True
    except RRuntimeError:
        return False


def go_enrich(gene, orgdb, pvalue):
    robjects.r('''
    enrich <- function(gene, orgdb){
        options(warn = - 1) 
        suppressMessages(library(orgdb, character.only = TRUE))
        suppressMessages(library(clusterProfiler))
        res <- data.frame()
        for(i in 1:nrow(gene)){
            tmp <- enrichGO(unlist(strsplit(gene[i,'genes'],' ')), keyType = 'GID', OrgDb = orgdb , ont='ALL', pvalueCutoff = 1, qvalueCutoff = 1)
            tmp_df <- tmp@result
            for(go_type in c('BP','MF','CC')){
                tmp_sub <- tmp_df[tmp_df$ONTOLOGY==go_type,]
                if((is.null(tmp_sub) || nrow(tmp_sub) == 0 || ncol(tmp_sub) == 0)){
                  next
                }
                tmp@result <- tmp_sub
                tmp@ontology <- go_type
                tmp_sim <- simplify(tmp, cutoff=0.7)
                tmp_sim_df <- tmp_sim@result
                #tmp_sim_df <- tmp@result
                tmp_sim_df <- tmp_sim_df[tmp_sim_df$Count>=3,]
                if(nrow(tmp_sim_df)==0){
                  next
                }
                tmp_sim_df$label <- i
                res <- rbind(res, tmp_sim_df)
            }
        }
        return(res)
    }
    ''')
    enrich = robjects.r['enrich']
    try:
        res = enrich(gene, orgdb)
        res = res.loc[res['p.adjust'] <= pvalue, :]
        return res
    except RRuntimeError:
        return None


def go_plot(res, plot_type, o):
    robjects.r('''
        p_go <- function(res, plot_type, o){
            options(warn = - 1)
            suppressMessages(library(clusterProfiler))
            suppressMessages(library(ggplot2))
            for(go_type in c('BP','MF','CC')){
                res_sub <- res[res$ONTOLOGY==go_type,]
                if((is.null(res_sub) || nrow(res_sub) == 0 || ncol(res_sub) == 0)){
                    next
                }
                x <- new("enrichResult",
                 result         = res_sub,
                 pvalueCutoff   = 0.05,
                 pAdjustMethod  = "BH",
                 qvalueCutoff   = 0.05,
                 gene           = unique(unlist(lapply(res_sub$geneID,function(x)unlist(strsplit(x,'/'))))),
                 universe       = 'Unknown',
                 geneSets       = list(),
                 organism       = 'Unknown',
                 keytype        = 'Unknown',
                 ontology       = go_type,
                 readable       = FALSE
                )
                for(p_type in plot_type){
                    if(p_type=='barplot')
                        p <- barplot(x)
                    else if(p_type=='dotplot')
                        p <- dotplot(x)
                    else if(p_type=='cnetplot')
                        p <- cnetplot(x)
                    else if(p_type=='heatplot')
                        p <- heatplot(x)
                    else if(p_type=='emapplot')
                        p <- emapplot(x)
                    else if(p_type=='upsetplot')
                        p <- upsetplot(x)
                    suppressMessages(ggsave(paste(o,go_type,p_type,'pdf',sep='.'), plot=p, device='pdf', width=9))
                }
            }
        }
    ''')
    p_go = robjects.r('p_go')
    try:
        p_go(res, plot_type, o)
        return True
    except Exception:
        return False
