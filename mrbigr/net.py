import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.stats import expon
import rpy2.robjects.packages as packages
from rpy2.robjects import pandas2ri
import rpy2.robjects as robjects
import subprocess
import sys, os
import re


pandas2ri.activate()


def get_weight(pvalue_list):
    pvalue_list = pvalue_list[['mTrait', 'pTrait', 'pvalue']]
    pvalue_matrix = pd.pivot_table(pvalue_list, index='mTrait', columns='pTrait', values='pvalue', fill_value=1)
    pvalue_matrix = -np.log10(pvalue_matrix) - 3
    pvalue_matrix[pvalue_matrix < 0] = 0
    weight = 1 - expon.pdf(pvalue_matrix)
    weight = (weight + weight.T) / 2
    weight_csr = csr_matrix(np.triu(weight, k=1))
    weight_pair = list()
    gene_id = pvalue_matrix.columns
    for row in range(weight.shape[0]):
        for col, w in zip(weight_csr.indices[weight_csr.indptr[row]:weight_csr.indptr[row + 1]], weight_csr.data[weight_csr.indptr[row]:weight_csr.indptr[row + 1]]):
            weight_pair.append([gene_id[row], gene_id[col], w])
    weight_pair = pd.DataFrame(weight_pair)
    return weight_pair


def module_identify(edge_weight_fn, module_size):
    prefix = edge_weight_fn.replace('.edge_list', '')
    cl1_path = ''
    for path in os.environ.get('PATH').split(':'):
        if re.search(r'MODAS/utils', path):
            cl1_path = path.rstrip('/')
    if not cl1_path:
        return None
    #print('java -jar ' + cl1_path + '/cluster_one-1.0.jar -s ' + str(module_size) + ' -f edge_list -F csv ' + edge_weight_fn + ' >' + prefix + '.cluster_one.result.csv 2>/dev/null')
    subprocess.call(
        'java -jar ' + cl1_path + '/cluster_one-1.0.jar -s ' + str(module_size) + ' -f edge_list -F csv ' + edge_weight_fn + ' >' + prefix + '.cluster_one.result.csv 2>/dev/null',
        shell=True)
    cluster_one_res = pd.read_csv(prefix + '.cluster_one.result.csv')
    cluster_one_res = cluster_one_res.loc[cluster_one_res['P-value'] <= 0.05, :]
    cluster_one_res = cluster_one_res[['Size', 'Members']].sort_values(by='Size', ascending=False)
    cluster_one_res.loc[:, 'module'] = np.arange(1, cluster_one_res.shape[0] + 1)
    cluster_one_res.columns = ['gene_num', 'genes', 'module']
    cluster_one_res = cluster_one_res[['module', 'gene_num', 'genes']]
    return cluster_one_res


def hub_identify(edge_weight, cluster_one_res):
    robjects.r('''
        hub_ide <- function(m_ew){
            g <- igraph::graph.data.frame(m_ew,directed = FALSE)
            hub <- igraph::hub_score(g)
            hub_df <- as.data.frame(hub$vector,stringsAsFactors=F)
            names(hub_df) <- 'hub_score'
            hub_df$gene_id <- rownames(hub_df)
            rownames(hub_df) <- NULL
            return(hub_df)
        }
    ''')
    hub_ide = robjects.r('hub_ide')
    hub_list = list()
    hub_res = pd.DataFrame()
    for index, row in cluster_one_res.iterrows():
        gene_list = row['genes'].split(' ')
        m_ew = edge_weight.loc[(edge_weight['row'].isin(gene_list)) & (edge_weight['col'].isin(gene_list)), :]
        hub = hub_ide(m_ew)
        hub_res = pd.concat([hub_res, hub])
        hub_list.append(' '.join(hub.loc[hub.hub_score >= 0.8, :].sort_values(by='hub_score', ascending=False).apply(lambda x: x['gene_id']+'('+str(x['hub_score'])+')', axis=1).values))
    cluster_one_res.loc[:, 'hub_gene'] = hub_list
    return cluster_one_res, hub_res


def module_network_plot(edge_weight, cluster_one_res, hub_res, prefix):
    robjects.r('''
        network_plot <- function(m_ew, m_hub, fn){
            options(warn = - 1) 
            suppressMessages(library(network))
            suppressMessages(library(ggplot2))
            net <- network(m_ew, matrix.type = 'edgelist', directed = F)
            m_hub <- m_hub[match(net%v%'vertex.names', m_hub$gene_id),]
            m_hub$hub_score_range <- '0-0.3'
            m_hub[m_hub$hub_score>=0.3 & m_hub$hub_score<0.8, 'hub_score_range'] <- '0.3-0.8'
            m_hub[m_hub$hub_score>=0.8, 'hub_score_range'] <- '0.8-1'
            m_hub$hub_score_range <- factor(m_hub$hub_score_range)
            m_hub$size <- 0.6
            m_hub[m_hub$hub_score>=0.8,'size'] <- 3
            col <- rev(scales::hue_pal()(3))[3-length(unique(m_hub$hub_score_range))+1:3]
            names(col) <- levels(m_hub$hub_score_range)
            suppressMessages(GGally::ggnet2(net, edge.size='weight',  mode = 'kamadakawai', edge.color = 'gray80',
                size = m_hub$size, label = m_hub[m_hub$hub_score>=0.8, 'gene_id'],
                color.legend = 'hub_score', color = m_hub$hub_score_range, palette=col) +
                theme(legend.position = 'bottom') +
                scale_size_discrete(guide = 'none'))
            suppressMessages(ggsave(fn, device='pdf'))
        }
    ''')
    network_plot = robjects.r('network_plot')
    for index, row in cluster_one_res.iterrows():
        gene_list = row['genes'].split(' ')
        m_ew = edge_weight.loc[(edge_weight['row'].isin(gene_list)) & (edge_weight['col'].isin(gene_list)), :]
        m_hub = hub_res.loc[hub_res['gene_id'].isin(gene_list), :]
        network_plot(m_ew, m_hub, prefix+'_module'+str(row['module'])+'.pdf')

