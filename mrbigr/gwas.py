import mrbigr.multiprocess as mu
from pandas_plink import read_plink1_bin
import pandas as pd
import numpy as np
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.rinterface_lib.embedded import RRuntimeError
import rpy2.robjects as robjects
from collections import Counter
import pyranges as pr
import warnings
import glob
import os
import shutil
import re

warnings.filterwarnings("ignore")


def read_genotype(geno_prefix):
    try:
        G = read_plink1_bin(geno_prefix + '.bed', geno_prefix + '.bim', geno_prefix + '.fam', ref='a0',verbose=False)
    except Exception:
        return None
    return G


def gwas_lmm(phe, geno, num_threads):
    geno_prefix = geno.split('/')[-1]
    related_matrix_cmd = 'gemma.linux -bfile {0}.link -gk 1 -o {1}'.format(geno_prefix,geno_prefix)
    gwas_cmd = 'gemma.linux -bfile {0}.link -k output/{0}.cXX.txt -lmm -n {1} -o {2}'
    fam = pd.read_csv(geno+'.fam', sep=r'\s+', header=None)
    fam[5] = 1
    fam = pd.merge(fam, phe, left_on=0, right_index=True, how='left')
    fam.to_csv(geno_prefix+'.link.fam', sep='\t', na_rep='NA', header=None, index=False)
    if os.path.exists(geno_prefix+'.link.bed'):
        os.remove(geno_prefix+'.link.bed')
    if os.path.exists(geno_prefix+'.link.bim'):
        os.remove(geno_prefix+'.link.bim')
    os.symlink(geno+'.bed', geno_prefix+'.link.bed')
    os.symlink(geno+'.bim', geno_prefix+'.link.bim')
    values = list()
    for _, p in enumerate(phe.columns):
        p = p.replace('/', '.')
        values.append((gwas_cmd.format(*[geno_prefix, _ + 2, str(p)]),))
    if not os.path.exists('output/'+geno_prefix+'.cXX.txt'):
        mu.run(related_matrix_cmd)
    s = mu.parallel(mu.run, values, num_threads)
    os.remove(geno_prefix+'.link.bed')
    os.remove(geno_prefix+'.link.bim')
    os.remove(geno_prefix+'.link.fam')
    return s

def gwas_lm(phe, geno, num_threads):
    geno_prefix = geno.split('/')[-1]
    gwas_cmd = 'gemma.linux -bfile {0}.link -lm -n {1} -o {2}'
    fam = pd.read_csv(geno+'.fam', sep=r'\s+', header=None)
    fam[5] = 1
    fam = pd.merge(fam, phe, left_on=0, right_index=True, how='left')
    fam.to_csv(geno_prefix+'.link.fam', sep='\t', na_rep='NA', header=None, index=False)
    if os.path.exists(geno_prefix+'.link.bed'):
        os.remove(geno_prefix+'.link.bed')
    if os.path.exists(geno_prefix+'.link.bim'):
        os.remove(geno_prefix+'.link.bim')
    os.symlink(geno+'.bed', geno_prefix+'.link.bed')
    os.symlink(geno+'.bim', geno_prefix+'.link.bim')
    values = list()
    for _, p in enumerate(phe.columns):
        p = p.replace('/', '.')
        values.append((gwas_cmd.format(*[geno_prefix, _ + 2, str(p)]),))
    s = mu.parallel(mu.run, values, num_threads)
    os.remove(geno_prefix+'.link.bed')
    os.remove(geno_prefix+'.link.bim')
    os.remove(geno_prefix+'.link.fam')
    return s


def generate_clump_input(dir):
    if os.path.exists('./clump_input'):
        shutil.rmtree('./clump_input')
    os.mkdir('./clump_input')
    for fn in glob.glob(dir.strip('/')+'/*.assoc.txt'):
        filename = fn.split('/')[-1]
        assoc = pd.read_csv(fn, sep='\t')
        assoc = assoc[['rs', 'p_wald']]
        assoc.columns = ['SNP', 'P']
        assoc.to_csv('./clump_input/' + filename.replace('.assoc.txt', '.assoc'), index=False, sep='\t')


def plink_clump(geno_path, p1, p2, num_threads):
    if os.path.exists('./clump_result'):
        shutil.rmtree('./clump_result')
    os.mkdir('./clump_result')
    cmd = 'plink --bfile {0} --clump {1}  --clump-p1 {2} --clump-p2 {3} --clump-kb {4} --clump-r2 0.2 --out {5} --clump-allow-overlap'
    cmds = list()
    ms = list()
    for fn in glob.glob('./clump_input/*'):
        phe_name = fn.split('/')[-1].replace('.assoc','')
        cmds.append((cmd.format(geno_path, fn, p1, p2, str(500), './clump_result/' + phe_name),))
        ms.append(phe_name)
    s = mu.parallel(mu.run, cmds, num_threads)
    if sum(s) != 0:
        print(','.join(list(np.array(ms)[s]))+' do not successfully generated clumped file.')
    return s


def merge_qtl_phe(qtl):
    qtl = qtl.sort_values(by=['CHR','qtl_start'])
    merged_phe_qtl = list()
    for index,row in qtl.iterrows():
        if not merged_phe_qtl:
            merged_phe_qtl.append(row)
        else:
            if row['CHR'] != merged_phe_qtl[-1]['CHR']:
                merged_phe_qtl.append(row)
            else:
                if row['qtl_start'] < merged_phe_qtl[-1]['qtl_end'] + 1:
                    if row['P'] < merged_phe_qtl[-1]['P']:
                        merged_phe_qtl[-1]['P'] = row['P']
                        merged_phe_qtl[-1]['SNP'] = row['SNP']
                    merged_phe_qtl[-1]['qtl_start'] = min(merged_phe_qtl[-1]['qtl_start'],row['qtl_start'])
                    merged_phe_qtl[-1]['qtl_end'] = max(merged_phe_qtl[-1]['qtl_end'], row['qtl_end'])
                    merged_phe_qtl[-1]['SP2_num'] += row['SP2_num']
                else:
                    merged_phe_qtl.append(row)
    merged_phe_qtl = pd.DataFrame(merged_phe_qtl)
    return merged_phe_qtl


def generate_qtl(clump_result_dir, cutoff):
    qtl_res = list()
    for fn in glob.glob(clump_result_dir.strip('/')+'/*clumped'):
        phe_name = fn.split('/')[-1].replace(".clumped","")
        clump_result = pd.read_csv(fn,sep='\s+')
        clump_result = clump_result.loc[clump_result.SP2!='NONE',:]
        qtl = clump_result[['CHR','BP','SNP','P','SP2']]
        qtl.loc[:,'SP2_num'] = qtl['SP2'].apply(lambda x: len(x.split(',')))
        qtl.loc[:,'log10P'] = -np.log10(qtl['P'])
        if (qtl['SP2_num'] >= cutoff).sum() > 0:
            qtl['qtl_start'] = qtl['SP2'].apply(lambda x:int(re.findall(r's_?(\d+)',x)[0]))
            qtl['qtl_end'] = qtl['SP2'].apply(lambda x:int(re.findall(r's_?(\d+)',x)[-1]))
            qtl['phe_name'] = phe_name
            mer_qtl_filter = merge_qtl_phe(qtl)
            mer_qtl_filter.loc[:,'qtl_length'] = mer_qtl_filter['qtl_end'] - mer_qtl_filter['qtl_start'] + 1
            qtl_res.append(mer_qtl_filter.loc[mer_qtl_filter['SP2_num']>=cutoff,['CHR','qtl_start','qtl_end','SNP','P','SP2_num','qtl_length','phe_name']])
    if not qtl_res:
        qtl_res = pd.DataFrame()
    else:
        qtl_res = pd.concat(qtl_res)
    return qtl_res


def gwas_plot(res, p, prefix, t):
    try:
        pandas2ri.activate()
        data_table = importr('data.table')
        base = importr('base')
        w = data_table.fread(res, data_table=base.getOption("datatable.fread.datatable", False))
        w_subset = w.loc[w.p_wald <= float(p), :]
        m = w_subset[['rs', 'chr', 'ps', 'p_wald']]
        q = w[['rs', 'chr', 'ps', 'p_wald']]
        m.columns = ['SNP', 'Chromosome', 'Position', prefix]
        q.columns = ['SNP', 'Chromosome', 'Position', prefix]
        thresholdi = robjects.FloatVector([1.0 / w.shape[0], 1e-6, 1e-5])
        lim = -np.log10(min(w_subset['p_wald'])) + 2
        base.sink('/dev/null')
        for path in os.environ.get('PATH').split(':'):
            if re.search(r'MRBIGR/utils',path):
                robjects.r('source("'+path+'/CMplot.r")')
        CMplot = robjects.r['CMplot']
        CMplot(m, plot_type='m', col=robjects.StrVector(["grey30", "grey60"]), ylim=robjects.FloatVector([2, lim]), threshold=thresholdi,
                cex=robjects.FloatVector([0.5, 0.5, 0.5]), signal_cex=robjects.FloatVector([0.5, 0.5, 0.5]),
                threshold_col=robjects.StrVector(['red', 'green', 'blue']), chr_den_col=robjects.rinterface.NULL, amplify=True,
                signal_pch = robjects.IntVector([19, 19, 19]), dpi=300,
                signal_col=robjects.StrVector(['red', 'green', 'blue']), multracks=False, LOG10=True, file=t)
        CMplot(q, plot_type='q', col='grey30', threshold=thresholdi[0],
               signal_cex=robjects.FloatVector([0.5, 0.5, 0.5]), signal_pch=robjects.IntVector([19, 19, 19]),
               conf_int_col='gray', signal_col='red', multracks=False, LOG10=True, file=t, dpi=300)
        base.sink()
    except RRuntimeError:
        return 0
    except ValueError:
        return 0
    else:
        return 1


def gwas_plot_parallel(p, threads, t, dir):
    values = list()
    for fn in glob.glob(dir.strip('/')+'/*.assoc.txt'):
        values.append((fn, p, fn.split('/')[-1].replace(".assoc.txt",""), t))
    s = mu.parallel(gwas_plot, values, threads)
    return s


def multi_trait_plot(gwas_dir, qtl, prefix, t):
    pandas2ri.activate()
    data_table = importr('data.table')
    base = importr('base')
    bk = pd.DataFrame()
    for fn in glob.glob(gwas_dir.strip('/')+'/*.assoc.txt'):
        d = pd.read_csv(fn, sep='\t')
        pn = fn.split('/')[-1].replace(".assoc.txt","")
        if bk.empty:
            bk = d.copy()
            bk.loc[bk.p_wald <= 1e-5, 'p_wald'] = 1e-5
            bk.index = bk['rs']
        d.index = d['rs']
        d = d.reindex(d.index.intersection(bk.index))
        bk = bk.reindex(d.index.intersection(bk.index))
        for index, row in qtl.loc[qtl.phe_name == pn, :].iterrows():
            peak_pos = int(re.findall(r's_?(\d+)',row['SNP'])[0])
            chrom = row['CHR']
            sig_tmp = pd.concat([bk.loc[(bk.chr.astype(str) == str(chrom)) & (bk.ps >= peak_pos-1000000) & (bk.ps <= peak_pos+1000000), 'p_wald'],
                                 d.loc[(bk.chr.astype(str) == str(chrom)) & (bk.ps >= peak_pos-1000000) & (bk.ps <= peak_pos+1000000), 'p_wald']], axis=1)
            sig_tmp.columns = ['bk', 'phe']
            bk.loc[(bk.chr.astype(str) == str(chrom)) & (bk.ps >= peak_pos-1000000) & (bk.ps <= peak_pos+1000000), 'p_wald'] = sig_tmp.apply(
                lambda x: x['bk'] if x['bk'] < x['phe'] else x['phe'], axis=1)
    bk = bk.loc[bk.p_wald <= 1e-2, :]
    bk.loc[bk.p_wald <= 1e-20, 'p_wald'] = 1e-20
    bk = bk[['rs', 'chr', 'ps', 'p_wald']]
    thresholdi = robjects.FloatVector([1.0 / d.shape[0], 1e-6, 1e-5])
    lim = -np.log10(min(bk['p_wald'])) + 2
    bk.columns = ['SNP', 'Chromosome', 'Position', prefix]
    base.sink('/dev/null')
    for path in os.environ.get('PATH').split(':'):
        if re.search(r'MODAS/utils', path):
            robjects.r('source("'+path+'/CMplot.r")')
    CMplot = robjects.r['CMplot']
    CMplot(bk, plot_type='m', col=robjects.StrVector(["grey30", "grey60"]), ylim=robjects.FloatVector([2, lim]),
           threshold=thresholdi,
           cex=robjects.FloatVector([0.5, 0.5, 0.5]), signal_cex=robjects.FloatVector([0.5, 0.5, 0.5]),
           threshold_col=robjects.StrVector(['red', 'green', 'blue']), chr_den_col=robjects.rinterface.NULL,
           amplify=True,
           signal_pch=robjects.IntVector([19, 19, 19]), dpi=300,
           signal_col=robjects.StrVector(['red', 'green', 'blue']), multracks=False, LOG10=True, file=t)
    base.sink()


def boxplot(phe, g, qtl,out):
    pandas2ri.activate()
    data_table = importr('data.table')
    base = importr('base')
    robjects.r('''box_plot <- function(d, phe, rs, level){
    library(ggplot2)
    library(ggsignif)
    d <- d[d$haplotype!=1,]
    d[d$haplotype==0,'haplotype'] <- level[1]
    d[d$haplotype==2,'haplotype'] <- level[2]
    d$haplotype <- factor(d$haplotype,levels = level)
    b <- as.numeric(formatC(max(d[,1],na.rm=T)*1.2/4,format = 'e',digits = 1))
    p <- ggplot(data = d,aes_string(x='haplotype',y=names(d)[1],fill='haplotype'))+
    theme_bw()+
    theme(legend.title = element_blank(),
          legend.background = element_blank(),
          legend.key = element_rect(colour = NA, fill = NA),
          legend.text = element_text(size = 4),
          legend.key.size = unit(3,'mm'),
          legend.position = 'none',
          panel.grid = element_blank(),
          axis.line = element_line(colour = 'black',size=0.4),
          axis.text = element_text(size = 6,color = 'black'),
          axis.ticks.length=unit(.1, 'cm'))+
    stat_boxplot(geom = 'errorbar', width = 0.2,size=0.1)+
    geom_boxplot(lwd=0.2,width=0.5,outlier.size = 0.2)+
    geom_signif(comparisons = list(level),map_signif_level = F,
                test= t.test, size=0.2 ,textsize=2, y_position = max(d[,1],na.rm=T)*1.1)+
    xlab('Group')+ylab('Value')+ggtitle(paste(phe,rs,spe=' '))+
    #scale_y_continuous(breaks=seq(0,4*b,by=b),labels = function(x) formatC(x, format = 'e',digits = 1), limits = c(0, max(d[,1], na.rm=T)*1.2))+
    scale_fill_manual(values=c('#E3FFE2', 'forest green'))
    ggsave(paste(phe,'_',rs,'_','boxplot','.pdf',sep=''),plot=p,width=3.5,height=4.5)
}''')
    robjects.r['options'](warn=-1)
    base.sink('/dev/null')
    box_plot = robjects.r['box_plot']
    g = g.where(g.snp.isin(qtl.SNP), drop=True)
    allele = pd.DataFrame([g.a0, g.a1], index=['a0', 'a1'], columns=g.snp.values)
    g = pd.DataFrame(g.values, index=g.sample, columns=g.snp.values)
    ril = g.index.intersection(phe.index)
    g = g.reindex(ril)
    phe = phe.reindex(ril)
    for index, row in qtl.iterrows():
        d = pd.concat([phe[row['phe_name']], g[row['SNP']]], axis=1)
        d.columns = ['trait.'+d.columns[0].replace('-', '.'), 'haplotype']
        genotype = str(allele[row['SNP']]['a1'].values)+str(allele[row['SNP']]['a0'].values)
        level = robjects.StrVector([allele[row['SNP']]['a1'].values*2, allele[row['SNP']]['a0'].values*2])
        box_plot(d, row['phe_name'], row['SNP'], level)
        d.to_csv(out+'/'+row['phe_name']+'_'+row['SNP']+'_'+genotype+'.csv', na_rep='NA')
    base.sink()


def qtl_anno(qtl, anno):
    anno = anno[['geneid', 'position']]
    anno.loc[:, 'chr'] = anno['position'].apply(lambda x: x.split(':')[0])
    anno.loc[:, 'start'] = anno['position'].apply(lambda x: x.split(':')[1].split('-')[0])
    anno.loc[:, 'end'] = anno['position'].apply(lambda x: x.split(':')[1].split('-')[1])
    anno = anno[['chr', 'start', 'end', 'geneid']]
    anno.columns = ['Chromosome', 'Start', 'End', 'geneid']
    qtl_range = qtl[['CHR', 'qtl_start', 'qtl_end', 'phe_name']]
    qtl_range.columns = ['Chromosome', 'Start', 'End', 'phe_name']
    qtl_range = pr.PyRanges(qtl_range)
    anno = pr.PyRanges(anno)
    qtl_anno_intersect = qtl_range.join(anno, how='left')
    qtl_anno = pd.DataFrame()
    for k in sorted(qtl_anno_intersect.dfs.keys()):
        qtl_anno = pd.concat([qtl_anno, qtl_anno_intersect.dfs[k]])
    qtl_anno = qtl_anno[['Chromosome', 'Start', 'End', 'phe_name', 'geneid']]
    qtl_anno.loc[:, 'Chromosome'] = qtl_anno.Chromosome.astype(str)
    qtl_anno.loc[:, 'Start'] = qtl_anno.Start.astype(int)
    qtl_anno.loc[:, 'End'] = qtl_anno.End.astype(int)
    qtl.loc[:, 'CHR'] = qtl.CHR.astype(str)
    qtl = pd.merge(qtl, qtl_anno, left_on=['CHR', 'qtl_start', 'qtl_end', 'phe_name'],
                   right_on=['Chromosome', 'Start', 'End', 'phe_name'])
    qtl = qtl.drop(['Chromosome', 'Start', 'End'], axis=1)
    qtl = qtl.groupby(['CHR', 'qtl_start', 'qtl_end', 'SNP', 'P', 'SP2_num', 'qtl_length', 'phe_name'])['geneid'].apply(';'.join).reset_index()
    qtl.columns = ['CHR', 'qtl_start', 'qtl_end', 'SNP', 'P', 'SP2_num', 'qtl_length', 'phe_name', 'qtl_all_gene']
    qtl.loc[qtl.qtl_all_gene == '-1', 'qtl_all_gene'] = 'nan'
    return qtl

def genetest(f1, f2, vcf, p, o):
    try:
        os.system('filtFromQtlAnno.pl {0} {1} >{2}'.format(f1, f2, o+'/matched.bed'))
        os.system('genoMatrix.pl {0} {1} |geneHaplotypeS1.pl - |geneHaplotypeS2.pl {0} - >{2}'.format(o+'/matched.bed', vcf, o+'/gene.haplotype'))
        os.system('tTestTrait.pl {0} {1} {2} {3}'.format(p, o+'/gene.haplotype', f1, o+'/gene_haplotest'))
    except Exception:
        return False
    else:
        return True
