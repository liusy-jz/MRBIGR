import mrbigr.geno as mg
import mrbigr.pheno as mp
import mrbigr.gwas as mw
import mrbigr.multiprocess as mu
import mrbigr.mr as mmr
import mrbigr.go as mgo
import mrbigr.net as mnet
import mrbigr.plot as mplt
import argparse
import sys, time, os, glob
import traceback
import pandas as pd
import numpy as np
import shutil, re


class Logger(object):

	def __init__(self, fn):
		self.log_file = open(fn, 'w')

	def log(self, msg):
		'''
		Print to log file and stdout with a single command.

		'''
		if isinstance(msg,list):
			msg = ' '.join(msg)
		else:
			msg = ('{T}: '+msg).format(T=time.ctime())
		print(msg,file=self.log_file)
		print(msg)

class ArgsError(Exception):
	def myexcepthook(type, value, tb):
		msg = ''.join(traceback.format_exception(type, value, tb))
		print(msg, end = '')

	sys.excepthook = myexcepthook


class FileTypeError(Exception):
	def myexcepthook(type, value, tb):
		msg = ''.join(traceback.format_exception(type, value, tb))
		print(msg, end = '')

	sys.excepthook = myexcepthook

def geno(args,log):

	if args.convert:
		if args.hapmap:
			log.log('Begin convert genotype file to plink bed format.')
			b = mg.hap2plink_bed(args.hapmap, args.o)
			if b:
				log.log('Successfully convert genotype file to plink bed format.')
			else:
				raise FileTypeError('Genotype format conversion raised an error!')
		if args.vcf:
			b = mg.vcf2plink_bed(args.vcf, args.o)
		if args.g:
			b = mg.plink_bed2vcf(args.g, args.o)

	elif args.anno:
		log.log('Begin to genotype data annotation using ANNOVAR...')
		if not (args.db and args.dbdir and args.vcf and args.o):
			raise ArgsError('Error: the parameters -db, -dbdir, -vcf and -o must be set!')
		if not os.path.exists(args.dbdir):
			os.makedirs(args.dbdir)
		if args.gtf and args.fa:
			try:
				os.system('annovardb.sh {0} {1} {2}/{3}'.format(args.gtf, args.fa, args.dbdir, args.db))
			except Exception:
				log.log('Error when building annotation database!')
			else:
				log.log('Successfully building the annotation database.')
		b = mg.geno_anno(args.vcf, args.dbdir, args.db, args.o)
		if not b:
			log.log('Error when genotype data annotation!')
		else:
			log.log('Successfully annotated the genotype data.')

	else:
		if not os.path.exists(args.g) and not os.path.exists(args.g+'.bed'):
			raise ArgsError('Genotype file {} is not exists.'.format(args.g))

		if args.qc:
			maf = 0.05
			mis = 0.1
			mind = 0.99
			if args.maf:
				maf = args.maf
			if args.mis:
				mis = args.mis
			if args.mind:
				mind = args.mind
			b = mg.snp_qc(args.g, args.o + '_qc', maf, mis, mind)

		if args.imput:
			method = 'mode'
			if args.method:
				method = args.method
			log.log('Begin genotype imputation...')
			if os.path.exists(args.o + '_imput.bed'):
				os.remove(args.o + '_imput.bed')
				os.remove(args.o + '_imput.bim')
				os.remove(args.o + '_imput.fam')
			b = mg.snp_imput(args.g+'.bed', args.o + '_imput.bed', method)
			if b:
				log.log('Genotype imputation is done.')
			else:
				log.log('Genotype imputation raised an error!')

		if args.clump:
			log.log('Begin snp clumping...')
			if os.path.exists(args.o + '_clump.bed'):
				os.remove(args.o + '_clump.bed')
				os.remove(args.o + '_clump.bim')
				os.remove(args.o + '_clump.fam')
			r2 = 0.8
			if args.r2:
				r2 = args.r2
			g_clump_list = mg.snp_clumping(args.g+'.bed', args.o + '_clump.bed', r2)
			log.log('SNP clumping is done.')
			log.log('There are {0} snps in input genotype file, after clumping, {1} snp left.'.format(g_clump_list[0],g_clump_list[1]))
			log.log('Clumped genotype file is saved as {0}_clump.bed, {0}_clump.bim and {0}_clump.fam.'.format(args.o))

		if args.kinship:
			log.log('Begin kinship/relatedness matrix generation...')
			ks = mg.kinship(args.g, args.o)
			if ks:
				log.log('Kinship analysis is done.')
			else:
				log.log('Kinship analysis raised an error!')

		if args.pca:
			log.log('Begin Principal Component Analysis...')
			dim = 10
			if args.dim:
				dim = int(args.dim)
			pc = mg.pca(args.g+'.bed', dim)
			pc.to_csv(args.o +'_pca.csv', na_rep='NA')
			log.log('PCA Analysis is done.')

		if args.tsne:
			log.log('Begin t-SNE Analysis...')
			dim = 3
			if args.dim:
				dim = int(args.dim)
			tsne = mg.tsne(args.g+'.bed', dim)
			tsne.to_csv(args.o +'_tsne.csv', na_rep='NA')
			log.log('t-SNE Analysis is done.')

		if args.tree:
			log.log('Begin to build a genotype ML-tree using FastTree...')
			ft = mg.fasttree(args.g, args.o)
			if os.path.exists(args.o + '.fasta'):
				os.remove(args.o + '.fasta')
			if os.path.exists(args.o + '.ped'):
				os.remove(args.o + '.ped')
			if ft:
				log.log('A nwk format ML-tree is built successfully.')
			else:
				log.log('Tree built raised an error!')



def pheno(args,log):

	if not os.path.exists(args.p):
		raise ArgsError('The phenotype file {} is not exists.'.format(args.p))
	phe = pd.read_csv(args.p, index_col=0)
	if phe.empty:
		log.log('The phenotype file {} is empty!'.format(args.phe))
		sys.exit()
	else:
		log.log('Read phenotype file finished, there are {} phenotypes in phenotype file.'.format(phe.shape[1]))

	if args.qc:
		mis = 0.5
		val = 0.1
		rout = 'zscore'
		if args.tr:
			log.log('Transpose the input data')
			phe = mp.transpose(phe)
			log.log('Transposition finished.')
		if args.rout:
			log.log('Perform outlier removal using {} method'.format(rout))
			phe = mp.outlier(phe,rout)
			log.log('Outlier data have been replaced by NA.')
		if args.mis:
			mis = args.mis
			if args.mis < 0 or args.mis > 1:
				raise ArgsError('The interval of missing rate is between 0 and 1.')
			log.log('Filtering phenotype missing ratio larger than {}'.format(mis))
			phe = mp.missing_filter(phe, mis)
			log.log('After filtering, there are {} traits left.'.format(phe.shape[1]))
		if args.val:
			val = args.val
			if args.val <= 0:
				raise ArgsError('Average value of phenotype should be larger than zero.')
			log.log('Filtering phenotype values smaller than {} in average'.format(val))
			phe = mp.abundance_filter(phe, val)
			log.log('After filtering, there are {} phenotypes left.'.format(phe.shape[1]))

	if args.imput:
		method = "mean"
		if args.method:
			method = args.method
			log.log('Imputation of missing values in phenotype data.')
			phe = mp.pheno_imputer(phe, method)
		log.log('Successfully imputed phenotype data.')

	if args.scale:
		if args.log2:
			log.log('log2 scale of the phenotype data')
			phe = mp.log2_scale(phe)
		if args.log10:
			log.log('log10 of the phenotype data')
			phe = mp.log10_scale(phe)
		if args.ln:
			log.log('ln of the phenotype data')
			phe = mp.ln_scale(phe)
		if args.boxcox:
			log.log('Box-Cox normalization of the phenotype data')
			phe = mp.boxcox_scale(phe)
		if args.qqnorm:
			log.log('Normal quantile transform of the phenotype data')
			phe = phe.apply(mp.qqnorm)
		if args.minmax:
			log.log('Min-Max normalization of the phenotype data')
			phe = mp.minmax_scale(phe)
		if args.zscore:
			log.log('Z-score normalization of the phenotype data')
			phe = mp.zscore_scale(phe)
		if args.robust:
			log.log('Scale phenotype data using statistics that are robust to outliers')
			phe = mp.robust_scale(phe)
		log.log('Successfully scaled phenotype data.')

	if args.correct:
		log.log('Correct phenotype values using pre-calculated PCA file.')
		if not os.path.exists(args.pca):
			raise ArgsError('The PCA file {} is not exists.'.format(args.pca))
		pca = pd.read_csv(args.pca, index_col=0)
		pca_phe_idx = pca.index.intersection(phe.index)
		phe = mp.trait_correct(pca.reindex(pca_phe_idx), phe.reindex(pca_phe_idx))
		log.log('Successfully corrected phenotype data using PCA file.')

	if args.merge:
		mm = 'mean'
		if args.mm:
			mm = args.mm
		log.log('Merge phenotype values from different enviroment using {} method.'.format(mm))
		if mm == 'mean':
			phe = mp.mean(phe)
		if mm == 'blup':
			phe = mp.blup(phe)
		if mm == 'blue':
			phe = mp.blue(phe)
		log.log('Successfully merged phenotype values from different enviroment.')

	phe.to_csv(args.o+'.phe.csv', na_rep='NA')
	log.log('phenotype processing is finished.')


def gwas(args, log):

	if args.gwas:
		log.log('Begin genome-wide association analysis and visualization')
		if not os.path.exists(args.p):
			raise ArgsError('The phenotype file {} is not exist.'.format(args.p))
		phe = pd.read_csv(args.p, index_col=0)
		if phe.empty:
			log.log('The phenotype file {} is empty and software is terminated.'.format(args.p))
			sys.exit()
		if not os.path.exists(args.g+'.bed'):
			raise ArgsError('the plink bed format genotype file {} is not exists.'.format(args.g))
		g = mw.read_genotype(args.g)
		if g is None:
			raise FileTypeError('{0}.bed is not a plink bed format file or {0}.bim file and {0}.fam file is not in the same folder with {0}.bed, or {0}.bed is empty.'.format(args.g))

		log.log('Begin GWAS analysis')
		thread = int(1)
		if args.thread:
			thread = args.thread
		model = 'lmm'
		if args.model:
			model = args.model
		if model == 'lmm':
			s = mw.gwas_lmm(phe, args.g, thread)
		if model == 'lm':
			s = mw.gwas_lm(phe, args.g, thread)
		if s is None:
			log.log('GWAS faild, please cheak your data.')
			sys.exit()
		else:
			status = np.array(s) == 0
			fail_phe = phe.columns[~status].astype(str)
			if not status.all():
				log.log('There are {0} GWAS is failed, they are {1}'.format(np.sum(~status), ','.join(fail_phe)))
			else:
				log.log('All phenotype GWAS done.')

	if args.qtl:
		log.log('Begin QTL analysis...')
		input_dir = 'output' # the default gwas output directory
		if args.i:
			input_dir = args.i
		if not os.path.exists(input_dir):
			raise ArgsError('The gwas output *.assoc file directory {} is not exist.'.format(input_dir))
		p1 = 1e-7
		if args.p1:
			p1 = args.p1
		p2 = 1e-5
		if args.p2:
			p2 = args.p2
		p2n = 5
		if args.p2n:
			p2n = args.p2n
		thread = int(1)
		if args.thread:
			thread = args.thread
		s = mw.generate_clump_input(input_dir)
		s = mw.plink_clump(args.g, p1, p2, thread)
		qtl_res = mw.generate_qtl('clump_result', p2n)
		qtl_res.index = np.arange(qtl_res.shape[0])
		qtl_res.to_csv(args.o+'.qtl_res.csv', index=False)
		log.log('Generate QTL result done.')

	if args.anno:
		log.log('Begin QTL annotation...')
		if not args.a and args.q:
			raise ArgsError('The parameters -q and -a must be set correctly!')
		if not os.path.exists(args.a):
			raise ArgsError('The gene annotation file {} is not exists.'.format(args.a))
		anno = pd.read_csv(args.a, sep='\t')
		if not os.path.exists(args.q):
			raise ArgsError('The qtl file {} is not exists.'.format(args.q))
		qtl = pd.read_csv(args.q, sep=',')
		qtl_anno = mw.qtl_anno(qtl, anno)
		qtl_anno.to_csv(args.o+'.qtl_anno.csv',index=False)
		log.log('QTL annotation is done.')

	if args.visual:
		log.log('Begin plot manhattan-plot and qqplot...')
		input = 'output'
		if args.i:
			input = args.i
		if not os.path.exists(input):
			raise ArgsError('The gwas output directory {} is not exists.'.format(input))
		thread = 1
		if args.thread:
			thread = args.thread
		if os.path.exists(args.o+'/manhattan_plot'):
			shutil.rmtree(args.o+'/manhattan_plot')
		if os.path.exists(args.o+'/qqplot'):
			shutil.rmtree(args.o+'/qqplot')
		os.makedirs(args.o+'/manhattan_plot')
		os.makedirs(args.o+'/qqplot')
		s = mw.gwas_plot_parallel(1e-2, thread, 'pdf', input)
		for fn in glob.glob('Manhattan.*pdf'):
			if fn != 'Manhattan.multi_trait.pdf':
				shutil.move(fn, args.o+'/manhattan_plot')
		for fn in glob.glob('QQplot.*pdf'):
			shutil.move(fn, args.o+'/qqplot')
		if args.multi:
			if not args.q:
				raise ArgsError('Parameters -q must be set!')
			log.log('Begin plot multi-trait manhattan plot...')
			qtl = pd.read_csv(args.q, sep=',')
			mw.multi_trait_plot(input, qtl, 'multi_trait', 'pdf')
			shutil.move('./Manhattan.multi_trait.pdf', args.o)
		log.log('GWAS Plot is done.')

	if args.peaktest:
		if not (args.p and args.q and args.g and args.o):
			raise ArgsError('Parameters -p, -g, -q and -o must be set!')
		log.log('Begin peak haplptype test and plot boxplot...')
		if os.path.exists(args.o+'/boxplot'):
			shutil.rmtree(args.o+'/boxplot')
		os.makedirs(args.o+'/boxplot')
		phe = pd.read_csv(args.p, index_col=0)
		qtl = pd.read_csv(args.q, sep=',')
		g = mw.read_genotype(args.g)
		mw.boxplot(phe, g, qtl, args.o)
		for fn in glob.glob('*boxplot.pdf'):
			shutil.move(fn, args.o+'/boxplot')
		log.log('Peak haplptype test done and boxplots are generated.')

	if args.genetest:
		if not (args.f1 and args.f2 and args.vcf and args.o):
			raise ArgsError('Parameters -f1, -f2, -vcf, -p and -o must be set!')
		log.log('Begin gene haplptype test and plot...')
		if os.path.exists(args.o+'/gene_haplotest'):
			shutil.rmtree(args.o+'/gene_haplotest')
		os.makedirs(args.o+'/gene_haplotest')
		b = mw.genetest(args.f1, args.f2, args.vcf, args.p, args.o)
		if not b:
			log.log('Error when perform gene haplotype test!')
		else:
			log.log('Successfully performed gene haplotype test.')


def mr(args, log):
	log.log('begin Mendelian Randomization analysis...')
	if args.tf is not None:
		if not os.path.exists(args.gene_exp):
			raise ArgsError('the gene expression file {} is not exists.'.format(args.gene_exp))
		gene_exp = pd.read_csv(args.gene_exp, index_col=0)
		if gene_exp.empty:
			log.log('the gene expression file {} is empty and software is terminated.'.format(args.gene_exp))
			sys.exit()
		tf = pd.read_csv(args.tf, sep='\t')
		target = pd.read_csv(args.target, sep='\t')
		mTrait = gene_exp.loc[:, gene_exp.columns.isin(tf.gene_id)]
		pTrait = gene_exp.loc[:, gene_exp.columns.isin(target.gene_id)]
	else:
		if args.pairwise:
			args.exposure = args.gene_exp
			args.outcome = args.gene_exp
		if not os.path.exists(args.exposure):
			raise ArgsError('the exposure file {} is not exists.'.format(args.exposure))
		mTrait = pd.read_csv(args.exposure, index_col=0)
		if mTrait.empty:
			log.log('the exposure file {} is empty and software is terminated.'.format(args.exposure))
			sys.exit()
		if not os.path.exists(args.outcome):
			raise ArgsError('the outcome file {} is not exists.'.format(args.outcome))
		pTrait = pd.read_csv(args.outcome, index_col=0)
		if pTrait.empty:
			log.log('the outcome file {} is empty and software is terminated.'.format(args.outcome))
			sys.exit()
	if not os.path.exists(args.qtl):
		raise ArgsError('the exposure QTL file {} is not exists.'.format(args.qtl))
	qtl = pd.read_csv(args.qtl)
	if qtl.empty:
		log.log('the exposure QTL file {} is empty and software is terminated.'.format(args.qtl))
		sys.exit()
	mTrait = mTrait[qtl.phe_name.unique()]
	if not os.path.exists(args.g+'.bed'):
		raise ArgsError('the plink bed format genotype file {} is not exists.'.format(args.g))
	if args.lm and args.mlm:
		log.log('please select linear model or mixed linear model for Mendelian Randomization analysis')
		log.log('software is terminated.')
		sys.exit()
	if args.lm:
		log.log('perform Mendelian Randomization through linear model')
		g = mmr.read_genotype(args.g)
		if g is None:
			raise FileTypeError('{0}.bed is not a plink bed format file or {0}.bim file and {0}.fam file is not in the same folder with {0}.bed, or {0}.bed is empty.'.format(args.g))
		g = g.where(g.snp.isin(qtl.SNP), drop=True)
		g = pd.DataFrame(g.values, index=g.sample, columns=g.snp.values)
		ril = g.index.intersection(mTrait.index.intersection(pTrait.index))
		g = g.reindex(ril)
		mTrait = mTrait.reindex(ril)
		pTrait = pTrait.reindex(ril)
		res = mmr.MR_parallel(qtl, mTrait, pTrait, g, args.thread, args.pvalue)
		if args.type == 'direct':
			coloc_df = mmr.target_type(qtl, tf, target)
			res_tmp = pd.DataFrame()
			for tf_id in coloc_df.phe_name.unique():
				res_tmp = pd.concat([res_tmp, res.loc[(res.mTrait == tf_id) & (res.pTrait.isin(coloc_df.loc[coloc_df.phe_name == tf_id, 'phe_name_b'])), :]])
			res = res_tmp
		res.to_csv(args.o + '.MR.csv', index=False)
		log.log('Successfully perform Mendelian randomization analysis using linear model')
	if args.mlm:
		log.log('perform Mendelian Randomization through mixed linear model')
		mmr.generate_geno_batch(qtl, mTrait, pTrait, args.g, args.thread, 'tmp_mr_bed', 'tmp_mr_rs')
		mmr.calc_MLM_effect('tmp_mr_bed', pTrait, args.thread, args.g)
		mTrait_effect, pTrait_effect, pTrait_se = mmr.get_MLM_effect_parallell('./output', mTrait, pTrait, args.thread)
		res = mmr.MR_MLM_parallel(qtl, mTrait_effect, pTrait_effect, pTrait_se, args.thread, args.pvalue)
		if args.type == 'direct':
			coloc_df = mmr.target_type(qtl, tf, target)
			res_tmp = pd.DataFrame()
			for tf_id in coloc_df.phe_name.unique():
				res_tmp = pd.concat([res_tmp, res.loc[(res.mTrait == tf_id) & (res.pTrait.isin(coloc_df.loc[coloc_df.phe_name == tf_id, 'phe_name_b'])), :]])
			res = res_tmp
		res.to_csv(args.o+'.MR.csv', index=False)
		log.log('Successfully perform Mendelian randomization analysis using mixed linear model')
	log.log('Mendelian Randomization analysis is successed')


def net(args, log):
	log.log('begin read pairwise Mendelian Randomization analysis file')
	if not os.path.exists(args.mr):
		raise ArgsError('the pairwise Mendelian Randomization analysis file {} is not exists.'.format(args.mr))
	mr = pd.read_csv(args.mr)
	if mr.empty:
		log.log('the pairwise Mendelian Randomization analysis file {} is empty and software is terminated.'.format(args.mr))
		sys.exit()
	log.log('begin calculate edge weight...')
	ew = mnet.get_weight(mr)
	ew.columns = ['row', 'col', 'weight']
	log.log('begin network analysis and identify module...')
	ew.to_csv(args.o+'.edge_list', sep='\t', index=False, header=None)
	cluster_one_res = mnet.module_identify(args.o+'.edge_list', args.module_size)
	if cluster_one_res is None:
		log.log('cluster_one-1.0.jar file is not in /pwd/MRBIGR/utils, network analysis can not be run.')
		sys.exit()
	log.log('begin hub node analysis of each module')
	cluster_one_res, hub_res = mnet.hub_identify(ew, cluster_one_res)
	cluster_one_res.to_csv(args.o+'.module.csv', index=False)
	if args.plot:
		log.log('begin plot network of each module')
		mnet.module_network_plot(ew, cluster_one_res, hub_res, args.o)
	log.log('Network analysis is done.')


def go(args, log):
	if not os.path.exists(args.gene_lists):
		raise ArgsError('the gene lists file {} is not exists.'.format(args.gene_lists))
	gene_lists = pd.read_csv(args.gene_lists)
	if gene_lists.empty:
		log.log('the gene lists file {} is empty and software is terminated.'.format(args.gene_lists))
		sys.exit()
	if args.go_anno is None:
		if not mgo.cheak_orgdb('org.Custom.eg.db'):
			log.log('No GO annotation database provided, please supply GO annotation file for database construction and software is terminated.')
			sys.exit()
		log.log('begin GO enrichment analysis...')
		res = mgo.go_enrich(gene_lists, 'org.Custom.eg.db', args.pvalue)
		if res is None:
			log.log('GO enrichment analysis is failed.')
			sys.exit()
		else:
			log.log('GO enrichment analysis is complete.')
	else:
		if not os.path.exists(args.go_anno):
			raise ArgsError('the GO annotation file {} is not exists.'.format(args.go_anno))
		go_anno = pd.read_csv(args.go_anno, sep='\t', header=None)
		if go_anno.empty:
			log.log('the GO annotation file {} is empty and software is terminated.'.format(args.go_anno))
			sys.exit()
		if not os.path.exists(args.gene_info):
			raise ArgsError('the gene information file {} is not exists.'.format(args.gene_info))
		gene_info = pd.read_csv(args.gene_info, sep='\t')
		if gene_info.empty:
			log.log('the gene information file {} is empty and software is terminated.'.format(args.gene_info))
			sys.exit()
		log.log('begin construct GO annotation database...')
		s = mgo.build_orgdb(go_anno, gene_info)
		if not s:
			log.log(
				'Failed to build GO annotation database, please check GO annotation file or gene information file for database construction and software is terminated.')
			sys.exit()
		log.log('Successfully constructed GO database.')
		log.log('begin GO enrichment analysis...')
		res = mgo.go_enrich(gene_lists, 'org.Custom.eg.db', args.pvalue)
		if res is None:
			log.log('GO enrichment analysis is failed.')
			sys.exit()
		else:
			log.log('GO enrichment analysis is complete.')
	res.to_csv(args.o+'.GO.csv', index=False)
	log.log('begin visualize GO enrichment result...')
	s = mgo.go_plot(res, args.plot_type.split(','), args.o)
	if s:
		log.log('Successful visualization of GO enrichment results')
	else:
		log.log('Visualization of GO enrichment results failed')


def plot(args, log):
	log.log('begin reading file used in plotting...')
	if args.plot_type == 'manhattan':
		log.log('begin manhattan plot')
		s = mplt.manhattan_plot(args.i, args.o)
	if args.plot_type == 'qqplot':
		log.log('begin QQ plot')
		s = mplt.qq_plot(args.i, args.o)
	if args.plot_type == 'SNPdensity':
		log.log('begin SNP density plot')
		s = mplt.snp_density_plot(args.i, args.o)
	if not os.path.exists(args.i):
		raise ArgsError('the file used in plotting {} is not exists.'.format(args.i))
	i = pd.read_csv(args.i)
	if i.empty:
		log.log('the file used in plotting {} is not exists'.format(args.i))
		sys.exit()
	i.columns = [c.replace('/', '.') for c in i.columns]
	i.columns = ['X'+c if re.match(r'\d+', c) else c for c in i.columns]
	if args.order is not None:
		log.log('begin reading order file...')
		if not os.path.exists(args.order):
			log.log('warnings:the order file {} is not exists, using default order to plot.'.format(args.order))
			args.order = None
		else:
			order = pd.read_csv(args.order, header=None)
			if order.empty:
				log.log('warnings:the order file {} is empty, using default order to plot.'.format(args.order))
				args.order = None
	if args.group is not None:
		log.log('begin reading group file...')
		if not os.path.exists(args.group):
			log.log('warnings:the group file {} is not exists, plotting without group classification.'.format(args.group))
			group = pd.DataFrame()
		else:
			group = pd.read_csv(args.group)
			if group.empty:
				log.log('warnings:the group file {} is empty, plotting without group classification.'.format(args.group))
	else:
		group = pd.DataFrame()
	if args.plot_type == 'forestplot':
		log.log('begin plotting Mendelian Randomization forest plot')
		for mTrait in i['mTrait'].unique():
			mr_sub = i.loc[i.mTrait == mTrait, :]
			for snp in mr_sub['snp'].unique():
				d = mr_sub.loc[mr_sub.snp == snp, :]
				if args.order is not None:
					d.index = d['pTrait']
					d = d.reindex(order.iloc[:, 0].values)
				s = mplt.forest_plot(d, mTrait, snp, args.o)
				if not s:
					log.log('mTrait {} with snp {} forest plot is failed'.format(mTrait, snp))
	if args.plot_type == 'scatter_mr':
		log.log('begin plotting Mendelian Randomization scatter plot')
		for mTrait in i['mTrait'].unique():
			mr_sub = i.loc[i.mTrait == mTrait, :]
			for snp in mr_sub['snp'].unique():
				d = mr_sub.loc[mr_sub.snp == snp, :]
				s = mplt.scatter_plot_mr(d, group, mTrait, snp, args.o)
				if not s:
					log.log('mTrait {} with snp {} scatter plot is failed'.format(mTrait, snp))
	if args.plot_type == 'hist':
		log.log('begin histogram plot')
		s = mplt.hist_plot(i, args.o)
	if args.plot_type == 'scatter_ps':
		log.log('begin population structure scatter plot')
		s = mplt.scatter_plot_ps(i, group, args.o)
	if args.plot_type == 'boxplot':
		log.log('begin boxplot plot')
		s = mplt.box_plot(i, group, args.o)
	if args.plot_type in ['barplot', 'dotplot']:
		s = mgo.go_plot(i, args.plot_type, args.o)
	if args.plot_type not in ['forest_plot', 'scatter_plot_mr'] and s:
		log.log('plotting is done')


def MRBIGR(args):
	if args.command == 'geno':
		log = Logger(args.o+'.geno.log.txt')
		log.log(sys.argv)
		geno(args, log)
	if args.command == 'pheno':
		log = Logger(args.o+'.phen.log.txt')
		log.log(sys.argv)
		pheno(args, log)
	if args.command == 'gwas':
		log = Logger(args.o+'.gwas.log.txt')
		log.log(sys.argv)
		gwas(args, log)
	if args.command == 'mr':
		log = Logger(args.o + '.mr.log.txt')
		log.log(sys.argv)
		mr(args, log)
	if args.command == 'net':
		log = Logger(args.o + '.net.log.txt')
		log.log(sys.argv)
		net(args, log)
	if args.command == 'go':
		log = Logger(args.o+'.go.log.txt')
		log.log(sys.argv)
		go(args, log)
	if args.command == 'plot':
		log = Logger(args.o + '.plot.log.txt')
		log.log(sys.argv)
		plot(args, log)


USAGE = ' MRBIGR: Mendelian Randomization-Based Inference of Genetic Regulation -- a toolkit for pre-GWAS, GWAS and especially, post-GWAS analysis '
parser = argparse.ArgumentParser(usage='%(prog)s command [options]', description=USAGE)
subparsers = parser.add_subparsers(title='command', metavar='', dest='command', prog=parser.prog)

parser_geno = subparsers.add_parser('geno', help='Genotype file based analyses', usage='%(prog)s [options]')
parser_geno.add_argument('-g', help='Genotype file of plink-bed format. Prefix of plink-bed file is needed, that is, when the input files are name.bed/name.bim/name.fam, -g should bed set to name. It should be used together with these parameters: –convert/-qc/-imput/-clump/-kinship/-pca/-tsne/-tree.')
parser_geno.add_argument('-o', default='geno_output', metavar='[geno_output]', help='Prefix of the output files, the default value is geno_output. It should be used together with these parameters: –convert/-qc/-imput/-clump/-kinship/-pca/-tsne/-tree/-anno.')
parser_geno.add_argument('-convert', action='store_true', help='Convert hapmap or vcf format genotype file to plink-bed format files, or convert plink-bed format genotype files to vcf format. It should be used together with the parameters –g/-hapmap/-vcf and –o.')
parser_geno.add_argument('-hapmap', help='The hapmap format genotype file to be converted to plink-bed format. It should be used together with the parameter –convert.')
parser_geno.add_argument('-vcf', help='The vcf format genotype file to be converted to plink-bed format. It should be used together with the parameter –convert.')
parser_geno.add_argument('-qc', action='store_true', help='Quality control of the original genotype data. SNPs and individuals pass the –maf, -mis and –mind cutoff would be kept. It should be used together with the parameters –g/-o/-maf/-mis/-mind.')
parser_geno.add_argument('-maf', type=float, metavar='[0.05]', help='Minimum Minor Allele Frequency (MAF) for a SNP to be kept, default is 0.05. It should be used together with the parameter –qc.')
parser_geno.add_argument('-mis', type=float, metavar='[0.1]', help='Maximum proportion of missing values for a SNP to be kept, default is 0.1. It should be used together with the parameter –qc.')
parser_geno.add_argument('-mind', type=float, metavar='[0.99]', help='Maximum proportion of missing values for a sample to be kept, default is 0.99. It should be used together with the parameter –qc.')
parser_geno.add_argument('-imput', action='store_true', help='Perform fast genotype imputation via mode, mean, sampling according to allele frequencies, or 0. It should be used together with the parameters –g/-o/–method.')
parser_geno.add_argument('-method', default='mode', metavar='[mode]', help='Genotype imputation method. Options: "mode" (most frequent call), "random" (sampling according to allele frequencies), "mean0" (rounded mean), "mean2" (rounded mean to 2 decimal places). It should be used together with the parameters –g/-o/–imput.')
parser_geno.add_argument('-clump', action='store_true', help='Clumping analysis is used to keep only one representative SNP per region of LD, to reduce the dimensionality of the original genotype file. It should be used together with the parameters –g/-o/-r2.')
parser_geno.add_argument('-r2', default=0.8, type=float, metavar='[0.8]', help=' r^2 value for the SNP clumping analysis, default is 0.8. It should be used together with the parameters –clump.')
parser_geno.add_argument('-kinship', action='store_true', help='Generate a kinship/relatedness matrix from plink-bed format genotype data. It should be used together with the parameters –g/-o.')
parser_geno.add_argument('-pca', action='store_true', help='Perform principal component analysis for genotype dimensionality reduction and population structure display. It should be used together with the parameters –g/-o/-dim.')
parser_geno.add_argument('-tsne', action='store_true', help='Perform t-SNE analysis to display the population structure. It should be used together with the parameters –g/-o/-dim.')
parser_geno.add_argument('-dim', type=int, metavar='[10|3]', help='Dimensions for the output of PCA or t-SNE analysis results, default is 10 for PCA and 3 for t-SNE. It should be used together with the parameters –pca/-tsne.')
parser_geno.add_argument('-tree', action='store_true', help='To build a ML-tree from plink-bed format genotype data. Clumping analysis or genotype pruning is recommended before this step for large dataset. It should be used together with the parameters –g/-o.')
parser_geno.add_argument('-anno', action='store_true', help='To annotate the vcf-format genotype data using ANNOVAR. It should be used together with the parameters -vcf /-o/–db/-dbdir/-gtf/-fa. ')
parser_geno.add_argument('-db', help='The generated or input annotation database name. It should be used together with the parameter -anno.')
parser_geno.add_argument('-dbdir', help='The generated/input annotation database directory name, with *_refGene.txt and *_refGeneMrna.fa files in this directory. It should be used together with the parameter -anno. ')
parser_geno.add_argument('-gtf', help='GTF format gene annotation file, used to generate the annotation database. It should be used together with the parameter -anno.')
parser_geno.add_argument('-fa', help='Fasta format genomic reference sequence file, used to generate the annotation database. It should be used together with the parameter -anno.')

parser_pheno = subparsers.add_parser('pheno', help='Phenotype file based analyses', usage='%(prog)s [options]')
parser_pheno.add_argument('-p', help='Phenotype file of *.csv format for processing and transformation, the first column and the first row are the names of samples and traits, respectively. It should be used together with the parameters –qc/-imput/-scale/-correct/-merge.')
parser_pheno.add_argument('-o', default='pheno_output', metavar='[pheno_output]', help='Prefix of the output files, default is pheno_output. It should be used together with the parameters –qc/-imput/-scale/-correct/-merge.')
parser_pheno.add_argument('-qc', action='store_true', help='Perform phenotype data quality control and filtration, this includes outlier data removal, missing rate filtering and small value filtering, each of the criterions is optional and can be skipped, while only the data points pass the cutoff would be kept. It should be used together with the parameters –p/-o/-rout/-mis/-val.')
parser_pheno.add_argument('-tr', action='store_true', help='Transposition of the input data')
parser_pheno.add_argument('-rout', metavar='[zscore]', help='Outlier removal of phenotype data, the default method is zscore. Options: zscore, boxplot. It should be used together with the parameter –qc.')
parser_pheno.add_argument('-mis', type=float, metavar='[0.5]', help='Filter traits with a missing rate larger than the cutoff, default is 0.5. It should be used together with the parameter –qc.')
parser_pheno.add_argument('-val', type=float, metavar='[0.1]', help='Filter traits with an average value smaller than the cutoff, default is 0.1. It should be used together with the parameter –qc.')
parser_pheno.add_argument('-imput', action='store_true', help='Perform imputation for missing phenotype values. It should be used together with the parameters –p/-o/-method.')
parser_pheno.add_argument('-method', default='mean', metavar='[mean]', help='Phenotype values imputation method. Options: mean, median, most_frequent. It should be used together with the parameter -imput.')
parser_pheno.add_argument('-scale', action='store_true', help='Scaling/Normalization/Standardization/Transformation of the phenotype data. We provide a variety of methods for phenotype data scaling. It should be used together with the parameters –g/-o and any one or more of the parameters –log2/-log10/-ln/-boxcox/-qqnorm/-minmax/-zscore/-robust.')
parser_pheno.add_argument('-log2', action='store_true', help='log2 scale of the phenotype data. It should be used together with the parameter -scale.')
parser_pheno.add_argument('-log10', action='store_true', help='log10 scale of the phenotype data. It should be used together with the parameter -scale.')
parser_pheno.add_argument('-ln', action='store_true', help='ln scale of the phenotype data. It should be used together with the parameter -scale.')
parser_pheno.add_argument('-boxcox', action='store_true', help='Box-Cox transformation of the phenotype data. It should be used together with the parameter -scale.')
parser_pheno.add_argument('-qqnorm', action='store_true', help='Normal quantile transformation of the phenotype data. It should be used together with the parameter -scale.')
parser_pheno.add_argument('-minmax', action='store_true', help='Min-Max normalization of the phenotype data. It should be used together with the parameter -scale.')
parser_pheno.add_argument('-zscore', action='store_true', help='Z-score Standardization of the phenotype data. It should be used together with the parameter -scale.')
parser_pheno.add_argument('-robust', action='store_true', help='Scale phenotype values using statistics that are robust to outliers. It should be used together with the parameter -scale.')
parser_pheno.add_argument('-correct', action='store_true', help='Correct phenotype data caused by population structure through a genotype-based pre-calculated PCA file. It should be used together with the parameters –p/-o/-pca.')
parser_pheno.add_argument('-pca', help='A pre-calculated PCA file used for population structure correction of the phenotype data. It should be used together with the parameter –correct.')
parser_pheno.add_argument('-merge', action='store_true', help='Merge the input phenotype values from different environment using either mean values, BLUP or BLUE methods. Only one trait is allowed in the input phenotype file. It should be used together with the parameters –p/-o/-mm.')
parser_pheno.add_argument('-mm', default='mean', metavar='[mean]', help='Merge method of phenotype values, default is mean. Options: mean, blup, blue. It should be used together with the parameter –merge.')

parser_gwas = subparsers.add_parser('gwas', help='GWAS and QTL related analyses', usage='%(prog)s [options]')
parser_gwas.add_argument('-g', help='Genotype file of plink-bed format. Prefix of plink-bed file is needed, that is, when the input files are name.bed/name.bim/name.fam, -g should bed set to name. It should be used together with the parameters -gwas/-peaktest.')
parser_gwas.add_argument('-p', help=' Phenotype file of *.csv format after processing and transformation, the first column and the first row are the names of samples and traits, respectively. It should be used together with the parameters –gwas/-peaktest/-genetest.')
parser_gwas.add_argument('-o', default='gwas_output', metavar='[pheno_output]', help='Prefix of the output files. It should be used together with the parameters -qtl/-anno/-visual/-peaktest/-genetest.')
parser_gwas.add_argument('-thread', default=1, type=int, metavar='[1]', help=' Number of threads for GWAS and QTL analysis, default is 1. It should be used together with the parameters –gwas/-qtl/-visual.')
parser_gwas.add_argument('-gwas', action='store_true', help='Perform GWAS using a linear mixed model (lmm) or a linear model (lm) implemented in GEMMA. All the result files will be generated in a newly-built output directory located in the working directory by default. It should be used together with the parameters –g/-p/–model/-thread.')
parser_gwas.add_argument('-model', default='lmm', metavar='[lmm]', help='Fit a linear mixed model (lmm) or a linear model (lm), default is lmm. Options: lmm, lm. It should be used together with the parameter –gwas.' )
parser_gwas.add_argument('-qtl', action='store_true', help='QTL detection from the GWAS results based on the PLINK-clump method. It should be used together with the parameter –i/-p1/-p2/-p2n/-thread.')
parser_gwas.add_argument('-i', default='output', metavar='[output]', help='Parent directory of the GWAS output files *.assoc.txt, default is output. It should be used together with the parameter –qtl.')
parser_gwas.add_argument('-p1', default=1e-7, type=float, metavar='[1e-7]', help='Significance threshold for index SNPs used to determine QTLs, default is 1e-7. It should be used together with the parameter –qtl.')
parser_gwas.add_argument('-p2', default=1e-5, type=float, metavar='[1e-5]', help='Secondary significance threshold for clumped SNPs used to determine the reliability of QTLs, default is 1e-5. It should be used together with the parameter –qtl.')
parser_gwas.add_argument('-p2n', default=5, type=int, metavar='[5]', help='Secondary significance SNP number in a QTL, default is 5. It should be used together with the parameter –qtl.')
parser_gwas.add_argument('-anno', action='store_true', help='Perform QTL annotation using a pre-formatted gene annotation file. It should be used together with the parameter –q/-a.')
parser_gwas.add_argument('-q', help='QTL result file in *.csv format generated by the -qtl subcommand. It should be used together with the parameter –anno.')
parser_gwas.add_argument('-a', help='Pre-formatted gene annotation file in used for QTL annotation. It should be used together with the parameter –anno.')
parser_gwas.add_argument('-visual', action='store_true', help='Generate Manhattan-plots and QQ-plots for visualization of the GWAS results. It should be used together with the parameters –i/-multi/-thread.')
parser_gwas.add_argument('-multi', action='store_true', help='Generate multiple-trait Manhattan-plot, -q must be set. It should be used together with the parameters –visual/-q.')
parser_gwas.add_argument('-peaktest', action='store_true', help='Perform peak-based haplotype test using lead SNP in each QTL. It should be used together with the parameters –p/-g/-q/-o.')
parser_gwas.add_argument('-genetest', action='store_true', help='Perform gene-based haplotype test for each trait. It should be used together with the parameters -f1/-f2/-vcf/-p/-o.')
parser_gwas.add_argument('-f1', help='QTL annotation file in *.csv format generated by the subcommand gwas –anno. It should be used together with the parameter –genetest.')
parser_gwas.add_argument('-f2', help='Genotype annotation file in *.bed format generated by the subcommand geno –anno. It should be used together with the parameter –genetest.')
parser_gwas.add_argument('-vcf', help='VCF format genotype file. It should be used together with the parameter –genetest.')

parser_mr = subparsers.add_parser('mr', help='perform Mendelian Randomization analysis', usage='%(prog)s [options]')
parser_mr.add_argument('-g', metavar='', help='Genotype file in plink bed format')
parser_mr.add_argument('-exposure', metavar='', default=None,
					   help='Exposure phenotype file, such as mTrait phenotype, the file format is csv format, the first column and the first row are the names of inbred lines and phenotypes respectively')
parser_mr.add_argument('-outcome', metavar='', default=None,
					   help='Outcome phenotype file, such as pTrait phenotype, the file format is csv format, the first column and the first row are the names of inbred lines and phenotypes respectively')
parser_mr.add_argument('-tf', metavar='', default=None, help='Transcription factors annotation file')
parser_mr.add_argument('-target', metavar='', default=None, help='targeted genes annotation file')
parser_mr.add_argument('-type', metavar='[All]', default='All', help='Regulatory relationship between transcription factors and targeted genes. direct, direct regulatory relationship between transcription factors and targeted genes. All, direct and indirect regulatory relationship between transcription factors and targeted genes')
parser_mr.add_argument('-pairwise', action='store_true', help='perform Mendelian randomization analysis between all gene pairs')
parser_mr.add_argument('-gene_exp', metavar='', help='gene expression file, the file format is csv format, the first column and the first row are the names of inbred lines and phenotypes respectively ')
parser_mr.add_argument('-qtl', metavar='', help='Exposure phenotype QTL file, generated by subcommand regiongwas')
parser_mr.add_argument('-lm', action='store_true', help='Perform Mendelian Randomization through linear model')
parser_mr.add_argument('-mlm', action='store_true', help='Perform Mendelian Randomization through mixed linear model')
parser_mr.add_argument('-pvalue', default=1, type=float, metavar='[default:1]',
					   help='The pvalue cutoff of Mendelian randomization analysis result output,default 1')
parser_mr.add_argument('-thread', default=1, type=int, metavar='[1]',
					   help='Number of threads for Mendelian Randomization analysis')
parser_mr.add_argument('-o', default='mr_out', metavar='[mr_out]',
					   help='The prefix of the output file, default mr_out')

parser_net = subparsers.add_parser('net', help='perform network analysis', usage='%(prog)s [options]')
parser_net.add_argument('-mr', metavar='', help='pairwise Mendelian Randomization analysis result')
parser_net.add_argument('-module_size', metavar='[5]', default=5, type=int, help='minimal size of genes in a module,default 5')
parser_net.add_argument('-plot', action='store_true', help='plot network of each module')
parser_net.add_argument('-o', default='net_out', metavar='[net_out]', help='The prefix of the output file, default go_out')

parser_go = subparsers.add_parser('go', help='perform GO enrichment analysis', usage='%(prog)s [options]')
parser_go.add_argument('-gene_lists', metavar='', help='gene lists file for GO enrichment analysis')
parser_go.add_argument('-go_anno', metavar='', default=None, help='GO annotation file for each gene used to construct GO annotation database')
parser_go.add_argument('-gene_info', metavar='', help='gene information file used to construct GO annotation database')
parser_go.add_argument('-pvalue', metavar='', default=0.05, type=float, help='adjusted pvalue cutoff on enrichment tests to report')
parser_go.add_argument('-plot_type', default='dotplot', metavar='[dotplot]', help='Visualization types of GO enrichment result, support barplot, dotplot, Gene-Concept Network(cnetplot), heatplot, Enrichment Map(emapplot) and upsetplot, one or more plot types can be specified at the same time, multiple types are separated by commas, such as barplot,emapplot,heatplot')
parser_go.add_argument('-o', default='go_out', metavar='[go_out]', help='The prefix of the output file, default go_out')

parser_plot = subparsers.add_parser('plot', help='Visualize the results generated by MRBIGR and the data used in MRBIGR by various types of plots')
parser_plot.add_argument('-i', metavar='', help='the result file generated by MRBIGR or the data file used in MRBIGR, the file format is csv format')
parser_plot.add_argument('-plot_type', metavar='', help='Visualization type of the result file generated by MRBIGR or the data file used in MRBIGR,including scatter_mr, forestplot, manhattan, qqplot, SNPdensity, hist, boxplot, scatter_ps, barplot, dotplot; support scatterplot and forest plot for Mendelian Randomization analysis, manhattan plot and QQ plot for GWAS result, SNP density plot for GWAS result and SNP info file, boxplot and histogram for phenotype data, barplot and dotplot for GO enrich result, scatter plot for population structure.')
parser_plot.add_argument('-order', metavar='', default=None, help='Visualizing the results of Mendelian randomization in the order of trait in the file, uesed in forest plot')
parser_plot.add_argument('-group', metavar='', default=None, help='group file for scatter plot and boxplot')
parser_plot.add_argument('-o', default='plot_out', metavar='[plot_out]', help='The prefix of the output file, default plot_out')


if __name__ == "__main__":
	if len(sys.argv) <= 2:
		sys.argv.append('-h')
	args = parser.parse_args()
	MRBIGR(args)
