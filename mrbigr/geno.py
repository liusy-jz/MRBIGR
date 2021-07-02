import os
import pandas as pd
import numpy as np
from sklearn.manifold import TSNE
from multiprocessing import cpu_count
from collections import Counter
import struct
import re

import warnings
warnings.filterwarnings("ignore")



def convert(g, n):
	c = [1, 4, 16, 64]
	s = 0
	for i in range(len(g)):
		if g[i][0] != g[i][1]:
			s += 2 * c[i]
		elif g[i] == 'NN':
			s += 1 * c[i]
		elif g[i][0] == n[0]:
			s += 0 * c[i]
		else:
			s += 3 * c[i]
	return s

def vcf2plink_bed(vcf, out):
	try:
		os.system('plink --vcf {0} --out {1} --double-id --biallelic-only strict '.format(vcf,out))
	except Exception:
		return False
	else:
		return True

def plink_bed2vcf(bed, out):
	try:
		os.system('plink --bfile {0} --recode vcf-iid --out {1} '.format(bed,out))
	except Exception:
		return False
	else:
		return True



def hap2plink_bed(hap, plink):
	try:
		with open(hap) as h, open(plink + '.bed', 'wb') as b, open(plink + '.bim', 'w') as bim, open(plink + '.fam','w') as fam:

			b.write(struct.pack('B', 108))
			b.write(struct.pack('B', 27))
			b.write(struct.pack('B', 1))
			out = list()
			for _, l in enumerate(h):
				l = l.strip().split('\t')
				if _ == 0:
					samples = l[11:]
					for s in samples:
						fam.write(' '.join([s, s, '0', '0', '0', '-9']) + '\n')
				else:
					g_seq = ''.join(l[11:])
					nlu_num = list(set(g_seq))
					if len(nlu_num) > 2 and 'N' not in nlu_num:
						print("Warning: There is more than two nucleotide letter snp in hapmap file. skip this snp")
						continue
					if len(nlu_num) == 1:
						print("Warning: There is one nucleotide letter snp in hapmap file. skip this snp")
						continue
					if 'N' in nlu_num:
						c = Counter(g_seq)
						del c['N']
						n = sorted(c, key=lambda x: c[x])
					else:
						g_seq_sort = ''.join(sorted(g_seq))
						major = g_seq_sort[len(g_seq) // 2]
						if g_seq_sort[0] == major:
							n = [g_seq_sort[-1], major]
						else:
							n = [g_seq_sort[0], major]
					bim.write('\t'.join([l[2], l[0], '0', l[3]] + n) + '\n')
					for num in range(11, len(l), 4):
						if num + 4 < len(l):
							out.append(convert(l[num:num + 4], n))
						else:
							out.append(convert(l[num:len(l)], n))
			b.write(struct.pack('B'*len(out),*out))
	except Exception:
		return False
	else:
		return True


def snp_qc(bed, out, maf, mis, mind):
	try:
		os.system('plink --bfile {0} --out {1} --maf {2} --geno {3} --mind {4} --make-bed'.format(bed,out, maf, mis, mind))
	except Exception:
		return False
	else:
		return True


def snp_imput(bed, out, meth):
	try:
		from rpy2.robjects.packages import importr
		from rpy2.robjects import pandas2ri
		import rpy2.robjects as robjects
		pandas2ri.activate()
		robjects.r['options'](warn=-1)
		base = importr('base')
		bigsnpr = importr('bigsnpr')
		g = bigsnpr.snp_readBed(bed, backingfile=base.tempfile())
		g = bigsnpr.snp_attach(g)
		g[0] = bigsnpr.snp_fastImputeSimple(g[0], method=meth, ncores=1)
		bigsnpr.snp_writeBed(g, out)
	except Exception:
		return False
	else:
		return True
		
def snp_clumping(bed, out, r2):
	from rpy2.robjects.packages import importr
	from rpy2.robjects import pandas2ri
	import rpy2.robjects as robjects
	pandas2ri.activate()
	robjects.r['options'](warn=-1)
	base = importr('base')
	bigsnpr = importr('bigsnpr')
	g = bigsnpr.snp_readBed(bed, backingfile=base.tempfile())
	g = bigsnpr.snp_attach(g)
	snp_keep = bigsnpr.snp_clumping(g[0], infos_chr=g[2]['chromosome'], infos_pos=g[2]['physical.pos'],thr_r2=r2, ncores=1)
	g_clump = bigsnpr.subset_bigSNP(g, ind_col=snp_keep)
	g_clump = bigsnpr.snp_attach(g_clump)
	bigsnpr.snp_writeBed(g_clump, out)
	return len(g[2]), len(snp_keep)

def kinship(bed, out):
	try:
		os.system('gemma.linux -bfile {0} -gk 1 -o {1}'.format(bed, out))
	except Exception:
		return False
	else:
		return True

def pca(bed, dim):
	try:
		from rpy2.robjects.packages import importr
		from rpy2.robjects import pandas2ri
		from rpy2.rinterface_lib.embedded import RRuntimeError
		pandas2ri.activate()

		base = importr('base')
		bigsnpr = importr('bigsnpr')
		bigstatsr = importr('bigstatsr')

		base.sink('/dev/null')
		g = bigsnpr.snp_readBed(bed, backingfile=base.tempfile())
		g = bigsnpr.snp_attach(g)
		svd = bigsnpr.snp_autoSVD(g[0], infos_chr=g[2]['chromosome'], infos_pos=g[2]['physical.pos'], thr_r2=np.nan, ncores=1, k=dim)
		base.sink()
	except RRuntimeError:
		return None
	else:
		pc = bigstatsr.predict_big_SVD(svd)
		pc = base.data_frame(pc, row_names=g[1]['sample.ID'])
		pc.columns = ['PC' + str(i) for i in range(1, dim+1)]
		return pc


def tsne(bed, dim):
	data = pca(bed, 50)
	tsne = TSNE(n_components=dim, learning_rate=100)
	redim_array = tsne.fit_transform(data)
	df = pd.DataFrame(redim_array)
	df.index = data.index
	df.columns = ['TSNE' + str(i) for i in range(1, dim+1)]
	return pd.DataFrame(df)


def ped2fasta(ped, out):
	fi = open(ped)
	fo = open(out,'w')
	for line in fi.readlines():
		tmp = line.strip().split(' ')
		name = tmp[1]
		seq = ''.join(tmp[6:])
		seq = re.sub("0", "N", seq)
		fo.write('>'+name+"\n"+seq+"\n")
	fi.close()
	fo.close()


def fasttree(bed, out):
	try:
		os.system('plink --bfile {0} --recode --out {1}'.format(bed, out))
		ped2fasta(out+'.ped', out+'.fasta')
		os.system('FastTree -nt -gtr -quiet {0}.fasta >{1}.tree.nwk'.format(out, out))
	except Exception:
		return False
	else:
		return True

def geno_anno(vcf, dbdir, db, o):
	try:
		os.system('tableAnnovar.sh {0} {1} {2} {3}'.format(vcf, dbdir, db, o))
	except Exception:
		return False
	else:
		return True