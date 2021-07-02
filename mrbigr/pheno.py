import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import norm
from sklearn.preprocessing import MinMaxScaler, StandardScaler, RobustScaler
from sklearn.impute import SimpleImputer
import warnings

def mean(d):
	return d.mean(axis = 1, skipna = True)

def transpose(d):
	return d.T


def abundance_filter(d, abundance):
	return d.loc[:, d.mean(skipna = True) >= abundance]


def isDigit(x):
	try:
		float(x)
		return True
	except ValueError:
		return False


def missing_filter(d, missing_ratio):
	if d.applymap(np.isreal).all().sum() == d.shape[1]:
		d = d.loc[:, d.applymap(np.isnan).sum() <= d.shape[0] * missing_ratio]
		#d = d.fillna(0)
	else:
		d = d.loc[:, d.applymap(lambda x:isDigit(x)).sum() >= d.shape[0] * missing_ratio]
		d_array = d.values
		d_array[~d.applymap(lambda x:isDigit(x))] = 'NA'
		d.loc[:, :] = d_array
		d = d.astype(float)
	return d


def log2_scale(d):
	return np.log2(d+1)


def ln_scale(d):
	return np.log(d+1)


def log10_scale(d):
	return np.log10(d+1)


def boxcox_scale(d):
	d = d.apply(lambda x: stats.boxcox(x)[0])
	return d

def minmax_scale(d):
	d.loc[:,:] = MinMaxScaler().fit_transform(d.values)
	return d

def zscore_scale(d):
	d.loc[:,:] = StandardScaler().fit_transform(d.values)
	return d

def robust_scale(d):
	d.loc[:,:] = RobustScaler().fit_transform(d.values)
	return d

def ppoints(n, a=None):
	try:
		n = np.float(len(n))
	except TypeError:
		n = np.float(n)
	if a is None:
		a = 3.0/8 if(n <= 10) else 1.0/2
	return (np.arange(n) + 1 - a)/(n + 1 - 2*a)


def qqnorm(y):
	ina = np.isnan(y)
	if ina.sum() > 0:
		yN = y
		y = y[~ina]
	n = y.shape[0]
	if n == 0:
		print('y is empty or has only NAs')
		return np.array([])
	x = np.around(norm.ppf(ppoints(n)[np.argsort(np.argsort(y))]), decimals=15)
	if ina.sum() > 0:
		y = x
		x = yN
		x[~ina] = y
	return x

def pheno_imputer(d, method):
	imputer = SimpleImputer(missing_values=np.nan, strategy=method)
	d.loc[:,:] = imputer.fit_transform(d.values)
	return d

def trait_correct(pc, y):
	pc1 = pd.concat([pd.DataFrame(np.ones((y.shape[0], 1)), index=pc.index), pc], axis=1)
	vhat = np.dot(np.linalg.pinv(np.dot(pc1.T, pc1)), np.dot(pc1.T, y))
	if len(vhat.shape) == 1:
		y_corr = y - np.dot(pc, vhat[1:])
	else:
		y_corr = y - np.dot(pc, vhat[1:, :])
	return y_corr

def blup(d):
	from rpy2.robjects.packages import importr
	from rpy2.robjects import pandas2ri
	from rpy2.rinterface_lib.embedded import RRuntimeError
	import rpy2.robjects as robjects
	pandas2ri.activate()
	data_table = importr('data.table')
	base = importr('base')
	robjects.r('''Blups <- function(dat=NULL){
  library("reshape")
  library("lme4")
  # Create lmer data
  dat <- data.frame(taxa=rownames(dat), dat)
  dat[,1] <- as.character(dat[,1])
  for(i in 2:ncol(dat)){
    dat[,i] <- as.numeric(dat[,i])
  }
  datnew <- reshape::melt(dat, id=c(colnames(dat)[1]))
  
  # Create lmer formula
  termlabels <- paste("(1|", colnames(datnew)[1:2], ")", sep = "")
  
  # fit the model
  lme <- lme4::lmer(formula=reformulate(termlabels=termlabels, response=colnames(datnew)[3]), data=datnew, REML=TRUE)
  ## BLUP
  modelblup <- lme4::ranef(lme)
  BLUP_out <- as.data.frame(modelblup[[1]] + summary(lme)$coefficients[1])
  BLUP_out <- data.frame(rownames(BLUP_out), BLUP_out)
  names(BLUP_out) <- c(colnames(dat)[1], "blups")
  
  # Results
  ID <- data.frame(ID=unique(dat[,1]))
  colnames(ID) <- colnames(dat)[1]
  blup <- merge(ID, BLUP_out, all.x=T)
  res <- blup[,-1]
  rownames(res) <- blup[,1]
  res <- res[match(rownames(dat), rownames(res)),]
  return(res)
}''')
	robjects.r['options'](warn=-1)
	base.sink('/dev/null')
	Blups = robjects.r['Blups']
	res.loc[:,:] = Blups(d)
	base.sink()
	return res

def blue(d):
	from rpy2.robjects.packages import importr
	from rpy2.robjects import pandas2ri
	from rpy2.rinterface_lib.embedded import RRuntimeError
	import rpy2.robjects as robjects
	pandas2ri.activate()
	data_table = importr('data.table')
	base = importr('base')
	robjects.r('''Blues <- function(dat=NULL){
  library("reshape")
  library("lme4")
  # Create lmer data
  dat <- data.frame(taxa=rownames(dat), dat)
  dat[,1] <- as.character(dat[,1])
  for(i in 2:ncol(dat)){
    dat[,i] <- as.numeric(dat[,i])
  }
  datnew <- reshape::melt(dat, id=c(colnames(dat)[1]))
  
  # Create lmer formula
  termlabels <- c(colnames(datnew)[1], paste("(1|", colnames(datnew)[2], ")", sep = ""), "-1")
  
  # fit the model
  lme <- lme4::lmer(formula=reformulate(termlabels=termlabels, response=colnames(datnew)[3]), data=datnew, REML=TRUE)
  ## BLUE
  modelblue <- lme4::fixef(lme)
  BLUE_out <- as.data.frame(modelblue)
  BLUE_out <- data.frame(sub(colnames(dat)[1], "", rownames(BLUE_out)), BLUE_out)
  names(BLUE_out) <- c(colnames(dat)[1], "blues")
  
  # Results
  ID <- data.frame(ID=unique(dat[,1]))
  colnames(ID) <- colnames(dat)[1]
  blue <- merge(ID, BLUE_out, all.x=T)
  res <- blue[,-1]
  rownames(res) <- blue[,1]
  res <- res[match(rownames(dat), rownames(res)),]
  return(res)
}''')
	robjects.r['options'](warn=-1)
	base.sink('/dev/null')
	Blues = robjects.r['Blues']
	res.loc[:,:] = Blues(d)
	base.sink()
	return res

def outlier(d,method):
	from rpy2.robjects.packages import importr
	from rpy2.robjects import pandas2ri
	from rpy2.rinterface_lib.embedded import RRuntimeError
	import rpy2.robjects as robjects
	pandas2ri.activate()
	data_table = importr('data.table')
	base = importr('base')
	robjects.r('''Routsi <- function(dat=NULL, method="zscore"){
  library("reshape")
  # Melt data
  dat <- data.frame(taxa=rownames(dat), dat)
  dat[,1] <- as.character(dat[,1])
  for(i in 2:ncol(dat)){
    dat[,i] <- as.numeric(dat[,i])
  }
  datnew <- reshape::melt(dat, id=c(colnames(dat)[1]))
  
  # Outliers Removal
  x <- datnew[,3]
  if(method=="boxplot"){y <- which(x %in% boxplot.stats(x)$out)}
  if(method=="zscore"){upi <- mean(x,na.rm=T)+3*sd(x,na.rm=T); lowi <- mean(x,na.rm=T)-3*sd(x,na.rm=T); y <- which((x>upi)|(x<lowi))}
  if(length(y)>0){
    outlier_dat <- datnew[y,]
    datnew[y,3] <- NA
    outlier_removed_dat <- reshape::cast(datnew, as.formula(paste(colnames(datnew)[1], ' ~ ', colnames(datnew)[2], sep="")), function(x) mean(x,na.rm=T))
  } else {
    outlier_dat <- NULL
    outlier_removed_dat <- dat
  }
  res <- outlier_removed_dat[,-1]
  rownames(res) <- outlier_removed_dat[,1]
  res <- res[match(rownames(dat), rownames(res)),]
  return(res)
}''')
	robjects.r['options'](warn=-1)
	base.sink('/dev/null')
	rout = robjects.r['Routsi']
	res.loc[:,:] = rout(d,method)
	base.sink()
	return res
