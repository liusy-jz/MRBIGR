# MRBIGR
MRBIGR: Mendelian Randomization-Based Inference of Genetic Regulation

MRBIGR is a multifunctional toolkit for pre-GWAS, GWAS and post-GWAS of both traditional and multi-omics data. MRBIGR provides all the components needed to build a complete GWAS pipeline, and integrated with rich post-GWAS analysis tools such as QTL annotation and haplotype analysis. In particular, Mendelian randomization (MR) analysis, MR-based network construction, module identification and gene ontology analysis are proposed for further genetic regulation studies. Additionally, it also produces rich plots for visualization of the analysis results and other formatted data.

## installation
### Installation using conda
```
git clone https://github.com/liusy-jz/MRBIGR.git
cd MRBIGR
conda create -n mrbigr python=3.7 -y
conda activate mrbigr
python setup.py build
python setup.py install

pip install pyranges
conda install -y -c conda-forge r-rcppeigen r-xml r-rsqlite r-europepmc r=3.6 rpy2
Rscript -e 'install.packages(c("data.table", "ggplot2", "ggsignif", "ggrepel","Matrix", "igraph", "network", "GGally", "sna","tidyr","ggraph"), repos="https://cloud.r-project.org")'
Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org");BiocManager::install(c("AnnotationForge","clusterProfiler"))'
Rscript -e 'install.packages("bigsnpr", dependence=T, repos="https://cloud.r-project.org")'

echo "export PATH=`pwd`/utils:\$PATH" >> ~/.bashrc
source ~/.bashrc
```

### Installation in local environment
```
#Depends: R(>=3.6), python(>=3.7)

git clone https://github.com/liusy-jz/MRBIGR.git
cd MRBIGR
python setup.py build
python setup.py install

pip3 install rpy2 pyranges
Rscript -e 'install.packages(c("data.table", "ggplot2", "ggsignif", "ggrepel","Matrix", "igraph", "network", "GGally", "sna","tidyr","ggraph"), repos="https://cloud.r-project.org")'
Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org");BiocManager::install(c("AnnotationForge","clusterProfiler"))'
Rscript -e 'install.packages("bigsnpr", dependence=T, repos="https://cloud.r-project.org")'

echo "export PATH=`pwd`/utils:\$PATH" >> ~/.bashrc
source ~/.bashrc
```

## Usage
Please refer to the User Mannual `MRBIGR_manual.pdf` for details.


