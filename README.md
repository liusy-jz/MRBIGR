# MRBIGR
MRBIGR: Mendelian Randomization-Based Inference of Genetic Regulation

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
Rscript -e 'install.packages(c("data.table", "ggplot2", "ggsignif", "ggrepel","Matrix", "igraph", "network", "GGally", "sna"), repos="https://cloud.r-project.org")'
Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org");BiocManager::install(c("AnnotationForge","clusterProfiler"))
Rscript -e 'install.packages("bigsnpr", dependence=T, repos="https://cloud.r-project.org")'

echo "export PATH=`pwd`/utils:\$PATH" >> ~/.bashrc
source ~/.bashrc
```

### Mannual Installation
```
#Depends: R(>=3.6), python(>=3.7)

git clone https://github.com/liusy-jz/MRBIGR.git
cd MRBIGR
python setup.py build
python setup.py install

pip3 install rpy2 pyranges
Rscript -e 'install.packages(c("data.table", "ggplot2", "ggsignif", "ggrepel","Matrix", "igraph", "network", "GGally", "sna"), repos="https://cloud.r-project.org")'
Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org");BiocManager::install(c("AnnotationForge","clusterProfiler"))
Rscript -e 'install.packages("bigsnpr", dependence=T, repos="https://cloud.r-project.org")'

echo "export PATH=`pwd`/utils:\$PATH" >> ~/.bashrc
source ~/.bashrc
```

