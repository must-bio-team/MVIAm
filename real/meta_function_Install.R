install.packages(c('foreach', 'dplyr', 'glmnet', 'beeswarm', 'reshape2', 'ROCR', 'cba', 'pheatmap',
						 'ggplot2', 'gridExtra', 'RColorBrewer', 'pamr', 'stringi'))

source('http://bioconductor.org/biocLite.R')
biocLite(c('GEOquery', 'affy', 'annotate', 'org.Hs.eg.db', 'impute', 'sva'))

# files were obtained from http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/18.0.0/entrezg.asp
install.packages('hgu95av2hsentrezgcdf_18.0.0.tar.gz', repos=NULL, type='source')
install.packages('hgu133plus2hsentrezgcdf_18.0.0.tar.gz', repos=NULL, type='source')
