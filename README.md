# Multi-view based integrative analysis of gene expression data for identifying biomarkers

This repository provides some codes for this paper.

The code Version 1.0.

If you are interested in it but have doubts about the content, please feel free to contact me.

Communication: [yangziyi091100@163.com](mailto:yangziyi091100@163.com)

## I. Abstract

The widespread applications in microarray technology have produced the vast quantity of publicly available gene expression datasets. However, analysis of gene expression data using biostatistics and machine learning approaches is a challenging task due to (1) high noise; (2) small sample size with high dimensionality; (3) batch effects and (4) low reproducibility of significant biomarkers. These issues reveal the complexity of gene expression data, thus significantly obstructing microarray technology in clinical applications. The integrative analysis offers an opportunity to address these issues and provides a more comprehensive understanding of the biological systems, but current methods have several limitations. This work leverages state of the art machine learning development for multiple gene expression datasets integration, classification and identification of significant biomarkers. We design a novel integrative framework, MVIAm - Multi-View based Integrative Analysis of microarray data for identifying biomarkers. It applies multiple cross-platform normalization methods to aggregate multiple datasets into a multi-view dataset and utilizes a robust learning mechanism Multi-View Self-Paced Learning (MVSPL) for gene selection in cancer classification problems. We demonstrate the capabilities of MVIAm using simulated data and studies of breast cancer and lung cancer, it can be applied flexibly and is an effective tool for facing the four challenges of gene expression data analysis. Our proposed model makes microarray integrative analysis more systematic and expands its range of applications.

## II. Introduce about code

### i. The repository can be divided into two parts:

1. Codes for simulated experiments.
2. Codes for real gene expression data experiments.

### ii . The compared methods:

We demonstrate the performance of the proposed MVSPL in simulation and real microarray experiments.

Four methods are compared with the MVSPL method: Sparse logistic regression with the Lasso penalty (L_1) [1], Sparse logistic regression with the Elastic Net penalty (L\_EN) [2], Self-Paced Learning (SPL) [3], and Ensemble-based Elastic Net (Ensemble\_EN) [4].

When MVIAm generates single-view data, it degenerates into traditional ``early stage'' data integration, and data analysis can be performed by L\_1, L\_EN and SPL. Ensemble\_EN constructs a prediction model on each view of data before combing the model predictions and obtains the final prediction result based on Equation (1).
$$
y_k=\mathop{argmin}_{\substack{y_k}}\sum_{j=1}^{m}{L_k\left(y_k,f^{(j)}\left(x_k^{(j)},\beta^{(j)}\right)\right)}         (1)
$$

### iii. The MVSPL model:

The objective function of MVSPL can be expressed as:
$$
\mathop{\min}_{\substack{\beta^{\left(j\right)},v^{\left(j\right)}\in\left[0,1\right]^n,\\j=1,2,\ldots,m}}E(\beta^{\left(j\right)},v^{\left(j\right)};\lambda^{(j)},\gamma^{\left(j\right)},\delta)=\sum_{j=1}^{m}\sum_{i=1}^{n}{v_i^{(j)}L\left(y_i,f^{\left(j\right)}\left(x_i^{\left(j\right)},\beta^{\left(j\right)}\right)\right)}+\sum_{j=1}^{m}{\lambda^{(j)}||\beta^{\left(j\right)}||_{1}} \\
-\sum_{j=1}^{m}\sum_{i=1}^{n}{\gamma^{\left(j\right)}v_i^{(j)}}-\delta\sum_{\substack{1\le k,j\le m,\\k\neq j}}{\left(v^{\left(k\right)}\right)^Tv^{\left(j\right)}}
$$

## III. Codes for simulated experiments

The code for simulated experiments contains three steps:

***Step1:*** Generating simulated data. 

Running "Simulation\_data.m" to generate three independent data sets. The rows of each data set represent samples, and the columns represent features:

```
$ matlab Simulation_data.m
```

***Step2:*** Using multiple cross-platform normalization methods to eliminate the batch effects.

You can use popular cross-platform normalization methods to eliminate the batch effects between different data sets, here we provide two cross-platform normalization methods ComBat [5] and Ber [6]. We use four functions ComBat_p, ComBat_n, ber and ber_bg to eliminate batch effects and generate view1, view2, view3 and view4 of the aggregated multi-view data, respectively. Setting the path and run the code "Simulation_aggregate.R".

```
$ Rscript Simulation_aggregate.R
```

***Step3:*** After integrating multiple homogenous data sets, using training data set to train the modal and obtain the best solution for the model. Need to modify the path before executing the code, e.g. "WhereAmI <- "D:/path".

Running single-view methods: L\_1, L\_EN, and SPL.

```
$ Rscript Simu_LR.R                 ## Using Logistic regression model with L1/Elatic Net penalty to train the model.
$ Rscript Simu_SPL.R                ## Using self-paced learning with L1 penalty to train the model.
```

Running multi-view methods: Ensemble\_EN, and MVSPL:

```
$ Rscript Simu_EnsembleEN.R        ## Using Logistic regression model with Elastic net penalty embeded into ensemble framework to train the model.
$ Rscript Simu_MVSPL.R             ## Using multi-view self-paced learning to train the model.
```



## IV. Codes for real gene expression data experiments

### Step1: Download all the publicly available microarray data sets and pre-processing data for conducting a meta-analysis.

#### (I) Create table of study metadata, with a row for each study. See Lung_main_study_metadata.csv for an example. At a minimum, the following  columns must exist: study, studyDataType, and platformInfo.

• study: Name of the study, which must be unique.
• studyDataType: Indicates how the expression data is stored. There are currently five options for studyDataType: affy_geo, affy_custom, affy_series_matrix, series_matrix, and eset_rds.
• platformInfo: Microarray platform, used for mapping probes to genes. 

#### (II) For each study, download the expression data. The form of the expression data will depend on the studyDataType.

• All the expression data should go in the same folder.
• If studyDataType is affy_geo or affy_custom, the expression data should be a folder containing cel or cel.gz files. The name of the folder should match the name of the study.
• If studyDataType is affy_series_matrix or series_matrix, the expression data should be a file ending in "_series_matrix.txt". The part of the file name prior to "_series_matrix.txt" should match the name of the study.
• If studyDataType is eset_rds, the expression data should be an RDS file containing a Bioconductor ExpressionSet. The part of the file name prior to ".rds" should match the name of the study.

Create fold named "data" in the current path, and download the dataset from GEO database. such as affy_series_matrix and affy_geo.

#### (III) Download all files and install all packages necessary for mapping probes to genes and for conducting a meta-analysis.

```
$ Rscript meta_function_Install.R
```

**For each study with studyDataType of affy_geo or affy_custom, download the corresponding custom CDF package.**

• Go to http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp.
• Click on the link for the latest version of "ENTREZG".
• Download the appropriate R source package(s), which should end in "cdf_nn.0.0.tar.gz", where nn refers to the version number.
• Install the custom CDF packages.
• Open R and set the working directory to the folder containing the 
downloaded custom CDF packages. Execute the following command for each 
custom CDF.
• > install.packages('file name of custom CDF', repos=NULL, type='source')
• Make sure that the part of the file name prior to "nn.0.0.tar.gz" matches the information in the platformInfo column.

**For each study with studyDataType of affy_series_matrix, download the corresponding custom CDF zip file.**

• Go to http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp.
• Click on the link for the latest version of "ENTREZG".
• Download the appropriate zip file(s) (under column "CDF Seq, Map, Desc").
• Unzip each zip file.
• Move the file that ends in "_mapping.txt" to the folder that contains the downloaded expression data.
• Make sure that the part of the file name that occurs prior to 
"_mapping.txt" matches the information in the platformInfo column.

#### (IV) Create table of sample metadata, with a row for each sample. See Lung_sample_metadata.csv for an example.

• At a minimum, the following columns must exist: study, sample, and a column for outcome or class.
• The values in the study column should match the values in the study column in the table of study metadata.
• For studies with a studyDataType of affy_geo or affy_custom, the sample names should be the names of the cel or cel.gz files (excluding the file extension).
• For studies with a studyDataType of affy_series_matrix or series_matrix, the sample names should be the GSM identifiers.
• For studies with a studyDataType of eset_rds, the samples names should be the colnames of the corresponding ExpressionSet.
• The name of each sample must be unique across all the samples.

#### (V) Import all the functions for pre-processing all the datasets.

```
source('meta_function.R')
source('meta_function_Pamr.R')
```

#### (VI) Merge expression data and perform cross-study normalization of training datasets, save the expression data and gene name.

```
$ Rscript preprocess_data.R
```



In this paper, we curated data from eight publicly available microarray studies: four breast cancer datasets (same platform) and four lung cancer datasets (disparate platform), as listed in Table 1 and Table 2. All these publicly available cancer gene expression datasets can be download from GEO (https://www.ncbi.nlm.nih.gov/geo/).

Table 1. Four publicly available breast cancer gene expression datasets used in the real data experiments.

| Dataset  | No. of Probes | Classes (Class1 / Class2) | No. of samples (Class1 / Class2) | Affymetrix Platform |
| :------: | :-----------: | :-----------------------: | :------------------------------: | :-----------------: |
| GSE1561  |     22215     |         -ve / +ve         |           49 (22 / 27)           |      HG-U133A       |
| GSE6532  |     22283     |         -ve / +ve         |          125 (40 / 85)           |      HG-U133A       |
| GSE20437 |     22283     |         -ve / +ve         |            18 (9 / 9)            |      HG-U133A       |
| GSE22093 |     22283     |         -ve / +ve         |           82 (41 / 41)           |      HG-U133A       |

Table 2. Four publicly available lung cancer gene expression datasets used in the real data experiments.

| Dataset  | No. of Probes | Classes (Class1 / Class2) | No. of samples (Class1 / Class2) | Affymetrix Platform |
| :------: | :-----------: | :-----------------------: | :------------------------------: | :-----------------: |
| GSE10072 |     22284     |      Normal / Tumor       |          107 (49 / 58)           |        U133A        |
| GSE19188 |     54675     |      Normal / Tumor       |          179 (88 / 91)           |    U133 Plus 2.0    |
| GSE19804 |     54676     |      Normal / Tumor       |          120 (60 / 60)           |    U133 Plus 2.0    |
| GSE43346 |     22283     |      Normal / Tumor       |           65 (42 / 23)           |        U133A        |



### Step2:  After integrating multiple gene expression data sets, using training data set to train the modal and obtain the best solution for the model. 

Need to modify the path before executing the code, e.g. "WhereAmI <- "D:/path".

Running single-view methods: L\_1, L\_EN, and SPL.

```
$ Rscript LR.R                 ## Using Logistic regression model with L1/Elatic Net penalty to train the model.
$ Rscript SPL.R                ## Using self-paced learning with L1 penalty to train the model.
```

Running multi-view methods: Ensemble\_EN, and MVSPL:

```
$ Rscript EnsembleEN.R        ## Using Logistic regression model with Elastic net penalty embeded into ensemble framework to train the model.
$ Rscript MVSPL.R             ## Using multi-view self-paced learning to train the model.
```

The result include the predicted labels of training samples, the predicted labels of test samples, and the selected genes by the model (selected gene and its coefficient). For example:

```
$valpredmatrix
          1  
GSM26804  "1"
GSM26867  "1"
GSM26868  "1"
GSM26869  "1"
GSM26870  "1"
GSM26871  "0"
GSM26875  "1"
GSM26876  "1"
GSM26877  "1"
GSM26878  "0"
GSM26879  "1"
GSM26883  "0"
GSM26884  "0"
...
$evlpredmatrix
          1  
GSM26872  "1"
GSM26873  "1"
GSM26874  "1"
GSM26880  "0"
GSM26881  "1"
GSM26882  "0"
GSM26886  "1"
GSM26889  "1"
GSM26890  "1"
GSM26893  "0"
GSM26894  "1"
GSM26896  "1"
GSM26898  "0"
GSM26900  "1"
...
coef.gene   coef.value          
 [1,] "KPNA5"     "-2.56373113493573" 
 [2,] "RNASE2"    "-1.46375323666761" 
 [3,] "BLCAP"     "1.02292385970561"  
 [4,] "ATXN7L3B"  "1.0087944164481"   
 [5,] "CCNC"      "-0.839655008878292"
 [6,] "ARID1A"    "-0.687994402204038"
 [7,] "TNPO1"     "0.541775628966189" 
 [8,] "SPOP"      "0.541041907153896" 
 [9,] "UBA5"      "-0.530869845792436"
[10,] "CACFD1"    "0.522049776533395" 
...
```

# Reference

[1] Tibshirani, R. Regression shrinkage and selection via the lasso. J. Royal Stat. Soc. Ser. B (Methodological) 267-288 (1996).

[2] Zou, H. & Hastie, T. Regularization and variable selection via the elastic net. J. Royal Stat. Soc. Ser. B (Statistical Methodol.)
67, 301-320 (2005).

[3]  Kumar, M. P., Turki, H., Preston, D. & Koller, D. Learning specific-class segmentation from diverse data. In Computer
Vision (ICCV), 2011 IEEE International Conference on, 1800-1807 (IEEE, 2011).

[4]   Günther, O. P. et al. A computational pipeline for the development of multi-marker bio-signature panels and ensemble classifiers. BMC bioinformatics 13, 326 (2012).

[5] Johnson, W. E., Li, C. & Rabinovic, A. Adjusting batch effects in microarray expression data using empirical bayes
methods. Biostatistics 8, 118-127 (2007).

[6] Giordan, M. A two-stage procedure for the removal of batch effects in microarray studies. Stat. Biosci. 6, 73-84 (2014).
