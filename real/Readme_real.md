**This is code for real gene expression experiments.**

**preprecess_data.R:** using this script to download the gene expression dataset from GEO database, and process each independent gene expression datasets into a multi-view aggregated dataset. (need to prepare two basic csv file: *_main_study*.csv and *_sample_metadata.csv).  

**LR.R:** This code is using logistic regression model with penalty function (Lasso, Elastic net et al.).  

**SPL.R:** This code is based on Self-Paced Learning using logistic regression model with penalty function (Lasso, Elastic net et al.).  

**MVSPL.R:** This code is based on Multi-view Self-Paced Learning using logistic regression model with penalty function (Lasso, Elastic net et al.).
