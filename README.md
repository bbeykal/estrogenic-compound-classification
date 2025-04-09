# Classification of estrogenic compounds by coupling high content analysis and machine learning algorithms

This repository contains the supplementary R Codes used in the analysis presented in the following research article:

> Mukherjee, R., Beykal, B., Szafran, A.T., Onel, M., Stossi, F., Mancini, M.G., Lloyd, D., Wright, F.A., Zhou, L., Mancini, M.A. and Pistikopoulos, E.N., 2020. Classification of estrogenic compounds by coupling high content analysis and machine learning algorithms. PLoS computational biology, 16(9), e1008191. DOI: [10.1371/journal.pcbi.1008191](https://doi.org/10.1371/journal.pcbi.1008191).

## Getting Started  :rocket:
You will need several libraries to perform this analysis: `dendextend`, `ggplot2`, `ggdendro`, `caret`, `boot`, `randomForest`, `statip`, `sm`, `reshape`, `pROC`, `stats`.

Load required libraries and their dependencies by creating a package list:
```r
package.list = c("dendextend","knitr","ggdendro","ggplot2","randomForest",
                  "caret","boot","statip","sm","kableExtra","formatR","reshape","pROC","stats")

if (length(setdiff(package.list, rownames(installed.packages()))) > 0) {
   install.packages(setdiff(package.list, rownames(installed.packages())), dependencies=T,
                        repos = "https://cran.rstudio.com")}
app <- lapply(package.list, require, character.only = TRUE, quietly=TRUE)
```

Alternatively, you can download required libraries one by one, for example:

```r
install.packages("dendextend")
```

or using the vector form:

```r
install.pacakges(c("dendextend","knitr","ggdendro","ggplot2","randomForest",
                    "caret","boot","statip","sm", "kableExtra","formatR","reshape",
                    "pROC","stats"))
```

## Accessing the Datasets & Setting the Directories :open_file_folder:
Datasets can be downloaded from the Supplementary Material of the original paper or the [Texas A&M repository](https://paroc.tamu.edu/Software/Mukherjee_etAl_2020_data.zip). Once unzipped, make sure the R Codes are housed under the same folder with replicate datasets housed in a different subfolder. I suggest renaming the replicates folder from `Reps4-21` to `EPA_Reps_CSV` to not to run into directory issues in the R Codes.

## Running the R Code  :computer:

### Check & Set the Working Directory in R

You can check the working directory by simply using this function in the R console:
```r
getwd()
```
You can set the working directory by simply using this function in the R console:
```r
setwd('path/to/your/desired/working_directory')
```
If you are using RStudio, you can also quickly set the working directory to your R Code by right clicking on the name of the R Script tab, and then clicking on the `Set Working Directory` option.

### Run the Analysis

Source `Project4_Analyze_v5.r` in the R Console. This will automatically do the entire analysis, generate plots, and save the plots to your working directory. It involves the following key steps:

1. Data Cleaning
2. Outlier Detection via Hierarchical Clustering
3. Data Normalization
4. Feature Selection
5. Classification Model Training
6. Data Visualization with PCA
7. Wilcoxon Rank Sum Test
8. Logistic Regression Model Training & Validation
9. Random Forest Classifier Training & Validation
10. Visualizing Model Performance Metrics
    - Boxplots
    - Feature Density Plots
   
## Final Notes :pencil:
The html RMarkdown of the original analysis can be found in this repository. This html file documents the session information with the specific versions of the R packages used in the analysis. 

Common errors faced running the codes: the packages or directories not found. Make sure the right pacakages are installed and the directories have the correct names.
