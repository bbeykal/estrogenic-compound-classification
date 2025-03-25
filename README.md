# Classification of estrogenic compounds by coupling high content analysis and machine learning algorithms

This repository contains the supplementary R codes used in the analysis presented in the following research article:

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
                  "caret","boot","statip","sm","kableExtra","formatR","reshape","pROC","stats"))
```
