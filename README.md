FuseSOM
======================================================

<img src=inst/FuseSOM.png height="200">

A correlation based `Multiview` Self Organizing Map for the characterization of 
cell types in `highly multiplexed in situ imaging cytometry assays` data.

Overview
--------

**FuseSOM** provides a pipeline for the clustering of `highly multiplexed in situ imaging cytometry assays`.
This pipeline uses the `Self Organizing Map` architecture coupled with `Multiview` hierarchical clustering.
We also provide functions for normalisation and estimation of the number of clusters.

Installation
--------

```r
# Install the development version from GitHub:
# install.packages("devtools")
# you will need the yasomi package for the Self Orgnanizing map
install.packages("yasomi", repos="http://R-Forge.R-project.org")
devtools::install_github("https://github.com/ecool50/FuseSOM")
library(FuseSOM)
```

### Submitting an issue or feature request

`FuseSOM` is still under active development. We would greatly appreciate any and 
all feedback related to the package.

* R package related issues should be raised [here](https://github.com/ecool50/FuseSOM/issues).
* For general questions and feedback, please contact us directly via [ewil3501@uni.sydney.edu.au](mailto:ewil3501@uni.sydney.edu.au).


## Author

* **Elijah Willie**
* **Ellis Patrick**  - [@TheEllisPatrick](https://twitter.com/TheEllisPatrick)
