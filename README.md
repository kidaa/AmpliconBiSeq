AmpliconBiSeq
========

This is an R package for amplicon bisulfite-sequencing.


Installation
---------
```R 
# dependencies
install.packages( c("data.table","devtools"))
source("http://bioconductor.org/biocLite.R")
biocLite(c("Gviz","QuasR","Rsamtools"))

# install from github
library(devtools)
install_github("methylKit", username = "BIMSBbioinfo",build_vignettes=FALSE)
```


Contact & Questions
-------
open an issue at github for bug reports

Contribute to the development
---------
You can contribute to the AmpiconBiSeq development via github by checking out "development" branch, making your changes and doing a pull request (all of these should be done on the "development" branch NOT on the "master" branch). You can first open an issue to see if we are planning or currently working on the issue.

In addition:
 * Bump up the version in the DESCRIPTION file on the 4th number. For example, the master branch has the version numbering as in "X.Y.Z".  If you make a change to the development branch you should bump up the version in the DESCRIPTION file to "X.Y.Z.1". If there are already changes done in the development just bump up the fourth number. 
 * Add your changes to the NEWS file as well under the correct version. Attribute the changes to yourself, such as "Contributed by X"

License
---------
Artistic License/GPL