bipartite
=========

repository for the bipartite R-package for network analysis


Hi,

this is the bipartite repository on github. Currently this package is maintained by Carsten Dormann, that's why the repository is hosted under this department.

It falls onto the maintainer to submit new versions to CRAN for release. All non-CRAN-versions are in development, which you can install using the following code (only if you have the tool-chain installed, see R instructions):

    library(devtools)
    install_github("biometry/bipartite/bipartite")

Since I do not use github everyday, please also send me an email if you add/change something (carsten.dormann@biom.uni-freiburg.de). Thanks!

If you want to add/change something, please add as many comments as possible into the R-code (to make other understand what you are doing), edit the .Rd help file accordingly (here you can use LaTeX-style comments using the %), and finally edit the file bipartite-package.Rd in the versioning section, making users aware of what you have added/changed.

For maximal confusion, the actual package content is in the "bipartite" subfolder, while other material (such as test files, wishlists, comments, the current working package compiled etc.) are directly in this folder.

Carsten

