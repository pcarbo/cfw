# cfw

R code implementing QTL mapping of physiological, behavioral and gene
expression phenotypes in Carworth Farms White (CFW) outbred mice. *A
few more details will go here.*

Some additional details: (1) Tested in R version 3.2.3 on a MacBook
Pro (2.53 GHz Intel Core 2 Duo, OS X 10.11). (2) Explain where to
retrieve data necessary to run these analyses/scripts.

###Citing this resource

*Details go here.*

### License

The *cfw* code repository is free software: you can redistribute it
and/or modify it under the terms of the
[GNU General Public License](http://www.gnu.org/licenses/gpl.html) as
published by the Free Software Foundation, either version 3 of the
License, or any later version.

This program is distributed in the hope that it will be useful, but
**without any warranty**; without even the implied warranty of
**merchantability** or **fitness for a particular purpose**. See file
[LICENSE](LICENSE) for the full text of the license.

###Credits

The *cfw* code was developed by
[Peter Carbonetto](http://www.cs.ubc.ca/spider/pcarbo) and
[Shyam Gopalakrishnan](http://www.google.com).
[Abraham Palmer](http://palmerlab.org) and
[Clarissa Parker](http://www.middlebury.edu/academics/neuro/faculty/node/454157)
also contributed to the development of this software.

### R scripts and modules

*Brief summary of scripts and modules goes here. This will be better
organized later.*

+ calc.pve.R: Estimate the proportion of variance in the phenotype,
after removing linear effects of covariates, that is explained
by the availablegenetic variants.

+ check.pheno.R: This is a small script to check the whether the
observed quantiles for each phenotype, conditioned on different sets
of covariates, match what we would expect under the normal
distribution.

+ examine.binary.covariates.R: A script to show scatterplots of
phenotype vs covariate, and calculate the proportion of variance in a
phenotype of interest that is explained by variance candidate binary
covariates (e.g., testing apparatus).

+ examine.bmd.R: A small script to compare the distribution of
bone-mineral density (BMD) against BMD data from other studies.

+ examine.covariates.R: A script to show scatterplots of phenotype
versus covariate, and calculate the proportion of variance in a
phenotype interest that is explained by various candidate covariates
(e.g., body weight).

+ gen.megamuga.snp.density.plot.R: Script to plot distribution of
MegaMUGA SNPs that are polymorphic in CFW mice across chromosomes
1-19.

+ misc.R: Defines several functions that don't fit anywhere else.

+ plotting.tools.R: Some functions for creating plots to summarize
results.

+ polygenic.R: This file contains functions that implement the
polygenic model for estimating the proportion of variance expained
by available genetic markers. Here is an overview of the functions
defined in this file:

+ qtl.analyses.R: Defines a data structure that provides information
about each QTL analysis.

+ read.data.R: Defines several functions for reading the QTL
experiment data from text files.

### To do items

+ Organize R scripts a bit better.

+ Fix loading of genotype data in calc.pve.R.
