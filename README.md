# cfw

This github repository contains code implementing QTL mapping of
physiological, behavioral and gene expression phenotypes, and other
data analyses to assess the viability of using Carworth Farms White
(CFW) mice for mapping genes and genetic loci underlying complex
traits relevant to the study of human disease and psychology. This
code accompanies the following publication:

  Parker, C.C., Gopalakrishnan, G., Carbonetto, P., Gonzales, N.M.,
  Leung, E, Park, Y.J., Aryee, E., Davis, J., Blizard, D.A.,
  Ackert-Bicknell, C.L., Lionikas, A., Pritchard, J.K., Palmer, A.A.
  Genome-wide association study of behavioral, physiological and gene
  expression traits in commercially available outbred CFW mice. To
  appear in *Nature Genetics*.

If you use this code for your research, please cite our paper
published in *Nature Genetics*.

*Some additional details to include in this README*: (1) Tested in R
version 3.2.3 on a MacBook Pro (2.53 GHz Intel Core 2 Duo, OS
X 10.11). (2) Explain where to retrieve data necessary to run these
analyses/scripts.

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

### Contact info

For questions and feedback, please contact:

Abraham Palmer<br>
Department of Psychiatry<br>
University of California, San Diego<br>
9500 Gilman Drive<br>
La Jolla, California, USA<br>
aapalmer@ucsd.edu

### Contents

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

+ Verify that map.qtls.gemma.R is working, and clean up code in that
  script.

+ Organize R scripts a bit better.

+ Fix loading of genotype data in calc.pve.R.

### Contributors

The code in this repository was developed by
[Peter Carbonetto](http://www.cs.ubc.ca/spider/pcarbo) and
[Shyam Gopalakrishnan](http://www.google.com).

Other contributors are Clarissa C. Parker, Natalia M. Gonzales, Emily
Leung, Yeonhee J. Park, Emmanuel Aryee, Joe Davis, David A. Blizard,
Cheryl L. Ackert-Bicknell, Arimantas Lionikas, Jonathan K. Pritchard
and Abraham A. Palmer.

### Acknowledgments

This project was funded by NIH R01GM097737 (A.A.P.), NIH T32DA07255
(C.C.P), NIH T32GM07197 (N.M.G.), NIH R01AR056280 (D.A.B.), NIH
R01AR060234 (C.A.B.), the Fellowship from the Human Frontiers Science
Program (P.C.), and the Howard Hughes Medical Institute (J.K.P.). The
authors wish to acknowledge technical assistance from Dana Godfrey,
Sima Lionikaite, Vikte Lionikaite, Ausra S. Lionikiene, and John
Zekos; as well as technical and intellectual input from Drs. Mark
Abney, Justin Borevitz, Karl Broman, Na Cai, Riyan Cheng, Nancy Cox,
Robert Davies, Jonathan Flint, Leo Goodstadt, Paul Grabowski, Bettina
Harr, Ellen Leffler, Richard Mott, Jerome Nicod, John Novembre, Alkes
Price, Matthew Stephens, Daniel Weeks, and Xiang Zhou.

