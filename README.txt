For more information see the official BRAVO site:

https://sites.google.com/site/bravotoolbox/

BRAVO is a Matlab toolbox for performing simple and nested regression
analysis on voxelwise observations in MRI data.  The key feature of
BRAVO is the use of "bootstrap" or permutation statistics to estimate
statistical significance.  While this toolbox designed with structural
MRI data in mind (e.g., fractional ani`sotropy, voxel based
morphometry), these routines can be applied to any voxelwise measure.
BRAVO is constructed to be a stand-alone program with minimal
dependencies except for the Matlab statistics toolbox and NifTI data
format.

BRAVO can perform multiple types of analyses:

Correlation: After regressing out covariate factors, BRAVO will
perform either a parametric (Pearson's) or non-parametric (Spearman's)
correlation between an external variable and the voxelwise measure.

Regression: An ordinary least squares regression on any set of
dependent variables and the voxelwise measure.  Allows for including
covariate factors to control for in the regression.

Mediation: Uses nested regression models to estimate causal mediating
pathways between a dependent variable, one or multiple potential
mediators, and the voxelwise measure, while controlling for nuisance
variables. The new release (Oct. 2014, BRAVO 2.0) includes options for
moderated mediation (on a-pathways) and 2-step mediation.

Performance

BRAVO was benchmarked and developed in Matlab R2013a on an Ubuntu
14.04 system with an Intel(R) Xeon(R) processor (4 cores, 3.10GHz)
system and 8Gb of RAM.  Minimum recommended system specifications are
dual core processor (at least 2Ghz) and 4Gb RAM (depending on dataset
size).  All code was tested using DTI and VBM data with a 91x109x91
matrix size, 2mm voxels and up to 141 participants with analysis times
ranging from 2-4hrs.  Larger datasets may significantly lengthen
processing time and system memory requirements.

Developers

Nearly all BRAVO source code was written by Timothy Verstynen, in
collaboration with Kirk Erickson, Peter Gianaros, and Andrea Weinstein
at the University of Pittsburgh.  BRAVO also uses some open source
software developed by Alexandros Leontitsis (Spearmans correlation
code) and Jimmy Shen (NIFTI loader functions).
