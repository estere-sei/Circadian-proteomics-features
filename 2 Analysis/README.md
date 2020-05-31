# Description

Please see the [referenced publication](https://doi.org/10.1101/2020.05.16.099556) for more details.

These analyses are split into "Independent" and "Combined" analyses, where the dataset is split into genotypes either before or after normalisation, respectively. The independent analysis allows the detection of all the proteins including those that would have been excluded because they are only present in one genotype. Naturally, this means that comparisons of absolute abundance between genotypes. The combined analysis allows this because batch correction occurs before splitting up the data. Therefore the combined analysis was used when comparing expression levels between genotypes, as well as when comparing the relationships between abundance deciles and proportion rhythmic.

Apart from the independent vs combined distinction, the workflow thence is standardised. First is the statistical test for rhythmicity - in this case [RAIN](https://journals.sagepub.com/doi/10.1177/0748730414553029?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%3dwww.ncbi.nlm.nih.gov) was used. This is followed by calculations of baseline (i.e. mean abundance) and relative amplitudes. Other analyses and code for plotting various types of graphs are included in subsequent notebooks.

Phosphoproteomics is handled in a similar way to proteomics. However, phosphopeptide levels are adjusted to take into account changes in the abundance of the associated protein. Therefore the normalised protein abundances are required for phosphoproteomics analysis. This was designed because we were interested in changes in phosphorylation, rather than changes in the abundance of the underlying protein. This means that we only include phosphopeptides for which we have proteome data, thus reducing our coverage as an unfortunate side effect.

### Directory structure

The root directory contains:
* "1 Raw files from Perseus" -> proteomics + phosphoproteomics .txt files.
* "2 Analysis" -> split into "Proteome" and "Phosphoproteome".
    
The "2 Analysis" folder contains:
* "Phosphoproteome".
* "Proteome".

Each of "Phosphoproteome" and "Proteome" folders contains:
* "Combined".
* "Independent" - R notebooks all begin with "Single" to denote that it is associated with the "Independent" analysis rather than "Combined".
    
Each of "Combined" and "Independent" folders contains:
* R notebooks - these are numbered in the order in which they are to be run, and each begins with a short description and a summary of input/output files. Numbers are appended with "p" when referring to phosphoproteomics analysis, e.g. "1p Normalising data.Rmd".
* The relevant raw file from Perseus (either proteomics or phosphoproteomics).
* Other relevant input files (e.g. manually annotated lists of interesting proteins).
* Output files (created by running the R notebooks in the specified order).

### Credits
* The half-life analysis was done by Estere Seinkmane. 
* The kinase inference analysis script was written in Python by Dr Tim Stevens.
* The group leader supporting this work was Dr John S. O'Neill, MRC Laboratory of Molecular Biology.
    
### VERSION INFORMATION

R 3.6.3

### LICENSE

This code is released with the [MIT License](LICENSE).
