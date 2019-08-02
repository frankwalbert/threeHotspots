This README describes the structure of the R code for analyzing plate reader data.

A) There are two source files, which define functions used in the indivudal plate reader analyses. These files are:

1. Maggie_MidLogRatio3_asFunction_withGrowthcurver.R

This was initially written by Maggie Kliebhan (hence the name) and later expanded by FWA. The functions in here take individual output files exported from the plate reader and extract descriptive variables from each well (OD & fluoresence values at a given set of timepoints, growth rates, capacity). For the latter two, we use the growthcurver package. For the former, it calculates the same statistics but at various time points, which are:
- "Midlog": the point at which the OD is the average of the min and the max, which is crudely similar to mid-log. In practice, we find this to be often too early in the plate reader run, such that fluorescent values are too low and not yet reliable
- "UserDef": the user can specify a target OD value at which to compute expression. Defaults to OD=0.25. In practice, this is often very close to the inflection point we use in the paper, but results in somewhat less clean results
- "Sat": "saturation", the highest OD point
- "TMid": the inflection point as reported by growthcurver as t_mid. We use this measure in the paper.

These functions require the plate reader output data to be in the format our plate reader produces. See the many provided export files as examples. Specifically, it requires the OD600 data (96 wells, 96 time points) to be located above the fluorescent data. Importantly, the code is NOT smart enough to deal with multiple fluorescent colors read on the same plate. It will only use the first data table under the OD table and assume that this is "GFP". It is up to the user to modify function calls and manage file output accordingly. The code for indiviual experiments has many examples of this.

These functions also require a second file for each exported plate that annotates the wells on that exported plate. It uses this information to label wells with genotypes. It is also needed to define blank wells for background correction. See the various files for examples. The columns in these files are:
- "Wells": self-explanatory
- "Blank": if the well is a media blank, this needs to say "blank". otherwise, leave the cell empty.
- "Desc_1": we usually use this to denote strain background (BY or RM)
- "Desc_2": we usually use this to denote which fluorophore is carried by the strain. for untagged strains, we've defined this to say "noGFP"
- "Desc_3": the various genotypes of interest during fine mapping or variant comparisons (wild type, edited to a certain chimera or variant, etc). Wildtypes should be called "wt". This gets used during fine mapping for comparisons of multiple edited strains to the same wildtype.
- "Desc_4": clone ID, i.e. identify different colonies picked when creating a given genotype. This is important because it is used to group clones within and across plates, which matters for statistical modeling
- "Desc_5": not really used; we sometimes denote strain IDs here

These functions create several output files for each plate:
- an "allGRowthCurvesCombined" pdf. This is an attempt to show all the raw data from this plate. It sometimes gives interpretable plots and sometimes struggles if there are too many different genotypes on the given plate. Can be useful for detecting samples that grew poorly or show other obvious artifacts
- an "output" pdf: one page per well, showing OD and fluoresence, raw and log transformed, with the four schemata for selecting time points highlighted. Useful for seeing if chosen time points make sense, and if there are other problems
- a growthcurver output pdf: shows the growthcurver fits for each well.
- output.txt: this table has the various computed values for each well. These values will be used in downstream statistical analyses.


2. functionsAcrossPlates.R

This file has plotting and statistical testing functions used in fine mapping. The figures in the paper were created from these plots, with some editing in Illustrator.

Note that many of the resulting plot pdf files have several pages, with one measure per page. E.g. the fluorescent measures at the inflection point are sometimes shown on page 4.


B) Each experiment has its own analysis script that calls the various functions above and filters, selects, and transforms data as needed. The structure of these files is similar across experiments, but they do differ in file names, factor orders, and other specifics to each experiment

Some of these files also contain modified versions of the statistical tests or plots, when the given experiments didn't neatly fold into the fine-mapping workflow for which the functions above were defined (e.g. RGT2 GxE).

For convenience, we provide the data and code files in one folder for each experiment.
