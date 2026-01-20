-------
Read Me
-------

In the following, the word "paper" refers to the publication 'Canonical host-pathogen tradeoffs subverted by mutations with dual benefits' by Beardmore, Meyer, Gudelj, Hewlett and Pena-Miller, accepted for publication by the American Naturalist in 2022.

This paper uses Matlab extensively. Matlab 2020 or later is recommended, earlier versions may still work but we have found they can generate different results due to changes in core algorithms.

This file contains instructions on the data files used in this paper and the Matlab scripts (m-files) and functions used to analyse them. More information is provided in the paper supplementary. Raw experimental or imaging data are held in CSV, MAT and XLS formats. CSV contain the same data as XLS file as comma-separated values and, note, NO data analysis takes place within any XLS files, they simply contain data. MAT files contain data that is, or has been, analysed using m-files and saved to disk.

Our scripts use a 3rd-party code `Functions for the rectangular assignment problem' by Markus Buehren. This must be installed into Matlab for our scripts to work correctly:

https://uk.mathworks.com/matlabcentral/fileexchange/6543-functions-for-the-rectangular-assignment-problem

The printing suite export_fig is not strictly necessary (it can be replaced in principle by the Matlab command 'exportgraphics') but it is used in our scripts:

https://uk.mathworks.com/matlabcentral/fileexchange/23629-export_fig/

Data and scripts have the following directory structure with . as the root:

./figures/

./figures/unused	(PDFs generated during analyses but not used in the paper) 

./src/			(NB: this must be in the Matlab path for scripts to function)

./src/spheretest/	(Illustrates the use of a matching algorithm described in the supplement, does not contain any data analysis scripts)

./data/			(Contains .MAT files saved following analyses)

./dataRepo/
./dataRepo/derivedData	(Contains .CSV and .XLS files saved following an analysis)
./dataRepo/rawData	(Contains .CSV and .XLS files prior to analysis)


The following files can be run with a single click of the Matlab "green arrow" in the EDITOR menu and are used to perform the core data production and analysis functions. They can also be run by typing their name into the command window:

1) ./validVariables.m - changing this file makes data quality control either more or less severe. Do not change it to keep the same QC as used in the paper. See the paper supplement for information on how to adjust this file, if required.

2) ./dataTasks.m - running this will produce all the "derived data" files needed for the paper, starting from raw biological data produced in the lab.

3) ./figureTasks.m - running this will produce all the figures in the paper, saving them with appropriate filenames to reflect their figure number.

If these scripts are run in the order 1-2-3, then PDF files will subsequently be produced. Scripts make use of the Matlab Parallel Toolbox, if it is available.

-------
Beware:
-------
Some of the figure window sizes have been chosen to be large. When these scripts are run on a display that cannot house these sizes, tests show this is likely to lead to a change in window size which will impinge on the quality of the figures that are subsequently produced. In this case, some editing of the scripts will be necessary to produce the figures as they appear in the paper.

A file not used in any Matlab scripts:
--------------------------------------

./dataRepo/rawData/relativeFitnesses/ColonyCountsForRelFitness.csv

This file is used in the computation of relative fitnesses in Figure 7C, these were not done in Matlab but by hand. The raw colony counts were used to calculate the values given in the files Fitnesses.xlsx (also Fitnesses.csv) which are then displayed in Figure 7C.

--

Zenodo Version 1.0 of these files was used during submission (https://zenodo.org/record/6883319#.YtqOyS1Q3aU), the most recent updates following publication will be at least Version 2.0.

