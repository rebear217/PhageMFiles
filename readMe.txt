-----------------------------------------------------------------------------------------------
# Title of Dataset: Canonical host-pathogen tradeoffs subverted by mutations with dual benefits
-----------------------------------------------------------------------------------------------

README author  : r.e.beardmore@exeter.ac.uk
2ndary contact : jrmeyer@ucsd.edu

--------------------------------------------------------------------------------------------
## Brief summary of dataset contents, contextualized in experimental procedures and results:
--------------------------------------------------------------------------------------------
- Throughout, "paper" or "the paper" refers to the publication 'Canonical host-pathogen tradeoffs subverted by mutations with dual benefits' by Beardmore, Meyer, Gudelj, Hewlett and Pena-Miller, accepted for publication by the American Naturalist in 2022. Detailed publication information, like DOI or weblinks to the paper, are not available at the time of writing.

- This paper uses Matlab extensively for its analysis. Matlab 2020 or later is recommended to runs the scripts below, earlier versions may still work in the sense of generating output but we have found they can generate erroneous results due to changes in core algorithms by Mathworks.

- This readme contains instructions concerning data files used in the paper and the Matlab scripts (m-files) and functions used to analyse them. Raw experimental data and imaging (i.e. photographical) data are held in JPG, CSV, MAT and XLS formats. All CSV files mirror the data held in XLS files using comma-separated values format and, we note, NO data analysis takes place within any XLS files, they simply hold data. MAT files contain data that is, or has been, analysed using m-files and subsequently saved.

- Thus there are two types of data referred to below: "raw data" and "derived data" and the directory structure reflects this. The former are data produced by an experimental device or procedure (that can include software processing for determining protein structures) whereas the latter refers to data that have been produced by an additional processing step that uses one of the Matlab scripts mentioned below.

- Our Matlab scripts use the 3rd-party code `Functions for the rectangular assignment problem' by Markus Buehren. This must be installed into Matlab for our scripts to work correctly and it can be accessed here:

	https://uk.mathworks.com/matlabcentral/fileexchange/6543-functions-for-the-rectangular-assignment-problem

- The printing suite export_fig is not strictly necessary as it can be replaced by other Matlab commands, like 'exportgraphics', but it is used for convenience in our scripts and should be installed, it can be accessed here:

	https://uk.mathworks.com/matlabcentral/fileexchange/23629-export_fig/

----
NB1:
----
- The figure window sizes have been chosen to be large in order to obtain PDF outputs that address publication requirements. When these scripts are run on displays that cannot house these sizes, this will lead to changes in window sizes which will impinge on figure quality. In this case, editing of the scripts will be necessary (for instance, changes in fontsize) to produce figures exactly as they appear in the paper.

----
NB2:
----
- Assume you have unpacked the files to a directory structure with "." as the root. The following files can be run with a single click of the Matlab "green arrow" in the EDITOR menu and they perform the core data production and analysis functions as described. They can also be run by typing their name into the command window:

1) ./validVariables.m

- Changing the entries in this file makes the data quality control procedures either more or less conservative, as discussed in the paper. Do not change these values if you wish to keep the same quality control values as used in the paper. Please read the paper's online supplement for information on how to adjust this file appropriately, if required.

2) ./dataTasks.m

- Running this script will produce all the "derived data" files needed for the paper, starting from "raw data" produced in the lab.

3) ./figureTasks.m

- Running this will produce all the figures in the paper, saving them with appropriate filenames to reflect their figure number. This file has sections that can be run one at a time by Matlab in order to produce just one figure.

- If these scripts are run in their numerical order, 1-2-3, then PDF files will subsequently be produced correctly, starting by processing the raw data. These scripts make use of the Matlab Parallel Toolbox to accelerate the algorithms, if it is available.

- Overview of the paper's figure contexts within the data structure:

Figure 1	:none, this uses no data
Figure 2	:A uses "bacteria-phage infectivities" held, for example, in "./dataRepo/derivedData/infectionMatrix" as discussed below
		:B uses "bacterial lawn and phage plaques" held, for example, in "./dataRepo/rawData/phageCalibration" as discussed below
Figure 3	:uses exactly the same data as Figure 2A, but is analysed and presented in a different way
Figure 4	:this illustrates the function of an algorithm using WT LamB protein shape data held in "./dataRepo/rawData/proteinCoords"
		:as discussed below
Figure 5A&B	:these use the algorithm illustrated in Figure 4 applied to the LamB proteins of several bacterial mutants, so again,
		:the relevant raw data is found in "./dataRepo/rawData/proteinCoords" discussed below
Figure 6	:A&C use raw optical density (OD) timeseries data found in, for example, folder "./dataRepo/rawData/OD-RYTO-runs" below
		:B uses yield data that are derived from OD data, found in "./dataRepo/derivedData/RKYield"
Figure 7	:A&B use growth rate and yield data found in "./dataRepo/derivedData/RKYield"
		:C uses relative fitness data found in "./dataRepo/rawData/relativeFitnesses" below
Figure 8AB&C	:all these use relative trait data normalised with respect to the WT strain, found in
		:"./dataRepo/derivedData/WTrelativeTraits" below
Figure 9AB&C	:all these use growth rate data determined from OD trajectory, as found in "./dataRepo/derivedData/RKYield" below

---------------------------------------------
## Description of the Data and file structure
---------------------------------------------

## This is a freeform section for you to describe how the data are structured and how a potential consumer might use them. Be as descriptive as necessary. Keep in mind that users of your data might be new to the field and unfamiliar with common terminology, metrics, etc.
## Describe relationship between data files, missing data codes, other abbreviations used. Be as descriptive as possible.

- Data and scripts have the following directory structure assuming . is the root where you unpacked these files into:

./readMe.txt			:this file
./dataTasks.m			:data curation script discussed above
./figureTasks.m			:figure production script discussed above
./validVariables		:quality control function discussed above (used by dataTasks.m)

./data/				:Contains 5 MAT files saved following data generation scripts applied
				:to spreadsheets i.e. these are all "derived data" & the files are ...
./data/allCarbonData.mat	:contains a variable cell "carbonData" that contains 3d coordinates of all the LamB protein structures
				:bacterial labels are contained within the same variable structure
./data/allData.mat		:contains a variable structure "data" that contains both the protein structures (in data.carbonData),
				:the imaging data that determine bacteria-phage infections (in data.infection) and the optical density
				:(in units of OD) growth data (in cell data.OD) that can be used to determine growth rates and yields
./data/phi-secondTry-large.mat	:contains the same imaging-derived bacteria-phage dataset as data.infection
				:thus - allData.mat packages datasets from other sources into a single location
./data/relativeFitnessData.mat	:contains relative fitnesses of different Ecoli strains when competed with the wild-type strain
./data/RKandLagData.mat		:contains a variable structure of cells "allRKLagData" which contain an analysis of the above OD
				:Ecoli population density data (determined from growth kinetic assays) which are used to derive
				:growth rates (r), maximal observed population densities (K) and various lag time estimates.

./dataRepo/
./dataRepo/derivedData/
./dataRepo/derivedData/infectionMatrix						:"matrixData.xlsx" hold bacteria-phage infectivities
										:in arbitrary imaging units
										:as a matrix with strain labels, also held in CSV format
./dataRepo/derivedData/RKYield							:"derivedRKYieldData.xlsx" has 1 sheet per bacterial strain that
										:contains growth rate, yield and lag (units are specified)
./dataRepo/derivedData/RKYield/derivedRKYieldData-csvSheets			:each of the sheets in "derivedRKYieldData.xlsx" has been
										:exported as a CSV file containing the same information,
										:resulting in 1 CSV per sheet in this directory

./dataRepo/derivedData/WTrelativeTraits						:"derivedRelativeTraitData.xlsx" has 3 sheets, one for each
										:ratio of 3 different phenotypes, thus
										:data are a unitless list of phenotypic ratios (see Figure 8)
./dataRepo/derivedData/WTrelativeTraits/derivedRelativeTraitData-csvSheets	:each of the sheets in derivedRelativeTraitData.xlsx's
										:has been exported as a CSV file to this directory
										:containing the same information, resulting in 1 CSV per sheet

./dataRepo/rawData
./dataRepo/rawData/OD-RYTO-runs		:holds XLS and CSV files (named "XX_RYTO.xlsx" and "XX_RYTO.csv" where XX is the strain label)
					:both containing the same data which are bacterial optical densities (OD) for the  
					:XX labels in each filename (e.g. 2b, 4a etc). Times in minutes are as indicated in the
					:first column, values in the first row are maltotriose concentrations (ug/mL) where "blank" indicates
					:a negative control containing no bacterial inoculate which should result in a constant readout.  
./dataRepo/rawData/phageCalibration	:contains 2 TIF images ("Lab 2022-02-18_09h10m46s.tif" and "Lab 2022-02-18_09h10m46s_inverted.tif")
					:of a bacterial lawn with phage plaques at different inoculation sizes where
					:one file's image is the inversion of the other
					:it also contains an XLS file ("RPM clearing data for Beardmore 2022.xlsx")
					:with 2 sheets which holds the sizes of the plaque clearing determined
					:using imaging algorithms (for examples of these, see standard software ImageJ - https://imagej.net)
					:The 1st sheet of this XLS has 8 replicates (A-H) and 12x 2-fold dilutions that yield
					:different phage plaque sizes.
					:The 2nd sheet labelled "Comparison" has 2 columns which are used to correlate plaque size (in PFU -
					:plaque forming units) with imaging units (Area - which are areas defined by the TIF format in pixel^2)
./dataRepo/rawData/phageCalibration/RPM clearing data for Beardmore 2022-csvSheets
					:This contains 2 CSV files ("Meyer clearing data for Beardmo.csv" and "Comparison.csv"),
					:one for each of the 2 sheets in the above file "RPM clearing data for Beardmore 2022.xlsx"
./dataRepo/rawData/phageImages		:These are bacterial lawn images in JPG format (named "XX.jpg" where XX is a bacterial strain label)
					:with phage plaques and they were used to determine the
					:matrix of infectivities in Figure 2.
./dataRepo/rawData/proteinCoords	:This contains a list of XLS files (named "XX coordinates.xls" where XX is a bacterial strain label)
					:that each contain data on the protein structures of the LamB
					:proteins referred to throughout the paper. Each XLS file is mirrored as a sub-directory that contains
					:CSV files whereby 1 CSV file is used per sheet in the corresponding XLS files. For example...

./dataRepo/rawData/proteinCoords/1a coordinates.xls contains the same information as

./dataRepo/rawData/proteinCoords/1a coordinates-csvSheets/1a coordinates.csv
./dataRepo/rawData/proteinCoords/1a coordinates-csvSheets/AA Pos XYZ.csv
./dataRepo/rawData/proteinCoords/1a coordinates-csvSheets/xyz AA only alpha.csv
./dataRepo/rawData/proteinCoords/1a coordinates-csvSheets/xyz cons AA only alpha C.csv
./dataRepo/rawData/proteinCoords/1a coordinates-csvSheets/xyz cons AA.csv

... where, instead of 1a, other strain labels such as 2b, 4a, etc, are used.

- The only data within file "1a coordinates.xls" (and analogously, 2b coordinates.xls etc) that are used in the paper are the 3d
(x,y,z) coordinates of the alpha carbons that we use as a backbone within the protein structure. This can be extracted from different sheets in each of the "XX coordinates.xls" files (recalling "XX" is the bacterial label). For example, the XLS sheet named "xyz cons AA only alpha C" and the corresponding CSV file ("xyz cons AA only alpha C.csv") have this information which has been imported into Matlab into file "allCarbonData.mat" discussed above.

- Alpha carbon units in all these XLS and CSV files are angstroms. The label "CA" in these files refers to a row of alpha carbon (x,y,z) coordinates.

./dataRepo/rawData/relativeFitnesses	:has "Fitnesses.xls" which contains the relative fitnesses found in Figure 7C. "606" referred to in the
					:XLS is the wild-type known as WT throughout the paper. Files Fitnesses.xlsx and Fitnesses.csv are
					:mirrors. "ColonyCountsForRelFitness.csv" holds the raw colony counts (units cell numbers) used to
					:determine those relative fitnesses. The colours referred to are found on the agar plate used to perform
					:the colony counting and can be ignored here.

./figures/		:all figures mentioned in the paper and its online supplement are held here as PDFs
			:The filename convention used can be explained by the example "Figure6BandCpartial_AndA4" which
			:means a part of this PDF is used in Fig 6B and it is also used in supplement figure 4
			:where the "A4" here denotes "Appendix Figure 4" because "Appendix" has been used as a synonym for the online supplement.
./figures/unused	:PDFs generated during analysis that are not used in the paper but are produced by analysis scripts

./src/			:holds all the Matlab m-files (scripts and functions) used by the paper.
			:NB - this must be in the Matlab path for scripts to function

./src/spheretest	:m-files and figures that, when run, illustrate the use of a geometric
			:matching algorithm described in the supplement
			:spheretest.m is the script that calls other functions and they only use
			:synthetic data, they do not contain any data analysis code
./src/figures		:empty legacy directory that was used to store figures, may still be
			:referred to by some m-files

-----------------------------------------------
NB3: A CSV file not used in any Matlab scripts:
-----------------------------------------------
./dataRepo/rawData/relativeFitnesses/ColonyCountsForRelFitness.csv

- This file is used in the computation of relative fitnesses in Figure 7C, these were not done in Matlab but by hand. The raw colony counts were used to calculate the values given in the files Fitnesses.xlsx (also Fitnesses.csv) which are then displayed in Figure 7C.


-----------------------------
## Sharing/access Information
-----------------------------
- Links to other publicly accessible locations of the data:

- Was data derived from another source?

--
NO
--

If yes, list source(s): N.A.


---------------
- END OF FILE -
---------------
