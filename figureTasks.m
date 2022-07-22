%In the following description of the article text, these two lines mean the same thing:

% Supp Figure x
% Figure Sx

%%

addpath('./src/')

%%

% This makes a clean start with all data loaded:
close all
clear all
clc

load('./data/allData.mat')
load('./data/RKandLagData.mat')

%%

% Figure 2
close all

% Show that PFUs correlate with the imaged phage plaque "blob" area:
analysePFUs;

%%

% Figure 2
close all
% Display the phage-bacteria interaction matrix:
plotInteractionMatrix(data)

%%

% Figure 3

close all
% Display the phage-bacteria interactions regressions:
phiLinearRegressions(data)

%%

% Figure 4
close all
% An example WT protein comparison with the strain 56a....
plotLamBDifferenceFromWt('56a',data);

%%

% Supp Figure S10 (similar to Figure 4):
close all
plotLamBDifferenceFromWt('2b',data);

%%

% Figure 5
close all
allHotspotsLamBMovers(data)

%%

% Figure 5
close all
plotHotSpotClusters;

%%

% Figure 6
close all
plotSeveralLagExamples(data,{'wt','19a'});

%%

% Figure 6
close all
plotAllRYTOs(data,allRKLagData);

%%

% Figure 7
close all
plotResistanceVersusTraitData_FAST(data,4,allRKLagData);
% this currently uses outlier removal at 3 x sigma:

figure(1)
export_fig('./figures/YieldAndRateCorrelation.pdf')
figure(2)
export_fig('./figures/YieldAndRateCorrelationAllData.pdf')

%%

%
% Figure 7
close all
load('./data/relativeFitnessData.mat');
plotRelFitness(lamBstrains,relativeFitnesses,data);

%%

% Figure 8
close all
KmData = plotRYTOResistanceOutliers(data,allRKLagData);
[p,h,stats] = signtest(KmData,1,'Tail','both')

%this produces the file ./data/derivedRelativeTraitData.xlsx

% -------------
% last output:
% -------------
% Unreliable data not used : 13a only has Monod fit adj R^2 = -0.0011583
% Unreliable data not used : 70b only has Monod fit adj R^2 = 0.12543
% Unreliable data not used : 94a only has Monod fit adj R^2 = -0.0012092

%%

% Figure 9
close all
plotGrowthRatesFor({'56a','26a'},data,allRKLagData)

%%

% Supp Figures S1-S5 are made by the above text because these are extensions of the figures in the main text.
% Supp Figure S6 is a pedagogical figure not produced by this script, as is Supp Figure S9 which places the greasy slide locations onto a LamB structure.

%% The synthetic 3d matching "proteins" example in the supplementary: Figure S7

close all
cd('./src/sphereTest')
sphereTest;
cd('../..')

%% Supp Figure network analysis: Figure S8

close all
nestedAnalysis;

%{

============
Last outputs:
============

Modularity:
	Used algorithm:             	        AdaptiveBrim
	N (Number of modules):      	                   4
	Qb (Standard metric):       	              0.0986
	Qr (Ratio of int/ext inter):	             -0.1181
Nestedness NODF:
	NODF (Nestedness value):    	              0.7977
	NODF (Rows value):          	              0.4967
	NODF (Columns value):       	              0.8911

===
AND
===

Modularity
	 Used algorithm:	                        LPBrim
	 Null model:    	       NullModels.EQUIPROBABLE
	 Replicates:    	                          1000
	 Qb value:      	                        0.0963
	     mean:      	                        0.1050
	     std:       	                        0.0101
	     z-score:   	                       -0.8685
	     t-score:   	                      -27.4631
	     percentil: 	                       26.7000
	 Qr value:      	                       -0.0493
	     mean:      	                       -0.1670
	     std:       	                        0.2821
	     z-score:   	                        0.4174
	     t-score:   	                       13.1995
	     percentil: 	                       66.0000
Nestedness
	 Used algorithm:	                NestednessNODF
	 Null model:    	       NullModels.EQUIPROBABLE
	 Replicates:    	                          1000
	 Nestedness value:	                        0.7977
	     mean:      	                        0.4364
	     std:       	                        0.0054
	     z-score:   	                       66.7637
	     t-score:   	                     2111.2538
	     percentil: 	                      100.0000

%}


