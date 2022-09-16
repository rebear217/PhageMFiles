
addpath('./src/')

%%

clear all
close all
clc

% generate raw Ecoli OD data from XLS files:
% It also curates protein data from a MAT file:

importAllTheData;

%%

close all
clear all
clc

disp('Please wait ...')
load('./data/allData.mat')

%%

% Generate derived growth rate and yield data from OD:
generateAllRKandLagData(data)

%% Datafile generation
% these make all necessary datafiles:

generateDerivedRKYieldXLSfile;
disp('done RK xls file ...')

%%

generateDerivedRRmaxKmXLSfile;
disp('done RRK xls file ...')

%%

disp('Please wait, making CSV files ... ')
csvTasks()

%%

disp('All done!')