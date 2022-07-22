%% this does that:

clc
clear all
close all

load('./data/RKandLagData.mat');
load('./data/allData.mat')

%%

plotRYTOResistanceOutliers(data,allRKLagData);