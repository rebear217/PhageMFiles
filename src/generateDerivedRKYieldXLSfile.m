%%export RK data to Excel

clc
clear all
close all

load('./data/RKandLagData.mat');

%%

for r = 1:46
  
    thisSheet = {};
    thisSheet{1,1} = ['bacterium label ',allRKLagData{r}.bacterium];
    
    thisSheet{2,1} = 'R (growth rate per h)';
    thisSheet{2,2} = allRKLagData{r}.Rdata;
    thisSheet{3,1} = 'R algorithmic standard error estimate';
    thisSheet{3,2} = allRKLagData{r}.Rsedata;

    thisSheet{4,1} = 'K (population size OD600nm)';
    thisSheet{4,2} = allRKLagData{r}.Kdata;
    thisSheet{5,1} = 'K algorithmic standard error estimate';
    thisSheet{5,2} = allRKLagData{r}.Ksedata;
    
    thisSheet{6,1} = 'Yield (OD600nm per ug/mL malto3)';
    thisSheet{6,2} = allRKLagData{r}.Yielddata;
    thisSheet{7,1} = 'Yield algorithmic standard error estimate';
    thisSheet{7,2} = allRKLagData{r}.Yieldsedata;
    
    thisSheet{8,1} = 'Lag estimate 1 (h)';
    thisSheet{8,2} = allRKLagData{r}.Lagdata;

	thisSheet{9,1} = 'Lag estimate 2 (h)';
    thisSheet{9,2} = allRKLagData{r}.Lagdata2;

	thisSheet{10,1} = 'malto3 concentration used for this datapoint (ug/mL)';
    thisSheet{10,2} = allRKLagData{r}.allSugars;

	thisSheet{11,1} = 'Rsquareds for goodness of fit in [0,1]';
    thisSheet{11,2} = allRKLagData{r}.Rsquareddata;
    
    writecell(thisSheet,'./data/derivedRKYieldData.xlsx','Sheet',allRKLagData{r}.bacterium);
end

