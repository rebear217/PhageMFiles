function [FitData,P,Pm] = logisticFitToODfromLabel(sugar,data,label,plotFlag)
    if nargin < 4
        plotFlag = 0;
    end
    
    I = getLocationLabels(label,data.bacteria);
    if isempty(I)
        error('Cannot find specified bacterial label');
    end
    
    T = 24*60;
    P = [];
    Pm = [];
    
    F = find(data.sugars{I} == sugar);
    for j = 1:length(F)
        [fitData{j},P(j),Pm(j)] = logisticFitToOD(T,data.OD{I}(:,F(j)),plotFlag);
        %[fitData{j},P(j)] = logisticFitToODNewLag(T,data.OD{I}(:,F(j)),plotFlag);
        
        FitData.Rsquared(j) = fitData{j}.Rsquared;
        FitData.growthRate(j) = fitData{j}.logisticCoefficients(3);
        FitData.growthRateSE(j) = fitData{j}.logisticCoefficientsSE(3);
        FitData.K(j) = fitData{j}.K;
        FitData.Kse(j) = fitData{j}.Kse;
        FitData.R(j) = fitData{j}.R;
        FitData.Rse(j) = fitData{j}.Rse;
        FitData.doublingTimeInhours(j) = fitData{j}.doublingTimeInhours;
        FitData.doublingTimeInhoursSE(j) = fitData{j}.doublingTimeInhoursSE;
        FitData.blank(j) = fitData{j}.blank;
        FitData.blankse(j) = fitData{j}.blankse;
        FitData.innoculum(j) = fitData{j}.innoculum;
        FitData.lagTimeMeasure1(j) = fitData{j}.lagTimeMeasure1;
        FitData.lagTimeMeasure2(j) = fitData{j}.lagTimeMeasure2;
        FitData.lagTimeMeasure3(j) = fitData{j}.lagTimeMeasure3;
        FitData.lagTimeMeasure4(j) = fitData{j}.lagTimeMeasure4;
        FitData.sugar = sugar;
        FitData.bacterium = label;
    end
end