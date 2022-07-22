function [validMaxGrowthRate , validMaxYield, ...
    validMaxCoV , keeplogisticRsquaredValue , RateYieldOutlierSigma , keepRsquaredValueInCostPlot ] = validVariables()

    %used in the review process as a simple check that noise in data does
    %not alter any conclusions substantively:
    
	%validMaxGrowthRate: the maximal per hour growth rate, anything above this is unphysical
	%validMaxYield: the maximal yield, anything above this is unphysical
	%validMaxCoV: the maximal coefficient of variation allowed for derived data
	%keeplogisticRsquaredValue: reject logistic fits to longitudinal OD data with $R^{2}$ below this
	%RateYieldOutlierSigma: reject yield values that are 3$\sigma$ from the mean as outliers
	%keepRsquaredValueInCostPlot: when analysing resistance costs, reject phenotypes with $R^{2}$ below this
	
    %keep almost all data, accept potentially erroneous data that may be
    %unphysical...this is the set used for figures in the article:
	validMaxGrowthRate = 2.5;
	validMaxYield = 1e-2;
	validMaxCoV = 1;
	keeplogisticRsquaredValue = 0.2;
	RateYieldOutlierSigma = 4;
	keepRsquaredValueInCostPlot = 0.2;
    
end