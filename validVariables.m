function [validMaxGrowthRate , validMaxYield, ...
    validMaxCoV , keeplogisticRsquaredValue , RateYieldOutlierSigma , keepRsquaredValueInCostPlot ] = validVariables()

	validMaxGrowthRate = 4; % the maximal per hour growth rate, anything above this is unphysical
	validMaxYield = 1e-2; % the maximal yield, anything above this is unphysical
	validMaxCoV = 0.2; % the maximal coefficient of variation allowed for derived data
	keeplogisticRsquaredValue = 0.7; % reject logistic fits to longitudinal OD data with $R^{2}$ below this
	RateYieldOutlierSigma = 3; % reject yield values that are 3$\sigma$ from the mean as outliers
	keepRsquaredValueInCostPlot = 0.5; % when analysing resistance costs, reject phenotypes with $R^{2}$ below this
	
end