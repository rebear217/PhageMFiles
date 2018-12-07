function lagPointData = getLagTimeIndex(Nj,pTest,slopeTest)
    %default parameters are used to test for non-constant data
    %in terms of slope and p-Value: so lagPointData = getLagTimeIndex(Nj) will execute
    if nargin < 2
        pTest = 0.01;
        slopeTest = 1e-4;
    end
    Nj = Nj(:);
    Nj = Nj/max(Nj); %normalise data to [0,1]
    [n,~] = size(Nj);
    %someData = FilterNoiseFromData(Nj); %un-comment to pre-filter data if need be
    %Fit a linear model to data:
    mdl = fitlm(1:n,Nj);
    slope = mdl.Coefficients.Estimate(2);
    pValue = mdl.Coefficients.pValue(2);
    %now process the regression slope data and so decide
    %if timeseries data are constant, or not:
    if slope > 0 && pValue < pTest
        notYetNonlinear = ones(n,1);
        lagPointData.J = NaN;
        for j = 3:n
            mdl = fitlm(1:j,Nj(1:j));
            slope = mdl.Coefficients.Estimate(2);
            pValue = mdl.Coefficients.pValue(2);
            notYetNonlinear(j) = ~((slope > slopeTest) && (pValue < pTest));
        end        
        if ~all(notYetNonlinear)
            lagPointData.J = find(notYetNonlinear == 1,1,'last');
            lagPointData.pValue = pValue;
            lagPointData.slope = slope;
        end
    else
        lagPointData.J = NaN;        
    end
end
