function lagPointData = getFirstNonlinearPoint2(someData,pTest)

    if nargin < 2
        pTest = 0.01;
        slopeTest = 1e-4;
    end

    someData = someData(:);
    someData = someData/max(someData);
    
    [n,~] = size(someData);
    someData = myFFTSmooth(someData,20);
    someData = mySmooth(someData,10);
    
    mdl = fitlm(1:n,someData);
    slope = mdl.Coefficients.Estimate(2);
    pValue = mdl.Coefficients.pValue(2);

    if slope > 0 && pValue < pTest
        notYetNonlinear = ones(n,1);
        lagPointData.j = NaN;

        for j = 3:n
            mdl = fitlm(1:j,someData(1:j));
            slope = mdl.Coefficients.Estimate(2);
            pValue = mdl.Coefficients.pValue(2);
            notYetNonlinear(j) = ~((slope > slopeTest) && (pValue < pTest));
        end
        
        if ~all(notYetNonlinear)
            lagPointData.j = find(notYetNonlinear == 1,1,'last');
            lagPointData.pValue = pValue;
            lagPointData.slope = slope;
        end
    else
        lagPointData.j = NaN;        
    end
end