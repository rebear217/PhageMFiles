function [mPD,sPD,D,sePD] = meanPairwiseDistance(v)

    v = v(:);
    O = ones(size(v));

    D = v*O' - O*v';    
    D = (D(:));
    absD = abs(D);
    
    mPD = mean(absD);
    sPD = std(absD);
    sePD = sPD/sqrt(length(v)^2 - 1 - length(v));

end