function plotMutantDistancesOnLamB(mutantLabel,data,pNorm,blobSize)

    if nargin < 3
        pNorm = 2;
    end
    if nargin < 4
        blobSize = 60;
    end
    
	doNotPlot = 1;

    x = getLocationLabels('wt',data.bacteria);
    y = getLocationLabels(mutantLabel,data.bacteria);
    
    pWt = data.carbonData{x};
    
    %we seem to need to limit the kind of proteins we work with...
    %doSomething = (length(pWt)==length(dataSet.carbonData{y}));    
    doSomething = (length(data.carbonData{y}) > 300);
    
    if doSomething && ~(y == x)
        [pNew,q,code,worstIndexList,bestKeepies,matchedDistance,totalDistance] = rectangularMatchingLabels('wt',data.bacteria{y},data,doNotPlot);
        
        if not(length(pWt) == length(pNew)) || not(length(pWt) == length(q))
            if length(q) == length(pWt)
                %so q is the (shorter) wt and pNew is an insertion mutation
                pNew = pNew(bestKeepies,:);
            else
                %so q is the (shorter) deletion mutant and pNew is the WT
                pNew = pNew(bestKeepies,:);
            end
            %in either case, do the same: reduce the long protein to the length of the shorter one
        end
        
        if isinf(pNorm)
            dVector = max(abs((pNew - q)'));
        else
            dVector = sum(abs((pNew - q)').^pNorm).^(1/pNorm);
        end
    elseif ~(y == x)
        error(['You have asked for a some strain I can''t find: ',mutantLabel])
    elseif (y == x)
        disp('You have asked for a mutant that is the WT')
        dVector = 0;
    end
    
    plotAProteinColourVector(pWt,dVector,blobSize);
    view(52,10)
end