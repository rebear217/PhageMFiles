function [alphaCarbonHistogramData,alphaList] = plotAllLamBChanges(data,minDistance,SNPflag)

    if nargin < 3
        SNPflag = 1;
    end

    grey = [1 1 1]*0.8;

    if nargin < 2
        minDistance = 0.01;
    end
    %fraction of total differences that are to be highlighted

    figure(1);
    set(1,'pos',[145         436        1146         910])
    view(54,2); 

    x=getLocationLabels('wt',data.bacteria);
    pWt = data.carbonData{x};
    plotAProtein2(pWt,grey);
    
    phi = data.infection.phi;
    S = sum(phi);
    s = 1 - S/max(S);

    alphaCarbonHistogramData = zeros(length(pWt),1);
    alphaList = [];
    
    for y = 1:length(data.bacteria)
        Y = getLocationLabels(data.bacteria{y},data.infection.bacteria_names);
        if ~isempty(Y)
            switch SNPflag
                case 2
                    doSomething = (length(pWt)==length(data.carbonData{y}));
                case 3
                    doSomething = not(length(pWt)==length(data.carbonData{y})) && ...
                        (length(data.carbonData{y}) > 400);
                        %this second condition will need changing to just
                        %correspond to mutants that are known to have a functional
                        %lamB
                case 4
                    doSomething = (length(data.carbonData{y}) > 400) && ...
                        (S(Y) == 0);
                otherwise
                    doSomething = true;
            end

            if doSomething && ~(y == x)
                [matchedDistance,totalDistance,pNew,q,code,worstIndexList] = rectangularMatching(pWt,data.carbonData{y},100,true,minDistance);
                alphaCarbonHistogramData(worstIndexList) = alphaCarbonHistogramData(worstIndexList) + 1;
                alphaList = [alphaList ; worstIndexList(:)];
                colour = [s(Y) 0 (1-s(Y))/2];
                plotAProtein1(pWt(worstIndexList,:),colour,s(Y));
                drawnow
            end
        end
    end   

    box off
    axis off    
    
end