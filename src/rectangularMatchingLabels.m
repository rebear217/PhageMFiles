function [pNew,q,code,worstIndexList,bestKeepies,matchedDistance,totalDistance] = rectangularMatchingLabels(a,b,data,noplot)

	maxIts = 100;

    if nargin < 4
        noplot = 0;
    end
    
    try
        x=getLocationLabels(a,data.bacteria);
        y=getLocationLabels(b,data.bacteria);
    catch problem
        error('One of your labels does not exist')
    end
    if isempty(x) || isempty(y)
        error('One of your labels does not exist')
    end        
    
    [matchedDistance,totalDistance,pNew,q,code,worstIndexList,bestKeepies] = rectangularMatching(data.carbonData{x},data.carbonData{y},maxIts,noplot);
    if not(noplot)
        disp(['Matched distance is: ',num2str(matchedDistance)])
        disp(['Total distance is: ',num2str(totalDistance)])

        figure(1)
        subplot(1,2,1)
        hold off
        subplot(1,2,2)
        hold off

        subplot(1,2,1);
        %title({['{\bf A comparison of two protein structures: between clones ',a,' and ',b,'}'],...
        %    'green rings denote non-matched \alpha-carbons'})
        text(-36,36,['comparing LamB structures of ',a,' and ',b]);
        subplot(1,2,2);

        set(gcf,'position',[63         216        2100         734])
        export_fig(['./figures/proteinShapes/proteinComparisonFigure',a,'-',b,'.pdf'])
    end
end