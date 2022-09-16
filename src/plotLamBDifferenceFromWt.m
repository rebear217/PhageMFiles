function plotLamBDifferenceFromWt(mutantLabel,data,minDistance,pNorm)
    blobSize = 16;
    close all
    
    if nargin < 4
        pNorm = 2;
    end
    if nargin < 3
        minDistance = 0.01;
    end
    
	doNotPlot = 0;

    x = getLocationLabels('wt',data.bacteria);
    y = getLocationLabels(mutantLabel,data.bacteria);
    pWt = data.carbonData{x};
    pMut = data.carbonData{y};

    %just to use SNPs of the same protein length as WT:
    %doSomething = (length(pWt)==length(dataSet.carbonData{y}));
    doSomething = (length(data.carbonData{y}) > 300);

    if doSomething
        [bestMatchedDistance,bestTotalDistance,pNew,q,code,worstIndexList,bestKeepies] = this_matchToWT();
        proteinChange = pNew - q;
        %D = norm(pNew(:) - q(:),pNorm);
        if isinf(pNorm)
            dVector = max(abs((proteinChange)'));
        else
            dVector = sum(abs((proteinChange)').^pNorm).^(1/pNorm);
        end

        figure(2)
        set(2,'pos',[742   571   887   774]);
        %proteinChange = proteinChange(:,1).*proteinChange(:,2).*proteinChange(:,3);
        P = plot(proteinChange,'linewidth',1);
        P = [P ; 1];
        hold on

        xlabel('\alpha carbon position')
        ylabel('change in distance from WT')
        axis tight
        ylim([-2 2])
        
        for j = 1:length(worstIndexList)
            P(end) = plot(worstIndexList(j),myMax(proteinChange(worstIndexList(j),:)),'o','markersize',8,'linewidth',4,'color',[0 0 0]);
            text(worstIndexList(j) + 3,0.02 + myMax(proteinChange(worstIndexList(j),:)),num2str(j))
            %P(end) = plot(worstIndexList(j),0,'s','markersize',15,'color',[1 0 0],'markerfacecolor',[1 0 0]);
        end
        %P = [q ; P];
        %legend(P,{['|\Deltax|+|\Deltay|+|\Deltaz| coordinate of ',mutantLabel],'\Deltax','\Deltay','\Deltaz','hotspot locations'});
        legend(P,{['\Deltax coordinate of ',mutantLabel],'\Deltay','\Deltaz','hotspot locations'});
        legend('boxoff')
        
        axes('pos',[0.565 0.06 0.4 0.4])
        plotAProteinColourVector(pWt,dVector,blobSize);
        for j = 1:length(worstIndexList)
            plot3(pWt(worstIndexList(j),1),pWt(worstIndexList(j),2),pWt(worstIndexList(j),3),'o','markersize',18,'linewidth',2,'color',[0 0 0]);
            text(3+pWt(worstIndexList(j),1),3+pWt(worstIndexList(j),2),1+pWt(worstIndexList(j),3),num2str(j),...
                'FontSize',18)

        end
        
        axis off
        view(52,10)
        export_fig(['./figures/',mutantLabel,'toWTmatch2.pdf']);
                
        figure(1)
        export_fig(['./figures/',mutantLabel,'toWTmatch1.pdf']);
        
    end
    
    function [bestMatchedDistance,bestTotalDistance,pNew,q,code,worstIndexList,bestKeepies] = this_matchToWT()

        [bestMatchedDistance,bestTotalDistance,pNew,q,code,worstIndexList,bestKeepies] = rectangularMatching(pWt,pMut,100,doNotPlot,minDistance);
        %[pNew,q,code,worstIndexList,bestKeepies] = rectangularMatchingLabels(dataSet.bacteria{x},dataSet.bacteria{y},dataSet,doNotPlot);
        %pNew is always longer than q (or the same length)

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

    end

end

function mm = myMax(v)

    mx = max(v);
    mi = min(v);
    
    mm = mx;
    if abs(mi) > abs(mx)
        mm = mi;
    end

end