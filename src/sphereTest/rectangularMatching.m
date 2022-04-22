function [bestMatchedDistance,bestTotalDistance,pNew,q,code,worstIndexList,bestKeepies] = rectangularMatching(p,q,maxIterations,noplot,minDistance)
    if nargin < 5
        minDistance = 0.01;%fraction of total differences that are to be highlighted
    end    
    if nargin < 4
        noplot = false;
    end
    if nargin < 3
        maxIterations = 100;
    end    
    [n,m] = size(p);
    [N,M] = size(q);
    if n < N
        r = q;
        q = p;
        p = r;
        
        x = n;
        n = N;
        N = x;
    end
    
    worstIndexList = [];
    
    %so p is longer (no shorter) than q at this point
    bestKeepies = 1:N;
    bestLoosies = N+1:n;
    bestDistance = inf;
    stop = false;
    J = 0;
	unmatchedDistance = 0;

    while not(stop) && (J < maxIterations)
    
        dm = distanceMatrix3d(p,q);
        
        %this part is the munkres (hungarian) algorithm that assigns points
        %to each other in pairs
        assignment = assignmentoptimal(dm);
        
        %some of the point can't be assigned if the proteins have different
        %sizes so let's lose the ones that are the worst matches (loosies):
        newLoosies = find(assignment == 0);
        newKeepies = find(assignment > 0);

        %let's find the procustes-like direct isometry that puts the proteins
        %we're keeping as the N best matches on top of the N alpha-carbons in q:
        [R,t] = matchPoints3d(p(newKeepies,:),q);
        pNew = zeros(n,3);
        %now let's shift all the old alpha carbons to their new positions
        for j = 1:n
            pNew(j,:) = R*(p(j,:)') + t;
        end
        %and compute the total distance of the new matches:
        newDistance = sum(totalEuclidean(pNew(newKeepies,:),q));
        %let's stop if the new matches were already the best we found as
        %our algorithm must now have found the best, or be in a loop that
        %cycles from the best one this algorithm can find:
        stop = isempty(setdiff(bestLoosies,newLoosies));
        
        if newDistance < bestDistance
            bestKeepies = newKeepies;
            bestLoosies = newLoosies;
            bestDistance = newDistance;
            if n > N
                unmatchedDM = dm(newLoosies,:);
                unmatchedDistance = (n-N)*mean(mean(abs(unmatchedDM)));
            end
        end
        p = pNew;
        J = J + 1;
        
        bestMatchedDistance = bestDistance;
        bestTotalDistance = bestDistance + unmatchedDistance;
    end
    
    code = J;
    if J >= maxIterations
        code = NaN;
    end
    
    worstIndexList = getWorstIndexList(pNew,q,bestKeepies,bestLoosies,minDistance);
    if not(noplot)
        plist = plotTwoProteins(pNew,q,bestKeepies,bestLoosies,minDistance);
        subplot(1,2,1);
        legend(plist,{['matched \alphaC (',num2str(length(q)),' \alphaCs)'],...
            ['original \alphaC (',num2str(length(p)),' \alphaCs)'],...
            'unmatched \alphaC (due to indel)',...
            ['matches greater than ',num2str(100*minDistance,'%10.1f'),'% total distance']});
        legend('boxoff');
        xlim([-40 40]);
        ylim([-40 40]);
        drawnow
    end
end

function worstIndexList = getWorstIndexList(p,q,pkeepies,ploosies,minDistance)
    
    pSub = p(pkeepies,:);
    worstIndexList = [];
        
	[n,~] = size(p);
	[m,~] = size(q);
    
    %normalise the centre of the proteins for plotting purposes
    pSubBar = mean(pSub,1);
    pSubBarq = repmat(pSubBar,m,1);
    pSubBarp = repmat(pSubBar,n,1);
    pSubBar = repmat(pSubBar,length(pkeepies),1);
    p = p - pSubBarp;
    pSub = pSub - pSubBar;
    q = q - pSubBarq;
    
    n = length(pkeepies);
    missMatches = n;
    [indexes,distances] = worstMatches(pSub,q,missMatches);
    
    for K = 1:missMatches
        j = indexes(K);
        if distances(j) > minDistance*sum(distances)
            worstIndexList = [worstIndexList j];
        end
    end
end

function plotList = plotTwoProteins(p,q,pkeepies,ploosies,minDistance)
    pSub = p(pkeepies,:);
    worstIndexList = [];
    
    fig = gcf;
    set(fig,'color','white')
    
	[n,~] = size(p);
	[m,~] = size(q);
    
    %normalise the centre of the proteins for plotting purposes
    pSubBar = mean(pSub,1);
    pSubBarq = repmat(pSubBar,m,1);
    pSubBarp = repmat(pSubBar,n,1);
    pSubBar = repmat(pSubBar,length(pkeepies),1);
    p = p - pSubBarp;
    pSub = pSub - pSubBar;
    q = q - pSubBarq;
    
    red = [1 0.1 0]*0.7;
    blue = [0.1 0.1 1]*0.7;
    grey = [1 1 1]*0.6;
    green = [0.1 1 0.1]*0.6;
    
    n = length(pkeepies);
    missMatches = n;
    [indexes,distances] = worstMatches(pSub,q,missMatches);
    
    plotList = [];
    for k = 2:-1:1
        subplot(1,2,k)
        
        plotList(1) = plot3(q(:,1),q(:,2),q(:,3),'o','linewidth',1,'markersize',12,'color',red);
        hold on
        plotList(2) = plot3(p(:,1),p(:,2),p(:,3),'.','markersize',16,'color',blue);
        if not(isempty(ploosies))
            plotList(3) = plot3(p(ploosies,1),p(ploosies,2),p(ploosies,3),'o','linewidth',3,'markersize',20,...
                'color',green);
        else
            plotList(3) = plot3(-100,-100,-100,'o','linewidth',6,'markersize',20,'color',green);
        end

        box off
        xlabel('x')
        ylabel('y')
        zlabel('z')
        
        if k <= 2
            for K = 1:missMatches
                j = indexes(K);
                if distances(j) > minDistance*sum(distances)
                    worstIndexList = [worstIndexList j];
                    l=line([q(j,1) pSub(j,1)],[q(j,2) pSub(j,2)],[q(j,3) pSub(j,3)]);
                    set(l,'linewidth',5)
                    set(l,'color','k')
                    plotList(4) = plot3(q(j,1), q(j,2), q(j,3),'ok','markersize',20,'linewidth',4);
                    plot3(pSub(j,1), pSub(j,2), pSub(j,3),'ok','markersize',20,'linewidth',4);
                end
            end
        end
        
        hold on
        plot3(p(:,1),p(:,2),p(:,3),'--','linewidth',1,'color',grey)
        if k == 1
            view(2)
            box on
        end
        if k == 2
            view(36,18)
            box on
        end
        axis([-30 30 -30 30 -30 30])
    end
    
	set(fig,'position',[63         267        1748         683])
    
end