clc
close all

col1 = [1 0.2 0.2];
col2 = [0.2 0.2 1];
Col1 = [0.25 0.25 1];
Col2 = [1 0 0];

%%

dataDistanceParameter = 0.01;
SNPflagTexts = {'allMutations','SNPsOnly','indelsOnly','completeResistanceMutants'};
xSNPflagTexts = {'all mutations','SNPs only','indels only','pan-resistant mutants'};

for SNPflag = [4 1]
    close all

    [alphaCarbonHistogramData,alphaList] = plotAllLamBChanges(data,dataDistanceParameter,SNPflag);
    set(gcf,'pos',[145         436        1146         910])
    view(54,2);

    figure(2)
    set(2,'pos',[742         957        1114         389]);
    bandwidthParameter = 0.05;

    allPoints = alphaList;
    minAllPts = min(allPoints);
    maxAllPts = max(allPoints);
    allPoints = (allPoints - min(allPoints))/(max(allPoints) - min(allPoints));

    bins = 60;%length(allPoints)/10;
    %xpts = length(allPoints);
    xpts = 60;
    [bandwidth,density,xmesh]=kde(allPoints,xpts,min(allPoints),max(allPoints));

    x = min(allPoints):(max(allPoints)-min(allPoints))/bins:max(allPoints);
    n_elements = histc(allPoints,x);

    bar(x,n_elements/sum(n_elements),'BarWidth',1,'facecolor',[0.4 0.15 0.99],'edgecolor','w')
    hold on

    xdata = allPoints;
    X = [zeros(size(xdata)) xdata]';
    [clustCent,point2cluster,clustMembsCell] = MeanShiftCluster(X,bandwidthParameter);

    clHeight = 0.175;
    p = plot(clustCent(2,:),clHeight*ones(size(clustCent(2,:))),'*r','markersize',15,'color',[1 1 1]*0.1);
    legend(p,[num2str(length(clustMembsCell)),' cluster centres \alphaC address']);
    legend('boxoff')

    De = density/sum(density);
    plot(xmesh,De,'-k','color',[1 1 1]*0.1);

    xlim([-0.025 1.025])
    ylim([0 0.22])
    xPts = [0,100,200,300,400,421];
    set(gca,'Xtick',xPts/421);
    set(gca,'Xticklabel',{'0','100','200','300','400','421'});

    alphaClustCent = round(minAllPts + (maxAllPts - minAllPts)*clustCent);
    xlabel('\alpha carbon address in identified WT LamB hotspot')
    ylabel('frequency')

    %%

    figure(3)

    Iw = getLocationLabels('wt',data.bacteria);
    [p1,p2] = plotAProtein2(data.carbonData{Iw});
    hold on

    set(3,'pos',[145         436        1146         910])

    set(p1,'linewidth',2)
    set(p1,'color',[1 1 1]*0.8)
    set(p2,'color',[1 1 1]*0.8)
    set(p2,'markersize',5)

    p = [];
    labels = {'cluster 1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'};
    labels1 = {'cluster centre 1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'};

    for j = 1:length(clustMembsCell)
        s = (j-1)/(length(clustMembsCell) - 1);
        colr = s*Col1 + (1-s)*Col2;
        clusterCarbons = data.carbonData{Iw}(alphaClustCent(2,j),:);
        if SNPflag == 1
            colr = s*col1 + (1-s)*col2;

            p(j) = plot3(clusterCarbons(1),clusterCarbons(2),clusterCarbons(3),'*','markersize',34,'color',colr);
            set(p(j),'color',colr);

            q = plotAProtein(data.carbonData{Iw}(alphaList(clustMembsCell{j}),:),colr);
            set(q,'markersize',30)

            figure(2)
            plot(clustCent(2,j),clHeight,'*','markersize',24,'color',colr);
            figure(3)

        else
            p(j)=plotAProtein(clusterCarbons,colr);
            set(p(j),'marker','.')
            set(p(j),'markersize',60)
            %set(p(j),'linewidth',2)

            q=plotAProtein(data.carbonData{Iw}(alphaList(clustMembsCell{j}),:),colr);
            set(q,'markersize',60)
            %set(q,'marker','.')
        end

    end

    if SNPflag == 1
        legend(p,labels1{1:length(p)},'fontsize',22);
    else
        legend(p,labels{1:length(p)},'fontsize',22);
    end

    legend1 = legend(gca,'show');
    title(legend1,xSNPflagTexts{SNPflag})
    if SNPflag == 1
        set(legend1,'Position',[0.3 0.61 0.15 0.1 + length(p)*0.02]);
    else
        set(legend1,'Position',[0.3 0.61 + length(p)*0.005 0.15 0.1 + length(p)*0.02]);
    end
    legend('boxoff');

    title('')

    box off
    axis off

    view(54,2);

    export_fig(['./figures/HotSpotClusters@',num2str(dataDistanceParameter),'-',SNPflagTexts{SNPflag},'.pdf']);

end