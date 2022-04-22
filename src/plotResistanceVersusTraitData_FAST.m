function Outdata = plotResistanceVersusTraitData_FAST(data,traitCode,allRKLagData,lagCode)

    warning('off')
    
    if nargin < 4
        lagCode = 1;
    end
    
    if nargin < 2
        traitCode = 1;
    end
    
    close all
    
    red = 0.5*[1 0 0];
    colourFunction = @(s)[1-s 0.5 1-s];
    
    [validMaxRate , validMaxYield , ~ , ~ , sig] = validVariables();
    
    MaxLag = 24;
    MaxPhi = 0.5;
    
	figure(1);
    set(1,'pos',[416         474        2085         872]);
    figure(2);
    set(2,'pos',[416     1   521   399]);
    
    sugars = unique(allRKLagData{1}.allSugars);
    for k = 1:length(sugars)
        totalXdataSet{k} = [];
        totalYdataSet{k} = [];
    end
    P = zeros(1,length(sugars));
    lFlag1 = P;
        
    for j = 1:length(data.bacteria)
        
        %RKLagdata = getRKdataFromLabel(data.bacteria{j},data,keepRsquaredValue);
        %this should be in the same order as the file loaded, but be careful:
        RKLagdata = allRKLagData{j};
        if ~strcmp(RKLagdata.bacterium,data.bacteria{j})
            error('plotResistanceVersusTraitData has used an incorrect bacterial label')
        end
        
        %sugars = unique(RKLagdata.allSugars);
        
        I = getLocationLabels(data.bacteria{j},data.infection.bacteria_names);

        for k = 1:length(sugars)
            las = log2(sugars);
            maxas = max(las);
            minas = min(las);
            
            F = RKLagdata.allSugars == sugars(k);
            
            infectionVector = data.infection.phi(:,I);
            %figure(k);
            
            switch traitCode
                case 1
                    if lagCode == 1
                        XdataSet = RKLagdata.Lagdata(F);
                    else
                        XdataSet = RKLagdata.Lagdata2(F);
                    end
                    YdataSet = mean(infectionVector);
                    YdataSet = ones(size(XdataSet))*mean(YdataSet);
                    xText = 'lag (h)';
                    yText = 'mean phage susceptibility';
                    axisV = [0 MaxLag 0 MaxPhi];
                case 2
                    XdataSet = RKLagdata.Rdata(F);
                    %XdataSet = RKLagdata.growthRateData(F);
                    YdataSet = mean(infectionVector);
                    YdataSet = ones(size(XdataSet))*mean(YdataSet);
                    xText = 'growth rate (per h)';
                    yText = 'mean phage susceptibility';
                    axisV = [0 validMaxRate 0 MaxPhi];
                case 3
                    XdataSet = RKLagdata.Kdata(F)/sugars(k);
                    YdataSet = mean(infectionVector);
                    YdataSet = ones(size(XdataSet))*mean(YdataSet);
                    xText = 'yield (OD_{600} per malto^3)';
                    yText = 'mean phage susceptibility';
                    axisV = [0 validMaxYield 0 MaxPhi];
                case 4
                    XdataSet = RKLagdata.Rdata(F);
                    %XdataSet = RKLagdata.growthRateData(F);
                    YdataSet = RKLagdata.Kdata(F)/sugars(k);
                    xText = 'growth rate (per h)';
                    yText = 'yield (OD_{600} per malto^3)';
                    axisV = [0 validMaxRate 0 validMaxYield];
                case 5
                    XdataSet = RKLagdata.Rdata(F);
                    %XdataSet = RKLagdata.growthRateData(F);
                    if lagCode == 1
                        YdataSet = RKLagdata.Lagdata(F);
                    else
                        YdataSet = RKLagdata.Lagdata2(F);
                    end
                    xText = 'growth rate (per h)';
                    yText = 'lag (h)';
                    axisV = [0 validMaxRate 0 MaxLag];
                case 6
                    if lagCode == 1
                        XdataSet = RKLagdata.Lagdata(F);
                    else
                        XdataSet = RKLagdata.Lagdata2(F);
                    end
                    YdataSet = RKLagdata.Kdata(F)/sugars(k);
                    xText = 'lag (h)';
                    yText = 'yield (OD_{600} per malto^3)';
                    axisV = [0 MaxLag 0 validMaxYield];
            end
            
            fx = find(XdataSet < axisV(2));
            XdataSet = XdataSet(fx);
            YdataSet = YdataSet(fx);
            fx = find(XdataSet >= 0);
            XdataSet = XdataSet(fx);
            YdataSet = YdataSet(fx);
            
            fy = find(YdataSet < axisV(4));
            XdataSet = XdataSet(fy);
            YdataSet = YdataSet(fy);
            fy = find(YdataSet >= 0);
            XdataSet = XdataSet(fy);
            YdataSet = YdataSet(fy);
            
            totalXdataSet{k} = [totalXdataSet{k} XdataSet];
            totalYdataSet{k} = [totalYdataSet{k} YdataSet];
                    
        end

    end

    for k = 1:length(sugars)
        
        s = (k-1)/(length(sugars)-1);
        thisColour = colourFunction(s);        

        XdataSet = totalXdataSet{k};
        YdataSet = totalYdataSet{k};        
        [XdataSet,YdataSet] = removeOutliers(XdataSet,YdataSet,sig);
        
        for f = 1:2
            if f == 1
                colr = [1 1 1]*0.7;
            else
                colr = thisColour;
            end

            figure(f);
            box on
            if f == 1
                subplot(2,4,k);
                box on
            end

            if f == 1
                thisP = plot(XdataSet,YdataSet,'.','color',colr,'MarkerSize',40);
                if ~isempty(thisP)
                    P(k) = thisP;
                    set(get(get(thisP,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                end
            else
                pl(k) = plot(XdataSet,YdataSet,'.','color',colr,'MarkerSize',40);
                set(get(get(pl(k),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end
            hold on

            if f == 1
                subplot(2,4,k);
                if ~isempty(thisP) && lFlag1(k) == 0
                    xlabel(xText);
                    ylabel(yText);
                    legend(['data for malto^3 @ ', num2str(sugars(k)),' \mug/ml']);
                    lFlag1(k) = 1;
                end
            end

            axis(axisV);
            drawnow
            
        end
    end
    
    figure(2)
    legend(pl,data.sugarString2{2:end});
    xlabel(xText);
    ylabel(yText);
                
    figure(1)
    for k = 1:length(sugars)
        subplot(2,4,k);
        box on
        axis tight
        xl = xlim;
        xlim([0 1.5*xl(2)]);
        yl = ylim;
        ylim([0 1.5*yl(2)]);
    end
    
    figure(1)
    QD = [];
    pci = [];

	demingmodel = @(b,x)(b(1) + b(2)*x);

    for k = 1:length(sugars)
        subplot(2,4,k);
        
        linearmdl = fitlm(totalXdataSet{k},totalYdataSet{k},'linear','RobustOpts','on');
        Coefficients = linearmdl.Coefficients.Estimate;
        Rsquared = linearmdl.Rsquared.Adjusted;
        pValue = linearmdl.Coefficients.pValue;
        X = [0 3];
        
        XdataSet = totalXdataSet{k}';
        YdataSet = totalYdataSet{k}';        
        [XdataSet,YdataSet] = removeOutliers(XdataSet,YdataSet,sig);
        
        [B,sigma2_x,x_est,y_est,stats] = deming(XdataSet,YdataSet);
        QD(k) = plot(X,demingmodel(B,X),'-k','color',[1 1 1]*0.25,'linewidth',3);        
        
        CI = ['[',num2str(stats.b_ci(2,1),2),',',num2str(stats.b_ci(2,2),2),']'];        
        legend([P(k),QD(k)],{['malto^3 @ ', num2str(sugars(k)),' \mug/ml'],...
        ['deming slope 95% CI ',CI]});

        
        [CrhoS,CpS] = corr(XdataSet,YdataSet,'type','Spearman');
        [CrhoP,CpP] = corr(XdataSet,YdataSet,'type','Pearson');
        
        switch k
            case 1
                xlim([0 0.15]);
                ylim([0 2e-3]);
            case 2
                xlim([0 0.18]);
                ylim([0 2e-3]);
            case 3
                xlim([0 0.3]);
                ylim([0 2e-3]);
            case 4
                xlim([0 0.4]);
                ylim([0 1e-3]);
            case 5
                xlim([0 0.5]);
                ylim([0 0.7e-3]);
            case 6
                xlim([0 1]);
                ylim([0 0.7e-3]);
            case 7
                xlim([0 2]);
                ylim([0 0.7e-3]);
            case 8                
                xlim([0 2.5]);
                ylim([0 0.7e-3]);
        end
        
        YL = ylim;
        XL = xlim;
        text(0.35*XL(2),0.11*YL(2),['Pearson \rho \approx ',num2str(CrhoP,2),', p \approx ',num2str(CpP,2)]);
        text(0.35*XL(2),0.06*YL(2),['Spearman \rho \approx ',num2str(CrhoS,2),', p \approx ',num2str(CpS,2)]);
        
    end
    
    Outdata.Xtrait = totalXdataSet;
    Outdata.Ytrait = totalYdataSet;
    Outdata.sugars = sugars;
    Outdata.Xmeaning = xText;
    Outdata.Ymeaning = yText;

end

function [Xd,Yd] = removeOutliers(Xd,Yd,sigma)

    mX = mean(Xd);
    mY = mean(Yd);
    sX = std(Xd);
    sY = std(Yd);        

    %outlier removal at sD x sigma:
    F = Xd < mX + sigma*sX;
    Xd = Xd(F);
    Yd = Yd(F);
    F = mX - sigma*sX < Xd;
    Xd = Xd(F);
    Yd = Yd(F);

    F = Yd < mY + sigma*sY;
    Xd = Xd(F);
    Yd = Yd(F);
    F = mY - sigma*sY < Yd;
    Xd = Xd(F);
    Yd  = Yd(F);        

end
