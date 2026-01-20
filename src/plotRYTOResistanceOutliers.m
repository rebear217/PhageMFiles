function KmVec = plotRYTOResistanceOutliers(data,allRKLagData)

    close all

    figure(77)
    set(77,'pos',[637         654        550         417]);
    figure(44)
    set(44,'pos',[637         654        550         417]);
    figure(33)
    set(33,'pos',[637         654        550         417]);
    figure(22)
    set(22,'pos',[637         654        550         417]);
    figure(13)
    set(13,'pos',[637         654        550         417]);
    figure(3)
    set(3,'pos',[637         654        550         417]);
    figure(2)
    set(2,'pos',[637         654        550         417]);
    figure(1)
    set(1,'pos',[637         654        550         417]);
    
    [validMaxRate , validMaxYield, ~ , ~ , ~ , keepRsquaredValueCostPlot ] = validVariables();
    
    red = 0.75*[1 0 0];
    betaMonod = [1 1];
    Bphi = (data.infection.phi > 0);

    susceptibleList = [];
    resistantList = [];
	P = [];
    
    NotConverged = zeros(1,length(allRKLagData));

	MonodModel = @(b,S)((abs(b(1))*S./(1+abs(b(2))*S)));

    Iwt = getLocationLabels('wt',data.bacteria);
	RKdataWT = allRKLagData{Iwt};
	sugars = RKdataWT.sugars;
    delta = (max(sugars)-min(sugars))/20;
	fineSugarsWT = [1:delta:max(sugars)];    
    RdataWT = RKdataWT.Rdata;
    allSugarsWT = RKdataWT.allSugars;
    options = statset('MaxIter',10000);
    MonodmdlWT = NonLinearModel.fit(allSugarsWT,RdataWT,MonodModel,betaMonod,'options',options);
    betaMonodWT = abs(MonodmdlWT.Coefficients.Estimate);    
    seWT = MonodmdlWT.Coefficients.SE;
    MonodRsquaredWT = MonodmdlWT.Rsquared.Adjusted;
    RcurveWT = MonodModel(betaMonodWT,fineSugarsWT);
    
    KmData = [];
    
    MphiVec = zeros(size(1:length(allRKLagData)));
    MrVec = zeros(size(1:length(allRKLagData)));
    MrmaxVec = zeros(size(1:length(allRKLagData)));
    KmVec = zeros(size(1:length(allRKLagData)));
    
    
    for j = 1:length(allRKLagData)
        label = allRKLagData{j}.bacterium;
        if strcmp(label,'wt')
            jWT = j;
        end
        RKdata = allRKLagData{j};

        %remove junk data
        for k = 1:2
            if k == 1
                F = RKdata.Rdata < validMaxRate;
            else
                F = RKdata.Kdata < validMaxYield*RKdata.allSugars;
            end
            
            RKdata.Rsquareddata = RKdata.Rsquareddata(F);
            RKdata.Rdata = RKdata.Rdata(F);
            RKdata.Kdata = RKdata.Kdata(F);
            RKdata.Rsedata = RKdata.Rsedata(F);
            RKdata.Ksedata = RKdata.Ksedata(F);
            RKdata.allSugars = RKdata.allSugars(F);
            RKdata.Yielddata = RKdata.Yielddata(F);
            RKdata.Yieldsedata = RKdata.Yieldsedata(F);
            RKdata.growthRateData = RKdata.growthRateData(F);
            RKdata.growthRateSEData = RKdata.growthRateSEData(F);
            RKdata.Lagdata = RKdata.Lagdata(F);
            RKdata.Lagdata2 = RKdata.Lagdata2(F);
            RKdata.Lagdata3 = RKdata.Lagdata3(F);
            
        end

        Rsquareddata = RKdata.Rsquareddata;
        Rdata = RKdata.Rdata;
        Kdata = RKdata.Kdata;
        Rsedata = RKdata.Rsedata;
        Ksedata = RKdata.Ksedata;
        allSugars = RKdata.allSugars;
        sugars = RKdata.sugars;

        las = log2(allSugars);
        %maxas = max(las);
        %minas = min(las);
        %delta = (max(sugars)-min(sugars))/20;
        fineSugars = fineSugarsWT;

        if j == 1
            meanRcurveDataZEROPhi = zeros(length(allRKLagData),length(fineSugars));
            RcurveData = zeros(length(allRKLagData),length(fineSugars));
            ms = max(sugars);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        figure(1);
        MonodModel = @(b,S)((abs(b(1))*S./(1+abs(b(2))*S)));
        options = statset('MaxIter',1000);
        Monodmdl = NonLinearModel.fit(allSugars,Rdata,MonodModel,betaMonod,'options',options);
        betaMonod = abs(Monodmdl.Coefficients.Estimate);    
        se = Monodmdl.Coefficients.SE;
        MonodRsquared = Monodmdl.Rsquared.Adjusted;
        Rcurve = MonodModel(betaMonod,fineSugars);
        
        % Note that in the usual Monod model rMax*S/(Km+S) the value of Km
        % corresponds to 1/b2 in the model used here which is b1*S/(1+b2*S)
        % This means KmData is actually 1/Km in the usual parlance and the
        % figures are labelled as such:

        KmData = [KmData ; betaMonod(2)];

        NotConverged(j) = (Monodmdl.Rsquared.Adjusted < keepRsquaredValueCostPlot);
        if NotConverged(j)
            disp(['Unreliable data not used : ',label,' only has Monod fit adj R^2 = ',num2str(Monodmdl.Rsquared.Adjusted)]);
        end
            
        %xFit = fineSugars;
        %alpha = 0.05;
        %[betaMonod,resid,J,Sigma] = nlinfit(allSugars, Rdata, MonodModel,betaMonod);
        %[yFit, delta] = nlpredci(MonodModel,xFit,betaMonod,resid,'Covar',Sigma,'Alpha',alpha,'PredOpt','curve');        
        %P = plotshaded(xFit,[yFit-delta ; yFit+delta],'k');
        %set(P,'FaceColor',[1 1 1]*0.8);

        thisColour = 0.75*[1 1 1];
        lw = 1;
        I = getLocationLabels(label,data.infection.bacteria_names);
        if strcmp(label,'wt')
            Iwt = I;
            phiWT = data.infection.phi(:,Iwt);
        end
        phi = data.infection.phi(:,I);
        
        RcurveData(j,:) = Rcurve;
        if sum(phi) == 0
            resistantList = [resistantList j];
        else
            susceptibleList = [susceptibleList j];
            P(1) = plot(fineSugars,Rcurve./RcurveWT,'-','linewidth',1,'color',[1 1 1]*0.75);
        end
        
        hold on

        ylim([0 validMaxRate]);

        ylabel('relative growth rate (unitless)')
        xlabel('malto^{3} (\mug/ml)')
    end
    
    for j = 1:length(resistantList)
    	P(2) = plot(fineSugars,RcurveData(resistantList(j),:)./RcurveWT,'-','linewidth',1,'color',[1 0 0]);
    end
    
	P(3) = plot(fineSugars,RcurveData(jWT,:)./RcurveWT,'-','linewidth',4,'color',[0 0 1]);
    
    box on
    axis tight
    ym = ylim;
    ylim([0 ym(2)*1.2])
    xlim([min(sugars) max(sugars)-1]);
    legend(P,{'susceptible','pan resistant','wild type'},'location','northwest');
    
    P = [];
	figure(2);
    S = 1:length(fineSugars);
    subset = S(1:3:end);
    P(1) = errorbar(fineSugars(subset),mean(RcurveData(susceptibleList,subset),1)./RcurveWT(subset),1.96*std(RcurveData(susceptibleList,subset),1)/sqrt(length(susceptibleList)-1)...
        ,'.','markersize',30,'linewidth',1,'color',[1 1 1]*0.25);
    hold on
	P(5) = plot(fineSugars,RcurveData(jWT,:)./RcurveWT,'-','linewidth',3,'color',[0 0 1]);

	figure(22);
    Q(1) = errorbar(fineSugars(subset),mean(RcurveData(susceptibleList,subset),1)./RcurveWT(subset),1.96*std(RcurveData(susceptibleList,subset),1)/sqrt(length(susceptibleList)-1),'.','markersize',30,'linewidth',1,'color',[1 1 1]*0.25);
    hold on
    
    susceptibleMean = mean(RcurveData(susceptibleList,:),1);
    
	figure(2);
    
    datumR = [];
    datumRmax = [];
    datumPhi = [];
    
    for j = 1:length(allRKLagData)
        label = allRKLagData{j}.bacterium;
        I = getLocationLabels(label,data.infection.bacteria_names);
        if ismember(j,resistantList)
            if all(RcurveData(j,:) > RcurveData(jWT,:)) %grow faster than WT
                text(1.01*ms,RcurveData(j,end)./RcurveWT(end),label,'fontsize',10);
                P(2) = plot(fineSugars,RcurveData(j,:)./RcurveWT,'-','linewidth',3,'color',[1 0 0]);
            else
                %P(4) = plot(fineSugars,RcurveData(j,:),'-','linewidth',3,'color',[1 1 1]*0.25);
            end
        elseif all(RcurveData(j,:) > RcurveData(jWT,:)) %grow faster than WT
            if sum(data.infection.phi(:,I)) < sum(phiWT) %more resistant than WT on average
                P(3) = plot(fineSugars,RcurveData(j,:)./RcurveWT,'-','linewidth',1,'color',[1 0.0 0.0]);
            end
        elseif all(RcurveData(j,:) < RcurveData(jWT,:)) %grow slower than WT
            if sum(data.infection.phi(:,I)) > sum(phiWT) %but more resistant than WT on average
                P(4) = plot(fineSugars,RcurveData(j,:)./RcurveWT,'-','linewidth',1,'color',[0.0 0.0 1]);
            end
        else
            if not(isempty(I))
                s = data.infection.phi(:,I)./phiWT;
                s = s(~isnan(s));
                s = s(~isinf(s));                
                s = mean(s);
                colour = [1-9*s/10 9*s/10 0.5];
                figure(22);
                Q(2) = plot(fineSugars,RcurveData(j,:)./RcurveWT,'-','linewidth',1,'color',[1 1 1]*0.5);
                figure(2);
            end
        end
        errorbar(fineSugars(subset),mean(RcurveData(susceptibleList,subset),1)./RcurveWT(subset),1.96*std(RcurveData(susceptibleList,subset),1)/sqrt(length(susceptibleList)-1),'.','markersize',30,'linewidth',1,'color',[1 1 1]*0.25);
        
        figure(33)        
        try
            if ~NotConverged(j)
                datumR = [datumR mean(RcurveData(j,:))/mean(RcurveData(jWT,:))];
                datumRmax = [datumRmax max(RcurveData(j,:))/max(RcurveData(jWT,:))];
                datumPhi = [datumPhi mean(data.infection.phi(:,I))/mean(phiWT)];

                Mphi = mean(data.infection.phi(:,I))/mean(phiWT);
                Mr = mean(RcurveData(j,:))/mean(RcurveData(jWT,:));
                MphiVec(j) = Mphi;
                MrVec(j) = Mr;

                colr = [1 1 1]*0.7;

                p1 = plot(Mr,Mphi,'.','markersize',40,'color',colr);
                hold on
                if strcmp(data.infection.bacteria_names(I),'26a') || ...
                        strcmp(data.infection.bacteria_names(I),'2b') || ...
                        strcmp(data.infection.bacteria_names(I),'56a')
                    p11 = plot(Mr,Mphi,'o','markersize',22,'color','r','linewidth',4);
                    text((1+rand/10)*Mr,0.05*(1+rand),data.infection.bacteria_names(I))
                end            
            end
            
        catch
            %just do nothing if an error occurs...
        end
        hold on

        figure(44)        
        try            
            if ~NotConverged(j)
                Mphi = mean(data.infection.phi(:,I))/mean(phiWT);
                Km = KmData(j)/KmData(jWT);
                KmVec(j) = Km;

                colr = [1 1 1]*0.7;

                p12 = plot(Km,Mphi,'.','markersize',40,'color',colr);
                hold on
                if strcmp(data.infection.bacteria_names(I),'56a')
                    p112 = plot(Km,Mphi,'o','markersize',22,'color','r','linewidth',4);
                    text((1+rand/10)*Km,0.05*(1+rand),data.infection.bacteria_names(I))
                end
                if strcmp(data.infection.bacteria_names(I),'2b')
                    p1122 = plot(Km,Mphi,'o','markersize',22,'color','b','linewidth',4);
                    text((1+rand/10)*Km,0.05*(1+rand),data.infection.bacteria_names(I))
                end
                if strcmp(data.infection.bacteria_names(I),'26a')
                    p11222 = plot(Km,Mphi,'o','markersize',22,'color','r','linewidth',4);
                    text((1+rand/10)*Km,0.05*(1+rand),data.infection.bacteria_names(I))
                end
            end
            
        catch
            %just do nothing if an error occurs...
        end
        hold on

        
        figure(77)        
        try
            if ~NotConverged(j)
                Mphi = mean(data.infection.phi(:,I))/mean(phiWT);
                Mrmax = max(RcurveData(j,:))/max(RcurveData(jWT,:));
                MrmaxVec(j) = Mrmax;

                colr = [1 1 1]*0.7;

                p71 = plot(Mrmax,Mphi,'.','markersize',40,'color',colr);
                hold on
                if strcmp(data.infection.bacteria_names(I),'26a') || ...
                        strcmp(data.infection.bacteria_names(I),'2b') || ...
                        strcmp(data.infection.bacteria_names(I),'56a')

                    p711 = plot(Mrmax,Mphi,'o','markersize',22,'color','r','linewidth',4);
                    text((1+rand/10)*Mrmax,0.05*(1+rand),data.infection.bacteria_names(I))
                end
            end
        catch
            %just do nothing if an error occurs...
        end
        hold on
        
        figure(2)
    end
    
	figure(33)
	plot([-3 3],[1 1],'-k','color',[1 1 1]*0.7,'linewidth',1)
    hold on
    plot([1 1],[-3 3],'-k','color',[1 1 1]*0.7,'linewidth',1)

    p2 = errorbar(mean(MrVec),mean(MphiVec),1.96*ste(MphiVec),'.','markersize',45,'linewidth',1,'color',[1 1 1]*0.1);
    hold on
    h = herrorbar(mean(MrVec),mean(MphiVec),1.96*ste(MrVec));
    set(h,'linewidth',1);
    set(h,'color',[1 1 1]*0.1);
    xlabel('\rho( r_{mean} )')
    ylabel('\rho ( \phi )')
    p3 = plot(1,1,'.b','markersize',45);
    legend([p1,p2,p3,p11],{'r_{mean} data','mean \pm 95% CI','wt','26a, 56a & 2b'},'location','northwest');
    xlim([-0.0 1.8])
    ylim([-0.02 1.2])       
    
	figure(77)
	plot([-3 3],[1 1],'-k','color',[1 1 1]*0.7,'linewidth',1)
    hold on
    plot([1 1],[-3 3],'-k','color',[1 1 1]*0.7,'linewidth',1)

    p72 = errorbar(mean(MrmaxVec),mean(MphiVec),1.96*ste(MphiVec),'.','markersize',45,'linewidth',1,'color',[1 1 1]*0.1);
    hold on
    h = herrorbar(mean(MrmaxVec),mean(MphiVec),1.96*ste(MrmaxVec));
    set(h,'linewidth',1);
    set(h,'color',[1 1 1]*0.1);
    xlabel('\rho( r_{max} )')
    ylabel('\rho ( \phi )')
    p73 = plot(1,1,'.b','markersize',45);
    legend([p71,p72,p73,p711],{'r_{max} data','mean \pm 95% CI','wt','26a, 56a & 2b'},'location','northwest');
    xlim([-0.0 1.8])
    ylim([-0.02 1.2])
    

	figure(44)    

    Km = KmData;
	plot([-4 4],[1 1],'-k','color',[1 1 1]*0.7,'linewidth',1)
    hold on
    plot([1 1],[-4 4],'-k','color',[1 1 1]*0.7,'linewidth',1)
    p22 = errorbar(mean(KmVec),mean(MphiVec),1.96*ste(MphiVec),'.','markersize',45,'linewidth',1,'color',[1 1 1]*0.1);
    hold on
    h = herrorbar(mean(KmVec),mean(MphiVec),1.96*ste(KmVec));
        
    set(h,'linewidth',1);
    set(h,'color',[1 1 1]*0.1);
    xlabel('\rho( 1/K_m )')
    ylabel('\rho ( \phi )')    
    p32 = plot(1,1,'.b','markersize',45);    
    legend([p12,p22,p32,p112,p1122,p11222],...
        {'K_{m} data','mean \pm 95% CI','wt','56a','2b','26a'},'location','northeast');

    xlim([-0.04 3])
    ylim([-0.04 1.2])   
     

	figure(2);
    axis tight
    ym = ylim;
    ylim([0.5 ym(2)*1.2])
    xlim([min(sugars) max(sugars)-1]);
    
    legend(P,{'susceptibles mean \pm 95% CI','pan resistant > WT growth',...
        '> resistant > WT growth ','> resistant < WT growth',...
        'wild type'},'location','northwest');

    ylabel('relative growth rate (unitless)')
    xlabel('malto^{3} (\mug/ml)')
    
	figure(22);
    Q(3) = plot(fineSugars,RcurveData(jWT,:)./RcurveWT,'-','linewidth',3,'color',[0 0 1]);
    axis tight
    ym = ylim;
    ylim([0 ym(2)*1.2])
    xlim([min(sugars) max(sugars)-1]);
    ylim([0.6 1.6])
    
    legend(Q,{'mean over susceptibles \pm 95% CI',...
        '> WT resistance',...
        'wild type'},'location','northwest');

    ylabel('relative growth rate (unitless)')
    xlabel('malto^{3} (\mug/ml)')
        
    figure(3);
    resistantAUCdata = RcurveData(resistantList,:);
    susceptibleAUCdata = RcurveData(susceptibleList,:);
    anovaData = nan(length(susceptibleList),2);
    anovaData(1:length(resistantList),1) = mean(RcurveData(resistantList,:),2);
    anovaData(:,2) = mean(RcurveData(susceptibleList,:),2);
    boxplot(anovaData,{'zero','positive'});

    box on
    ylabel('mean growth rate (per h)')
    xlabel('phage susceptibility group')
    %p = anova1(anovaData);
    p = kruskalwallis(anovaData);
            
	figure(3)
    ylim([0 ym(2)*1.2])

    text(0.605,1.1,['Kruskal-Wallis p = ',num2str(p,3)]);
    hold on
    plot(2.12,mean(RcurveData(jWT,:)),'*','markersize',10,'linewidth',1);
    text(2.2,mean(RcurveData(jWT,:)),'wt','fontsize',10);
    
    for j = 1:length(allRKLagData)
        label = allRKLagData{j}.bacterium;
        if ismember(j,resistantList)
            if sum(RcurveData(j,:)) > sum(susceptibleMean)
                text(1.2,mean(RcurveData(j,:)),label,'fontsize',10);
                P(2) = plot(1.12,mean(RcurveData(j,:)),'*','markersize',10,'color',red,'linewidth',1);
            end
            plot(1.0 + 0.05*(rand-0.5),mean(RcurveData(j,:)),'.','markersize',20,'color',[1 1 1]*0.75,'linewidth',1);
        else
            plot(2.0 + 0.05*(rand-0.5),mean(RcurveData(j,:)),'.','markersize',20,'color',[1 1 1]*0.75,'linewidth',1);
        end
    end
    ylim([0 1.2])

    figure(13);
    resistantAUCdata = RcurveData(resistantList,:);
    susceptibleAUCdata = RcurveData(susceptibleList,:);
    anovaData = nan(length(susceptibleList),2);
    anovaData(1:length(resistantList),1) = max(RcurveData(resistantList,:),[],2);
    anovaData(:,2) = max(RcurveData(susceptibleList,:),[],2);
    boxplot(anovaData,{'zero','positive'});

    box on
    ylabel('max growth rate (per h)')
    xlabel('phage susceptibility group')
    %p = anova1(anovaData);
	p = kruskalwallis(anovaData);
        
    figure(13);
    ylim([0 ym(2)*1.2])

    text(0.605,2.2,['Kruskal-Wallis p = ',num2str(p,3)]);
    hold on
    plot(2.12,max(RcurveData(jWT,:)),'*','markersize',10,'linewidth',1);
    text(2.2,max(RcurveData(jWT,:)),'wt','fontsize',10);
    
    for j = 1:length(allRKLagData)
        label = allRKLagData{j}.bacterium;
        if ismember(j,resistantList)
            if sum(RcurveData(j,:)) > sum(susceptibleMean)
                text(1.2,max(RcurveData(j,:)),label,'fontsize',10);
                P(2) = plot(1.12,max(RcurveData(j,:)),'*','markersize',10,'color',red,'linewidth',1);
            end
            plot(1.0 + 0.05*(rand-0.5),max(RcurveData(j,:)),'.','markersize',20,'color',[1 1 1]*0.75,'linewidth',1);
        else
            plot(2.0 + 0.05*(rand-0.5),max(RcurveData(j,:)),'.','markersize',20,'color',[1 1 1]*0.75,'linewidth',1);
        end
    end
    ylim([0 2.5])    
    
    close(1)
    close(2)
    close(3)
    close(13)
    close(6)
    close(4)
    close(7)
	close(5)
    close(22)
    
    figure(33)
    export_fig('./figures/findResistantGrowthRateTOsAndTUs.pdf');
        
    figure(44)
    export_fig('./figures/findResistantKmTOsAndTUs.pdf');
    
    figure(77)
    export_fig('./figures/findResistantRmaxTOsAndTUs.pdf');

    
    thisSheet = {};
	thisSheet{1,1} = 'rho(R) v rho(phi)';
    thisSheet{2,1} = 'x axis data (rho-R)';
    thisSheet{3,1} = 'y axis data (rho-phi)';
    thisSheet{2,2} = MrVec(:)';
    thisSheet{3,2} = MphiVec(:)';
    writecell(thisSheet,'./dataRepo/derivedData/WTRelativeTraits/derivedRelativeTraitData.xlsx','Sheet','rho(R)Vrho(phi)');
    
    thisSheet = {};
	thisSheet{1,1} = 'rho(Rmax) v rho(phi)';
    thisSheet{2,1} = 'x axis data (rho-Rmax)';
    thisSheet{3,1} = 'y axis data (rho-phi)';
    thisSheet{2,2} = MrmaxVec(:)';
    thisSheet{3,2} = MphiVec(:)';
    writecell(thisSheet,'./dataRepo/derivedData/WTRelativeTraits/derivedRelativeTraitData.xlsx','Sheet','rho(Rmax)Vrho(phi)');
    
    thisSheet = {};
	thisSheet{1,1} = 'rho(Km) v rho(phi)';
    thisSheet{2,1} = 'x axis data (rho-Km)';
    thisSheet{3,1} = 'y axis data (rho-phi)';
    thisSheet{2,2} = KmVec(:)';
    thisSheet{3,2} = MphiVec(:)';
    writecell(thisSheet,'./dataRepo/derivedData/WTRelativeTraits/derivedRelativeTraitData.xlsx','Sheet','rho(Km)Vrho(phi)');
    
end
