function [coefficients,se,AIC,RKdata,RYplotLine,outGuesses] = plotJustRYTOfromLabel_FAST(label,allRKLagData,plotFlag,guesses)

    opts = statset('nlinfit');
    %opts.RobustWgtFun = 'cauchy';
    opts.RobustWgtFun = 'fair';
    %opts.RobustWgtFun = 'welsch';

    warning('off')
    if nargin < 3
        plotFlag = 0;
    end
    if nargin < 4
        guesses = [];
    end
    
    red = 0.5*[1 0 0];
    colourFunction = @(s)[1 1 1]*(0.6*s + 0.1);

    J = NaN;
    for j = 1:46
        if strcmp(allRKLagData{j}.bacterium,label)
            J = j;
            RKdata = allRKLagData{j};
        end
    end
    
    [validMaxRate , validMaxYield , validMaxCoV] = validVariables();
    
    texts = {'high yield CoV','high K CoV','max Yield breach','max Rate Breach'};
    for j = 1:4
        
        switch j
            case 1
                %option to remove noisy data based on Yield CoV:        
                F = RKdata.Yieldsedata ./ RKdata.Yielddata < validMaxCoV;
            case 2
                %option to remove noisy data based on K CoV:        
                F = RKdata.Ksedata ./ RKdata.Kdata < validMaxCoV;
            case 3
                F = RKdata.Yielddata < validMaxYield;                
            case 4
                F = RKdata.Rdata < validMaxRate;                
        end
        
        if not(all(F))
            disp([label,' - ',texts{j},' : ignoring ',...
                num2str(length(find(F == 0))),' / ',num2str(length(F)),' datapoint(s)']);                
        end
        
        RKdata.Rsquareddata = RKdata.Rsquareddata(F);
        RKdata.Rdata = RKdata.Rdata(F);        
        RKdata.Kdata = RKdata.Kdata(F);
        RKdata.Rsedata = RKdata.Rsedata(F);
        RKdata.Ksedata = RKdata.Ksedata(F);
        RKdata.allSugars = RKdata.allSugars(F);
        RKdata.Yielddata = RKdata.Yielddata(F);
        RKdata.Yieldsedata = RKdata.Yieldsedata(F);
        RKdata.Lagdata = RKdata.Lagdata(F);
        RKdata.Lagdata2 = RKdata.Lagdata2(F);
        RKdata.Lagdata3 = RKdata.Lagdata3(F);
        RKdata.growthRateData = RKdata.growthRateData(F);
        RKdata.growthRateSEData = RKdata.growthRateSEData(F);        
    end
        
    Rsquareddata = RKdata.Rsquareddata;
    Rdata = RKdata.Rdata;
    Kdata = RKdata.Kdata;
    Rsedata = RKdata.Rsedata;
    Ksedata = RKdata.Ksedata;
    Lagdata = RKdata.Lagdata2;
    allSugars = RKdata.allSugars;
    sugars = RKdata.sugars;
    
	las = log2(allSugars);
    maxas = max(las);
    minas = min(las);
    
    Yieldsedata = RKdata.Yieldsedata;
    Yielddata = RKdata.Yielddata;
    
    cHiGuess = max(Yielddata);
    cLoGuess = min(Yielddata);

    RYTOmodelY = @(b,S)((abs(b(1))+abs(b(2))*abs(b(3)).*S)./(1+abs(b(3))*S));
    if isempty(guesses)
        RYTOmdl = NonLinearModel.fit(allSugars,Yielddata,RYTOmodelY,[cHiGuess cLoGuess 1],'Options',opts);
    else
        RYTOmdl = NonLinearModel.fit(allSugars,Yielddata,RYTOmodelY,guesses.betaRYTOguess,'Options',opts);
    end
    betaRYTO = abs(RYTOmdl.Coefficients.Estimate);
    seRYTO = RYTOmdl.Coefficients.SE;
    
    outGuesses.betaRYTOguess = betaRYTO;

    %delta = (max(sugars)-min(sugars))/1000;
    %fineSugars = [0.1:delta:(1.5*max(sugars))];
    fineSugars = 0.1:0.01:500;

    RYTORsquared = RYTOmdl.Rsquared.Adjusted;
    
	linearmodel = @(b,x)(b(1)*ones(size(x)));
    linearmdl = NonLinearModel.fit(allSugars,Yielddata,linearmodel,[cLoGuess],'Options',opts);
    constantYield = linearmdl.Coefficients.Estimate(1);
    ConstantRsquared = linearmdl.Rsquared.Adjusted;
    
    AIC.constantmodel = linearmdl.ModelCriterion.AICc;
    AIC.RYTOmodel = RYTOmdl.ModelCriterion.AICc;
    
    coefficients.cLoSugar = betaRYTO(1);
    coefficients.cHiSugar = betaRYTO(2);
    coefficients.p = betaRYTO(3);

    se.cLoSugar = seRYTO(1);
    se.cHiSugar = seRYTO(2);
    se.p = seRYTO(3);
    
	AIC.RL = exp(-abs(AIC.RYTOmodel - AIC.constantmodel)/2);        

    if plotFlag
        [betaYield2,resid,J,Sigma] = nlinfit(allSugars, Yielddata, RYTOmodelY,betaRYTO,'Options',opts);
        %xFit = fineSugars;
        alpha = 0.05;
        [yieldFit2, Yielddelta] = nlpredci(RYTOmodelY,fineSugars,betaYield2,resid,'Covar',Sigma,'Alpha',alpha,'PredOpt','curve');            
        
        %{
        figure(1)
        P = [];

        semilogx(allSugars,Yielddata,'.','markersize',10,'color',0.75*[1 1 1]);
        hold on
        P(4) = plotshaded(xFit,[yFit-delta ; yFit+delta],'k');
        set(P(4),'FaceColor',[1 1 1]*0.8); 
        semilogx(allSugars,Yielddata,'.','markersize',10,'color',0.75*[1 1 1]);
                
        for j = 1:length(allSugars)
            s = (las(j)-minas)/(maxas - minas);
            thisColour = colourFunction(s);
            if j == length(allSugars)
                P(1) = errorbar(allSugars(j),Yielddata(j),Yieldsedata(j),'.','markersize',30,'color',thisColour,'linewidth',1);
            else
                errorbar(allSugars(j),Yielddata(j),Yieldsedata(j),'.','markersize',30,'color',thisColour,'linewidth',1);
            end
        end
        
        P(2) = plot(fineSugars,coefficients.cHiSugar*ones(size(fineSugars)),'-','linewidth',1,'color',0.75*[1 1 1]);
        P(3) = plot(fineSugars,RYTOmodelY(betaRYTO,fineSugars),'-','linewidth',2,'color',0.25*[1 1 1]);
        P(5) = plot(fineSugars,constantYield*ones(size(fineSugars)),'--','linewidth',1,'color',red);
        
        axis tight
        ym = ylim; ym = max(Yielddata);
        ylim([0 1.6*ym]);

        ylabel('yield (OD_{600} per malto^{3})')
        xlabel('malto^{3} (\mug/ml)')
        
        legend(P,{[label,' data'],'high sugar asymptote',['nonlinear regression (adj R^2 \approx ',num2str(RYTORsquared,3),')'],...
            '95% CIs',['constant regression (adj R^2 \approx ',num2str(ConstantRsquared,3),', RL \approx ',num2str(AIC.RL,3),')']},'location','NorthWest');
        legend('boxoff')
        %}
        
        %{
        figure(2)
        P = [];
        plot(allSugars,Rdata,'.','markersize',10,'color',0.75*[1 1 1]);
        hold on
        %}
        
        RYTOmodelR = @(b,S)((abs(b(1))*S./(1+abs(b(2))*S)).*(1+abs(b(3))*S)./(1+(abs(b(3))+abs(b(4)))*S));
        if isempty(guesses)
            MonodModel = @(b,S)((abs(b(1))*S./(1+abs(b(2))*S)));
            Monodmdl = NonLinearModel.fit(allSugars,Rdata,MonodModel,[max(Rdata) 1]);
            betaMonod = abs(Monodmdl.Coefficients.Estimate);    
            RYTOmdlR = NonLinearModel.fit(allSugars,Rdata,RYTOmodelR,[betaMonod(1) betaMonod(2) 1 0]);
        else
            RYTOmdlR = NonLinearModel.fit(allSugars,Rdata,RYTOmodelR,guesses.betaGuess,'Options',opts);            
        end
        beta = abs(RYTOmdlR.Coefficients.Estimate);
        se = RYTOmdlR.Coefficients.SE;
        
        outGuesses.betaGuess = beta;

        %{
        MonodRsquared = Monodmdl.Rsquared.Adjusted;
        RYTORsquared = RYTOmdlR.Rsquared.Adjusted;
        
        AICMonodmodel= Monodmdl.ModelCriterion.AICc;
        AICRYTOmodel= RYTOmdlR.ModelCriterion.AICc;
        RL = exp(-abs(AICRYTOmodel - AICMonodmodel)/2);
        %}
        
        [betaRate2,resid,J,Sigma] = nlinfit(allSugars, Rdata, MonodModel,betaMonod,'Options',opts);
        [rateFit2, Ratedelta] = nlpredci(MonodModel,fineSugars,betaRate2,resid,'Covar',Sigma,'Alpha',alpha,'PredOpt','curve');        
        
        %{
        P(3) = plotshaded(xFit,[yFit-RateDelta ; yFit+RateDelta],'k');
        set(P(3),'FaceColor',[1 1 1]*0.8);
        
        for j = 1:length(allSugars)
            s = (las(j)-minas)/(maxas - minas);
            thisColour = colourFunction(s);
            if j == length(allSugars)
                P(1) = errorbar(allSugars(j),Rdata(j),Rsedata(j),'.','markersize',30,'color',thisColour,'linewidth',1);
            else
                errorbar(allSugars(j),Rdata(j),Rsedata(j),'.','markersize',30,'color',thisColour,'linewidth',1);
            end
        end
        
        P(2) = plot(fineSugars,RYTOmodelR(beta,fineSugars),'-','linewidth',1,'color',0.5*[1/3 1/3 1]);
        P(4) = plot(fineSugars,MonodModel(betaMonod,fineSugars),'-','linewidth',1,'color',0.5*[1 1/3 1/3]);
        
        axis tight
        ym = ylim; ym = max(Rdata);
        ylim([0 1.6*ym]);

        ylabel('growth rate (per h)')
        xlabel('malto^{3} (\mug/ml)')
        legend(P,{[label,' data'],['Monod regression (adj R^2 \approx ',num2str(MonodRsquared,3),')'],...
            '95% CIs',['RYTO-adjusted regression (adj R^2 \approx ',num2str(RYTORsquared,3),', RL \approx ',num2str(RL,3),')']},'location','NorthWest');
        legend('boxoff')
        
        figure2Data.allSugars = allSugars;
        figure2Data.Rdata = Rdata;
        figure2Data.fineSugars = fineSugars;
        figure2Data.RYTOmodelR = RYTOmodelR(beta,fineSugars);
        %}
        
        P = [];
        semilogx(Rdata,Yielddata,'.','markersize',24,'color',0.75*[1 1 1]);
        hold on
        heb = herrorbar(Rdata,Yielddata,1.96*Rsedata,'.');
        errorbar(Rdata,Yielddata,1.96*Yieldsedata,'.','markersize',24,'color',0.5*[1 1 1],'linewidth',1);
        set(heb,'color',0.5*[1 1 1])
        set(heb,'linewidth',1);
        
        %E=errorbarxy(Rdata,Yielddata,Rsedata,Yieldsedata,{'k.', 0.5*[1 1 1], 0.5*[1 1 1]});        
        %set(E.hErrorbar(1),'linewidth',1);
        %set(E.hErrorbar(2),'linewidth',1);
        %set(E.hErrorbar(3),'linewidth',1);
        %set(E.hErrorbar(4),'linewidth',1);
        %set(E.hErrorbar(5),'linewidth',1);
        %set(E.hErrorbar(6),'linewidth',1);
        
        for j = 1:length(allSugars)
            s = (las(j)-minas)/(maxas - minas);
            thisColour = colourFunction(s);
            if j == length(allSugars)
                P(1) = plot(Rdata(j),Yielddata(j),'Marker','.','MarkerSize',34,'color',thisColour);
            else
                plot(Rdata(j),Yielddata(j),'Marker','.','MarkerSize',34,'color',thisColour);
            end
        end
        
        RY = RYTOmodelR(beta,fineSugars);
        RYY = RYTOmodelY(betaRYTO,fineSugars);        
        %P(2) = plot(RY,RYY,'-','linewidth',2,'color',0.25*[1 1 1]);                
        
        P(2) = plot(rateFit2,yieldFit2,'-','color',0.25*[1 1 1],'linewidth',2);
        
        set(gca,'Xtick',[0.01 0.1 1 10]);
        axis tight
        ym = max(Yielddata);
        ylim([0 1.6*ym]);
        xlim([0 max(RY)])

        plot(rateFit2 + Ratedelta,yieldFit2 + Yielddelta,'--','color',[1 1 1]/2,'linewidth',1);
        pd = plot(rateFit2 - Ratedelta,yieldFit2 - Yielddelta,'--','color',[1 1 1]/2,'linewidth',1);
                
        legend([P pd],{[label,' data'],'RY theory','95% CIs'},'location','NorthEast');
        %legend('boxoff')
        
        RYplotLine.x = RY;
        RYplotLine.y = RYY;
        
        %ylabel('yield')
        %xlabel('growth rate (per h)')
        
        %{
        figure(4)
        P = [0 0 0];
        
    	lagModel = @(b,S)(b(1) + b(2)*S.^(abs(b(3))));
        lagMdl = NonLinearModel.fit(allSugars,Lagdata,lagModel,[1 1 1]);
        lagCoeff = [lagMdl.Coefficients.Estimate(1:2) ; abs(lagMdl.Coefficients.Estimate(3))];
        
        [beta,resid,J,Sigma] = nlinfit(allSugars, Lagdata, lagModel,lagCoeff);
        [yFit, delta] = nlpredci(lagModel,xFit,beta,resid,'Covar',Sigma,'Alpha',alpha,'PredOpt','curve');        
        
        P(1) = plotshaded(xFit,[yFit-delta ; yFit+delta],'k');
        hold on
        set(P(1),'FaceColor',[1 1 1]*0.8);
        P(2) = plot(xFit,yFit,'-k');
        
        hold on
        axis tight
        yl = ylim;
        ylim([0 1.05*yl(2)])
        xlim([0 253])
        xlabel('malto^{3} concentration (\mug per ml)')
        ylabel('estimated lag (h)')
        
        for j = 1:length(sugars)
            s = (j-1)/(length(sugars) - 1);
            thisColour = colourFunction(s);
            F = sugars(j) == allSugars;
            if j == 1
                P(3) = errorbar(sugars(j),mean(Lagdata(F)),ste(Lagdata(F)),'.','markersize',30,'color',thisColour,'linewidth',1);
            else
                errorbar(sugars(j),mean(Lagdata(F)),ste(Lagdata(F)),'.','markersize',30,'color',thisColour,'linewidth',1);
            end
        end
        lagR2 = lagMdl.Rsquared.Adjusted;
        
        se = lagMdl.Coefficients.SE(3);
        
        S1 = ['power regression (adj R^2 \approx ',num2str(lagR2,2),', power\approx ',num2str(lagCoeff(3),2),'\pm ',num2str(se,2),')'];
        S2 = [label,' data'];
        legend(P,{'95% CIs',S1,S2},'location','NorthWest');
        legend('boxoff')
        box on
        %}
        
    end
    warning('on')
end
