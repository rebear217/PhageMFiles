function [coefficients,se,AIC,RKdata,figure2Data] = plotRYTOfromLabel_FAST(label,data,allRKLagdata,plotFlag)
    warning('off')
    if nargin < 4
        plotFlag = 0;
    end
    
    red = 0.5*[1 0 0];
    %colourFunction = @(s)[s 0 1-s];
    colourFunction = @(s)[0 0 0];
    
    I = getLocationLabels(label,data.bacteria);
    RKdata = allRKLagdata{I};    
    %RKdata = getRKdataFromLabel(label,data,[],1);
    
    [validMaxRate , validMaxYield, validMaxCoV , keepRsquaredValue] = validVariables();
	T = length(RKdata.Rdata);
    for j = 1:5
        switch j
            case 1
                F = RKdata.Rdata < validMaxRate;
            case 2
                F = RKdata.Kdata < validMaxYield*RKdata.allSugars;
            case 3
                F = RKdata.Rsquareddata > keepRsquaredValue;
            case 4
                F = RKdata.Yieldsedata ./ RKdata.Yielddata < validMaxCoV;
            case 5
                F = RKdata.Rsedata ./ RKdata.Rdata < validMaxCoV;                
        end
        
        RKdata.Rsquareddata = RKdata.Rsquareddata(F);
        RKdata.Rdata = RKdata.Rdata(F);        
        RKdata.growthRateData = RKdata.growthRateData(F);
        RKdata.growthRateSEData = RKdata.growthRateSEData(F);
        RKdata.Kdata = RKdata.Kdata(F);
        RKdata.Rsedata = RKdata.Rsedata(F);
        RKdata.Ksedata = RKdata.Ksedata(F);
        RKdata.allSugars = RKdata.allSugars(F);
        RKdata.Yielddata = RKdata.Yielddata(F);
        RKdata.Yieldsedata = RKdata.Yieldsedata(F);
        RKdata.Lagdata = RKdata.Lagdata(F);
        RKdata.Lagdata2 = RKdata.Lagdata2(F);
        RKdata.Lagdata3 = RKdata.Lagdata3(F);
    end
    L = length(RKdata.Rdata);
    disp([label,' data QC: kept ',num2str(L),' out of ',num2str(T),' data points']);
    
    if L > 1
        Rsquareddata = RKdata.Rsquareddata;
        Rdata = RKdata.Rdata;
        Kdata = RKdata.Kdata;
        Rsedata = RKdata.Rsedata;
        Ksedata = RKdata.Ksedata;
        Lagdata = RKdata.Lagdata;
        %Lagdata = RKdata.Lagdata2;
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
        RYTOmdl = NonLinearModel.fit(allSugars,Yielddata,RYTOmodelY,[cHiGuess cLoGuess 1]);
        betaRYTO = abs(RYTOmdl.Coefficients.Estimate);
        seRYTO = RYTOmdl.Coefficients.SE;

        delta = (max(sugars)-min(sugars))/1000;
        fineSugars = [1:delta:(1.5*max(sugars))];

        RYTORsquared = RYTOmdl.Rsquared.Adjusted;

        linearmodel = @(b,x)(b(1)*ones(size(x)));
        linearmdl = NonLinearModel.fit(allSugars,Yielddata,linearmodel,[cLoGuess]);
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
            [betaRY2,resid,J,Sigma] = nlinfit(allSugars, Yielddata, RYTOmodelY,betaRYTO);
            xFit=fineSugars;
            alpha = 0.05;
            [yFit, delta] = nlpredci(RYTOmodelY,xFit,betaRY2,resid,'Covar',Sigma,'Alpha',alpha,'PredOpt','curve');    

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
            %legend('boxoff')

            figure(2)
            P = [];
            plot(allSugars,Rdata,'.','markersize',10,'color',0.75*[1 1 1]);
            hold on

            MonodModel = @(b,S)((abs(b(1))*S./(1+abs(b(2))*S)));
            Monodmdl = NonLinearModel.fit(allSugars,Rdata,MonodModel,[max(Rdata) 1]);
            betaMonod = abs(Monodmdl.Coefficients.Estimate);    

            RYTOmodelR = @(b,S)((abs(b(1))*S./(1+abs(b(2))*S)).*(1+abs(b(3))*S)./(1+(abs(b(3))+abs(b(4)))*S));
            RYTOmdlR = NonLinearModel.fit(allSugars,Rdata,RYTOmodelR,[betaMonod(1) betaMonod(2) 1 0]);
            beta = abs(RYTOmdlR.Coefficients.Estimate);
            se = RYTOmdlR.Coefficients.SE;

            MonodRsquared = Monodmdl.Rsquared.Adjusted;
            RYTORsquared = RYTOmdlR.Rsquared.Adjusted;

            AICMonodmodel= Monodmdl.ModelCriterion.AICc;
            AICRYTOmodel= RYTOmdlR.ModelCriterion.AICc;
            RL = exp(-abs(AICRYTOmodel - AICMonodmodel)/2);

            [betaRY3,resid,J,Sigma] = nlinfit(allSugars, Rdata, MonodModel,betaMonod);
            [yFit, delta] = nlpredci(MonodModel,xFit,betaRY3,resid,'Covar',Sigma,'Alpha',alpha,'PredOpt','curve');        

            P(3) = plotshaded(xFit,[yFit-delta ; yFit+delta],'k');
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
            %legend('boxoff')

            figure2Data.allSugars = allSugars;
            figure2Data.Rdata = Rdata;
            figure2Data.fineSugars = fineSugars;
            figure2Data.RYTOmodelR = RYTOmodelR(beta,fineSugars);

            figure(3)
            P = [];
            plot(Rdata,Yielddata,'.','markersize',10,'color',0.75*[1 1 1]);
            hold on

            for j = 1:length(allSugars)
                s = (las(j)-minas)/(maxas - minas);
                thisColour = colourFunction(s);
                if j == length(allSugars)
                    P(1) = plot(Rdata(j),Yielddata(j),'Marker','.','MarkerSize',40,'color',thisColour);
                else
                    plot(Rdata(j),Yielddata(j),'Marker','.','MarkerSize',40,'color',thisColour);
                end
            end        
            RY = RYTOmodelR(beta,fineSugars);
            P(2) = plot(RY,RYTOmodelY(betaRYTO,fineSugars),'-','linewidth',2,'color',0.25*[1 1 1]);        
            axis tight
            ym = ylim; ym = max(Yielddata);
            ylim([0 1.6*ym]);
            xlim([0 max(RY)])
            xlabel('growth rate (per h)')
            ylabel('yield (OD_{600} per malto^{3})')
            legend(P,{[label,' data (blue is low sugar, red is high sugar)'],'theoretical RYTO'},'location','NorthWest');
            %legend('boxoff')

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
            %legend('boxoff')
            box on


        end
    else
        disp('stopped analysis, not enough high quality data');
        
        coefficients = NaN;
        se = NaN;
        AIC = NaN;
        RKdata = NaN;
        figure2Data = NaN;
    end
    warning('on')
end
