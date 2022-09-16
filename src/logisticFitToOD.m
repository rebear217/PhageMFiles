function [fitData,P,Pm] = logisticFitToOD(T,OD,plotFlag)
    %T is length of time in minutes
    warning('off');
    n = length(OD);
    
    x = (T/60)*(0:(n-1))/(n-1);
    xpred = 24*(0:(n-1))/(n-1);
    
    y = OD;
    
    Coefficients = NaN;
    SE = NaN;
    flag = -1;
        
    DeltaOD = max(OD)-min(OD);
    rguess = 1;
    
    beta0 = [DeltaOD 10000 rguess mean(OD(1:5))];
    modelfun = @(b,X)(abs(b(4)) + abs(b(1))./(1 + abs(b(2))*exp(-X*abs(b(3)))));
	modelfunDeriv = @(b,X)(abs(b(1)).*(abs(b(2))*abs(b(3))*exp(-X*abs(b(3))))./(1 + abs(b(2))*exp(-X*abs(b(3)))).^2);

    %modelfunR = @(b,X)(abs(b(4)) + abs(b(1))*exp(X*abs(b(3)))./(abs(b(2))));
    
    
    if plotFlag
        mdl = NonLinearModel.fit(x,y,modelfun,beta0)
    else
        mdl = NonLinearModel.fit(x,y,modelfun,beta0);
    end
    Coefficients = abs(mdl.Coefficients.Estimate);
    SE = mdl.Coefficients.SE;
    blank = Coefficients(4);
    
    P = -1;
    Pm = -1;
    
    if plotFlag
        P = plot(x,OD-blank,'-','linewidth',1,'color',[1 1 1]*0.25);
        hold on
        Pm = plot(xpred,modelfun(Coefficients,xpred)-blank,'-','linewidth',2,'color',[1 1 1]*0.5);
        %plot(xpred,60*modelfunDeriv(Coefficients,xpred),'-','linewidth',1,'color',[1 1/3 1/3]*0.75);
        axis tight
        xlabel('time (h)','fontsize',38);
        ylabel('OD_{600}','fontsize',38);
        set(gca,'Fontsize',34)
        %legend('OD_600','fit','location','northwest');
    end
    
    fitData.Rsquared = mdl.Rsquared.Ordinary;
    fitData.logisticCoefficients = Coefficients;
    fitData.logisticCoefficientsSE = SE;
    
    fitData.K = Coefficients(1);
    fitData.Kse = SE(1);

    RmaxPerH = 60*Coefficients(1)*Coefficients(3)/4;
    dRmaxPerH = 60*(SE(1)*Coefficients(3) + Coefficients(1)*SE(3))/4;
    
    fitData.R = RmaxPerH;
    fitData.Rse = dRmaxPerH;
    
    %I think this is wrong, it's not the doubling time...maybe...
    fitData.doublingTimeInhours = log(2)/RmaxPerH;
    fitData.doublingTimeInhoursSE = dRmaxPerH*log(2)/RmaxPerH^2;
    
    fitData.blank = blank;
    fitData.blankse = SE(4);
    
    fitData.innoculum = 1/(1+Coefficients(2));    
    lagTimeMeasure = log(Coefficients(2))/Coefficients(3);
    fitData.lagTimeMeasure1 = lagTimeMeasure;
        
    dydtAtLag = modelfunDeriv(Coefficients,lagTimeMeasure);
    yAtLag = modelfun(Coefficients,lagTimeMeasure);
    %solve yAtLag + (t-lagTimeMeasure)*dydtAtLag = blank
    %linearLagFun = @(t)(yAtLag + (t-lagTimeMeasure)*dydtAtLag);
    %fsolve(@(x)(linearLagFun(x)-blank),lagTimeMeasure)
    t = lagTimeMeasure+(blank-yAtLag)/dydtAtLag;
    fitData.lagTimeMeasure2 = t;

    %lagPointData = getFirstNonlinearPoint(OD);
    lagPointData = getFirstNonlinearPoint2(OD);
    if ~isnan(lagPointData.j)
        fitData.lagTimeMeasure3 = x(lagPointData.j);
    else
        fitData.lagTimeMeasure3 = NaN;
    end
    
    fitData.lagTimeMeasure4 = NaN;
    
    warning('off');
    
    if plotFlag == 2
        t = fitData.lagTimeMeasure2;
        %plot(t,0.0,'.k','markersize',10,'color',[1 1 1]*0.5);
        %plot([t t],[0 1],':','color',[1 1 1]*0.5,'linewidth',1);
    end
    
    if plotFlag == 3
        plot(xpred(end),fitData.blank+fitData.K,'.k','markersize',30,'color',[1 1 1]*0.5);
        plot(xpred,fitData.blank+ones(size(xpred))*fitData.K,':','color',[1 1 1]*0.5,'linewidth',1);
    end
    
    %if plotFlag == 4
    %    plot(xpred,modelfunR(Coefficients,xpred),':','color',[1 1 1]*0.5,'linewidth',1);
    %end
end