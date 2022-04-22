function plotOneLagExample(data,label,sugar)

    if nargin < 2
        label = 'wt';
        sugar = 80;
    end
    %labels = {labels{1}};
    
	sugars = unique(data.sugars{1});
    sugars = sugars(~isnan(sugars));
    
    %figure(j)
    %set(j,'pos',[259   878   573   438]);

    k = find(sugars == sugar);
    if isempty(k)
        k = 4;
    end
    
    [FitData,P,Pm] = logisticFitToODfromLabel(sugars(k),data,label,2);
    R2 = min(FitData.Rsquared);
    %hold on
    
    %delete(P(5:end));
    %P = [P(1) P(2) P(3) P(4)];
    
    %for J = 1:12
    %    unplot;
    %end
    
    s = (k-1)/(length(sugars) - 1);
    clr = [1 1 1]/2;

    %axis tight
    ylim([-0.001 0.1]);
    xlim([0 24]);
    
    for J = 1:length(P)
        set(P(J),'color',clr);
        set(P(J),'linewidth',3);
    end

    %spaceplots;        
    disp(['bacterium ',label]);
    %text(0.1,0.103,'population density')
    %annotation(gcf,'arrow',[0.58 0.51],[0.54 0.68]);

    %KP = plot([0 24],[M.K(1) M.K(1)],'-','color',[0 0 0],'linewidth',1);
    %plot([0 24],[M.blank(1) M.blank(1)],':','color',[1 1 1]*0.5,'linewidth',1);

    legend([P(1),Pm(1)],{[label,' OD @ ',data.sugarString2{k+1}],['logistic fit (R^2 \approx ',num2str(R2,2),')']});
    %legend('boxoff');
    %drawnow
    
    %export_fig('./figures/plotOneLagExample.pdf');
    %axis tight

end