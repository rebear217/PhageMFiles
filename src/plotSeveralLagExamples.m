function plotSeveralLagExamples(data,labels)

    close all
    if nargin < 2
        labels = {'wt'};
    end
    
	sugars = unique(data.sugars{1});
    sugars = sugars(~isnan(sugars));
            
    for j = 1:length(labels)
        
        disp('closing all figures in plotSeveralLagExamples')
        
        close all;
        
        figure(1);
        set(1,'pos',[58   722   834   616]);
        label = labels{j};
                
        Plots = zeros(length(sugars)+1,1);
        R2 = zeros(length(sugars),1);
        L = {};
        
        for k = 1:length(sugars)
            %subplot(1,length(sugars),k);
            [FitData,P,Pm] = logisticFitToODfromLabel(sugars(k),data,label,2);
            R2(k) = mean(FitData.Rsquared);
            delete(P(2:end));
            delete(Pm(2:end));
            
            P = P(1);
            Pm = Pm(1);
            
            Plots(k) = P(1);            
            
            s = (k-1)/(length(sugars) - 1);
            clr = [1-s 0 s];            
            
            ylim([0.0 0.1]);
            xlim([0 24]);
            
            for J = 1:length(P)
                set(P(J),'color',clr);
                set(P(J),'linewidth',4);
                set(Pm(J),'linewidth',1);
                set(Pm(J),'color',[0 0 0]);
            end
            
            drawnow
            
            if k == 1
                L{k} = [data.sugarString2{k+1},' ',label,' OD data (R^2 \approx ',num2str(R2(k),2),')'];
            else
                L{k} = [data.sugarString2{k+1},' (R^2 \approx ',num2str(R2(k),2),')'];
            end
        end
        
        Plots(k+1) = Pm;
        L{k+1} = 'logistic datafits';
        
        legend(Plots,L);
        legend('location','northwest','FontSize',24);
        
        disp(['bacterium ',labels{j}]);
      
        export_fig(['./figures/plotOneLagExample',label,'.pdf']);
        
    end


end