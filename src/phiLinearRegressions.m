function phiLinearRegressions(data)
    
    Phi = data.infection.phi';
    binaryPhi = (Phi > 0);

    figure(1)
    set(1,'pos',[55           1        2478        1345])

    [n,~] = size(Phi);
    f = @(b,x)(b(1) + b(2)*x);
    
    sp = 1;
    
	xData = sum(binaryPhi);
    for j = 1:n-1
        subplot(6,7,1)
        yData = Phi(j,:);
        h = (max(xData) - min(xData))/3;
        xFine = min(xData):h:max(xData);

        if sum(yData) > 0
            s = (j-1)/(n-1);
            colour = [0 0 0];
            colour2 = [1 1 1]*0.7;

            %mdl = NonLinearModel.fit(xData,yData,f,[1 1]);
            %b0 = [mdl.Coefficients.Estimate(1) mdl.Coefficients.Estimate(2)];
            
            [ b0 sigma2_x x_est y_est stats] = deming(xData(:),yData(:));
            [rho,pval] = corr(xData(:),yData(:));
                
            if pval > 0.05
                colour2 = [1 0.5 0.5];
            end
                        
            subplot(6,7,sp)
            plot(xData,yData,'.','color',colour2,'markersize',30);
            axis tight
            hold on
            
            leg = [data.infection.bacteria_names{j},...
                ': p \approx ',num2str(pval,3),...
                ', \rho \approx ',num2str(rho,2)];
            
            leg1 = [data.infection.bacteria_names{j},...
                ': p \approx ',num2str(pval,3)];
            leg2 = ['\rho \approx ',num2str(rho,2)];
                        
            %legend(leg,'location','northeast','fontsize',12)
            %legend('boxoff')                
            
            text(1,0.9*max(yData),leg1);
            text(1,0.8*max(yData),leg2);            
            
            p(j) = plot(xFine,f(b0,xFine),'-','linewidth',3,'color',colour);
            hold on
            l{j} = [data.infection.bacteria_names{j},' (p \approx',num2str(pval,2),')'];
                    
            drawnow

            xlabel('# bacteria phage can kill','fontsize',14)
            ylabel(['phage infectivity to ',data.infection.bacteria_names{j}],'fontsize',14)
            
            sp = sp + 1;
        end
    end
    
    %L = legend(p,l);
    %LEG = findobj(L,'type','text');
    %set(LEG,'FontSize',8)

    export_fig('./figures/phiLinearRegressions.pdf')
    
end