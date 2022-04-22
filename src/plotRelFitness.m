function plotRelFitness(lamBstrains,relativeFitnesses,data)

    WTRF = relativeFitnesses(end-2:end,:);
    wtRF = mean(WTRF(:));
    wtSE = STE(WTRF);

    lamBstrains{end} = 'wt';
    lamBstrains{end-1} = 'wt';
    lamBstrains{end-2} = 'wt';

    %remove 28a:
    B = data.infection.bacteria_names;
    Ij = getLocationLabels('28a',data.infection.bacteria_names);
    B1 = {data.infection.bacteria_names{1:Ij-1}};
    B2 = {data.infection.bacteria_names{Ij+1:end}};
    B = {B1{:},B2{:}};
    data.infection.bacteria_names = B;
    
    [strains,li,di] = intersect({lamBstrains{:}},{data.infection.bacteria_names{:}});

    lamBstrains = {strains{1:end-1},'wt','wt','wt'};
    relativeFitnesses = [relativeFitnesses(li(1:end-1),:) ; WTRF];

    figure(1)
    set(1,'pos',[42         394        1399         412]);
    
    N = length(lamBstrains);
    col = @(s)[s 0.25 1-s];

    Ste = ste(relativeFitnesses);
    Mean = mean(relativeFitnesses,2);
        
    [~,K] = sort(Mean);
    PhiMap = zeros(1,N);
    Col = zeros(N,3);
    
    for J = 1:N
        
        X = getLocationLabels(lamBstrains{J},data.infection.bacteria_names);
        if ~isempty(X)
            PhiMap(J) = X;
        end
        
        s = (J-1)/(N-1);
        j = K(J);

        plot([J J],[Mean(j) - 1.96*Ste(j),Mean(j) + 1.96*Ste(j)],'color',[1 1 1]*0.2,'linewidth',1);
        hold on
        p = plot(J,Mean(j),'.','markersize',30);
                
        Col(J,:) = col(s);
        set(p,'color',Col(J,:));
    end
    
    lB = lamBstrains(K);
    ylim([0 1.2]);
    xlim([1 N]);
    set(gca,'Xtick',1:N);
    set(gca,'XtickLabel',lB,'FontSize',9);
    %set(gca,'FontSize',9);
    ylabel('relative fitness','FontSize',16);
    xlabel('E.coli strain','FontSize',16)
    
    %rotateXLabels( gca(), 45 );
    %xtickangle(45);
    
    axes('Position',[0.3138      0.23301      0.57327       0.5267]);
    set(gca,'FontSize',11);
    for J = 1:N
        X = PhiMap(J);
        j = K(J);
        if X > 0
            phi = data.infection.phi(:,X);
            plot(myLog([mean(phi) mean(phi)]),[Mean(j) - 1.96*Ste(j),Mean(j) + 1.96*Ste(j)],...
                '-','color',[0.2 0.2 0.2],'linewidth',1);
            hold on
            plot(myLog(mean(phi)),Mean(j),'.','markersize',30,'color',Col(J,:));
        end
    end

    axis tight
    xlim([-3.1 0.1])
    ylim([0.65 1.175])
    xl = xlim;
    plot([xl(1) xl(2)],[wtRF wtRF],':','color',[0.6 0.2 0.2],'linewidth',1);
    set(gca,'Xtick',[-3 -2 -1 0]);
    set(gca,'Xticklabel',{'-inf','-2','-1','0'});
    xlabel('log (mean phage susceptibility)','FontSize',14);
    ylabel('relative fitness','FontSize',14);
    
    export_fig('./figures/relativeFitness.pdf')
    
end

function s = ste(v)
	[~,n] = size(v);
	s = std(v,[],2)/sqrt(n-1);
end

function s = STE(v)
	v = v(:);
    n = length(v);
	s = std(v)/sqrt(n-1);
end

function I = getLocationLabels(name,nameSet)
    I = find(cellfun(@(s)strcmpi(s,name),nameSet));
end

function y = myLog(x)
    %this is for inset plotting purposes:
    if x > 0
        y = log10(x);
    else
        y = log10(1e-3*ones(size(x)));
    end
end
