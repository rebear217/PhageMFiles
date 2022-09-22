function plotInteractionMatrix(data)

    bacteria = data.infection.bacteria_names;
    phage = data.infection.phage_names;
    phi = data.infection.phi';
    
    %28a and 28b are the same, no need to plot 2x :
    removeBacList = {'28a','wt(first)'};
    for i = 1:length(removeBacList)
        f = find(strcmp(bacteria,removeBacList{i}));
        if f > 1 && f < length(bacteria)
            phi = phi([1:f-1 f+1:end],:);
            B1 = {bacteria{1:f-1}};
            B2 = {bacteria{f+1:end}};
            B = {B1{:},B2{:}};
        end
        if f == 1
            phi = phi([2:end],:);
            B = {bacteria{2:end}};
        end
        if f == length(bacteria)
            phi = phi([1:f-1],:);
            B = {bacteria{1:f-1}};
        end
        bacteria = B;
    end

    removePhaList = {'blank','blank','blank'};
    for i = 1:length(removePhaList)
        f = find(strcmp(phage,removePhaList{i}));
        f = f(1);
        if f > 1 && f < length(phage)
            phi = phi(:,[1:f-1 f+1:end]);
            P1 = {phage{1:f-1}};
            P2 = {phage{f+1:end}};
            P = {P1{:},P2{:}};
        end
        if f == 1
            phi = phi(:,[2:end]);
            P = {phage{2:end}};
        end
        if f == length(bacteria)
            phi = phi(:,[1:f-1]);
            P = {phage{1:f-1}};
        end
        phage = P;        
    end
    
    N = 10;
    col = zeros(N,3);
    for j = 1:N
        s = (j-1)/N;
        c = [1 1 1]*(1-s) + s*[0 1 0.8]*0.8;
        col(j,:) = c;
    end

    close all
    clc

    %phi = PHI;
    bacteria_names = B;
    phage_names = P;

    phi = phi/max(max(phi));
    
    [~,J] = sort(sum(phi,2));
    [~,I] = sort(sum(phi,1));
    bacteria_names = bacteria_names(J);
    phage_names = phage_names(I);
    phi = phi(J,I);
    
    figure(1)
    set(1,'color','white')
    set(1,'position',[59         168        1622         961])
    imagesc(phi);
        
    colormap(col);
    colorbar

    myLabels();
    hold on
    plot([0 94],[12.5 12.5],'-k','linewidth',1);
    text(1,11,'pan resistant','fontsize',22)
    text(1,14,'susceptible to phage','fontsize',22)
    
	export_fig('figures/infectionMatrix-imageSc.pdf')
	%exportgraphics(gca,'figures/infectionMatrix-imageSc.png','Resolution',1000)
    %export_fig('figures/infectionMatrix-imageSc.tif')
    
    function myLabels()
        box on
        axis([-0.0 1+length(phage_names) -0.25 length(bacteria_names)+1]);
        set(gca,'YTick', -0.5 + 1.5:1:(length(bacteria_names)+0.5))
        set(gca,'YTickLabel',bacteria_names,'fontsize',17)   
        set(gca,'TickLength', [0 0]);
        set(gca,'XTick',(1.5:1:(length(phage_names)+0.5)) - 0.5)
        set(gca,'XTickLabel',phage_names,'fontsize',17)  
        xticklabel_rotate
        %xtickangle(-45)
        xlabel('Phage','fontsize',24);
        ylabel('Bacteria','fontsize',24);
    end
    
end