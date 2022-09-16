function plotGrowthRatesFor(strains,data,allRKLagData)

    [~,~,~,~,WTfigure2Data] = plotRYTOfromLabel_FAST('wt',data,allRKLagData,1);
    
    figure(1)
    export_fig('./figures/unused/RYTOplotForWT.pdf')
    figure(2)
    export_fig('./figures/GrowthRatePlotForWT.pdf')
    figure(3)
    export_fig('./figures/unused/RY_RYTOplotForWT.pdf')
    figure(4)
    export_fig('./figures/unused/LagTOplotForWT.pdf')
    
    close all

    for j = 1:length(strains)
        close all
        try
            plotRYTOfromLabel_FAST(strains{j},data,allRKLagData,1);

            figure(1)
            export_fig(['./figures/unused/RYTOplotFor',strains{j},'.pdf'])

            figure(2)
            if ~strcmp('wt',strains{j})
                h = plot(WTfigure2Data.fineSugars,WTfigure2Data.RYTOmodelR,'-','linewidth',2,'color',0.5*[1/3 1/3 1]);
                text(200,1,'WT fit');
                set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end
            export_fig(['./figures/GrowthRatePlotFor',strains{j},'.pdf'])

            figure(3)
            export_fig(['./figures/unused/RY_RYTOplotFor',strains{j},'.pdf'])

            figure(4)
            export_fig(['./figures/unused/LagTOplotFor',strains{j},'.pdf'])

        catch
            disp([strains{j},' has faulty data']);
        end
    end
    
end
