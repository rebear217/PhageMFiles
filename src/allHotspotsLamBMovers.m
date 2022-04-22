function allHotspotsLamBMovers(data)
        
    pNorm = 2;
    blobSize = 20;
    n = length(data.bacteria);
    
    figure(1)
    set(1,'pos',[65           1        1850        1344])
    
    sp = 1;
    for j = 1:n
        try
            subplot(6,7,sp);
            plotMutantDistancesOnLamB(data.bacteria{j},data,pNorm,blobSize);
            title(data.bacteria{j});
            if j == n
                axis tight
                box on
            else
                axis tight
                box on
                axis off
            end
            drawnow
            sp = sp + 1;
        catch
            disp([data.bacteria{j},' has an issue - protein too short?']);
        end
    end
    
    figure(1)
    export_fig('./figures/allHotspotsLamBMovers.pdf')
    
end
