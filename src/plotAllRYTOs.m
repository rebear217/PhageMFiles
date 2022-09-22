function plotAllRYTOs(data,allRKLagData)

    figure(1)
    set(1,'pos',[1           1        2490        1344])

    yieldIssueList = {'13a','19a','28b','94a'};
    sugarsList = [250,160,250,250];
    sL = 1;
    outGuesses = [];
    
    for j = 1:46
        subplot(6,8,j);
        if ~ismember(data.bacteria{j},yieldIssueList)
            if j <= length(data.bacteria)
                if strcmp(data.bacteria{j},'70b')
                    [~,~,~,~,~,outGuesses] = plotJustRYTOfromLabel_FAST(data.bacteria{j},allRKLagData,1,outGuesses);
                else
                    [~,~,~,~,~,outGuesses] = plotJustRYTOfromLabel_FAST(data.bacteria{j},allRKLagData,1);
                end
                ylabel('yield','FontSize',20)
                xlabel('growth rate (per h)','FontSize',18)

                xlim([5e-3 2])
                ylim([0 1.75e-3])
            else
                box on
            end
        else
            plotOneLagExample(data,data.bacteria{j},sugarsList(sL));
            legend('location','west')
            sL = sL + 1;
        end

        if ismember(j,[1 1+8 1+2*8 1+3*8 1+4*8 1+5*8])
            drawnow
        end

    end

    export_fig('./figures/allRYTOplots.pdf')

end