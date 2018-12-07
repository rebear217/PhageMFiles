function endstate = Go(state)

    x = 0.0:0.001:1;
    m = length(x);

    if nargin == 0
        iN = -300*(x-0.25).^2;%-rand(size(x));%
        iM = -300*(x-0.25).^2;%-rand(size(x));%
        inS = 0.5;
    else
        iN = state(1:m);
        iM = state(m+1:end-1);
        inS = state(end);
    end

    initialTime = 500;
    stepTime = 1000;
    [iNe,iMe,inS,fig]=run(initialTime,iN,iM,inS,0);
    timeMin = initialTime;
    figlist = [ fig ];
    while input('1 to go on, 0 to stop ... ')
        [iNe,iMe,inS,fig] = run(stepTime,iNe,iMe,inS,timeMin);
        timeMin = timeMin + stepTime;
        figlist = [figlist ; fig];
    end

    endstate = [iNe iMe inS];
    
    for f = 1:length(figlist)
        figure(figlist(f))
        export_fig(['./figures/coevolutionFigure',num2str(f),'.pdf']);
    end

    F = figure();
    set(F,'pos',[524   267   670   438]);
    p1=plot(x,exp(iNe),'-b');
    hold on
    p2=plot(x,exp(iMe),'-r');
    legend([p1,p2],{'bacteria','phage'})
    legend('boxoff')
    xlabel('recognition phenotype');
	ylabel('densities');
    export_fig(['./figures/coevolutionEndFigure',num2str(F.Number),'.pdf']);
    
end