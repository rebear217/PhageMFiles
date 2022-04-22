function P = plotAProteinColourVector(p,colourVector,blobSize)
    if nargin < 3
        blobSize = 60;
    end
    
    if colourVector == 0
        colourVector = zeros(size(p));
        colourVector = colourVector(:,1);
        useBlack = 1;
    else
        colourVector = (colourVector - min(colourVector))/(max(colourVector) - min(colourVector));
        useBlack = 0;
    end
    
    grey = [1 1 1]*0.95;
    
    P = plot3(p(:,1),p(:,2),p(:,3),'-','linewidth',0.5,'color',grey);
    hold on
    box off
    xlabel('x')
    ylabel('y')
    zlabel('z')
	view(36,18)
    
    for j = 1:length(colourVector)
        if useBlack
            col = [0 0 0];
        else
            s = colourVector(j);
            col = s*[1 0 0] + (1-s)*[1 1 1];
        end
        plot3(p(j,1),p(j,2),p(j,3),'.','markersize',blobSize,'color',col);
    end
    
end