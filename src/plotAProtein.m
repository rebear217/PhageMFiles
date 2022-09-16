function p=plotAProtein(p,col)
    c = [1 1 1]*0.7;
    if nargin == 1
        col = c;
    end
    
    %plot3(p(:,1),p(:,2),p(:,3),'-','linewidth',0.5,'color',c)
	p=plot3(p(:,1),p(:,2),p(:,3),'.','markersize',30,'color',col);
    hold on
    box on
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title('A singe protein structure')
end