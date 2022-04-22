function p=plotAProtein1(p,col,size)    
	p=plot3(p(:,1),p(:,2),p(:,3),'o','markersize',6+15*size,'color',col,'linewidth',1);
    hold on
    box on
    xlabel('x')
    ylabel('y')
    zlabel('z')
end