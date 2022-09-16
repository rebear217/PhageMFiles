function [p1,p2]=plotAProtein2(p,col)
    grey = [1 1 1]*0.8;
    if nargin == 1
        col = grey;
    end
    
    p1 = plot3(p(:,1),p(:,2),p(:,3),'-','linewidth',1,'color',grey);
    hold on    
	p2 = plot3(p(:,1),p(:,2),p(:,3),'.','markersize',40,'color',col);
    box on
    xlabel('x')
    ylabel('y')
    zlabel('z')
	view(36,18)
    plot3(p(239:263,1),p(239:263,2),p(239:263,3),'-','linewidth',4,'color',[1 1 0.5]);%L6
    plot3(p(149:165,1),p(149:165,2),p(149:165,3),'-','linewidth',4,'color',[0 0 1]);%L4
    plot3(p(196:208,1),p(196:208,2),p(196:208,3),'-','linewidth',4,'color',[0.5 0.75 1]);%L5
    plot3(p(165:195,1),p(165:195,2),p(165:195,3),'--','linewidth',4,'color',[0.5 0.75 1]);%L5 connection
    plot3(p(379:401,1),p(379:401,2),p(379:401,3),'-','linewidth',4,'color',[1 0 0]);%L9
    plot3(p(106:123,1),p(106:123,2),p(106:123,3),'-','linewidth',4,'color',[1 1 1]*0.5);%L3
    plot3(p(69:79,1),p(69:79,2),p(69:79,3),'-','linewidth',4,'color',[0 0 0]);%L2
end