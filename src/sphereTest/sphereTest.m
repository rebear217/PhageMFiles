theta = (0:0.2:2*pi)';
phi = (0:0.2:2*pi)';

x = cos(theta).*sin(phi);
y = sin(theta).*sin(phi);
z = cos(phi);

sphere1 = [x y z];

sphere2 = 0.25*rand(size(sphere1)) + sphere1 + ones(size(sphere1))*0.2;
sphere2 = [1.1 1.1 1.01 ; sphere2 ; -0.0125 -0.0125 -0.102];

c = [1 0.1 0]*0.7;
d = [0.1 0.1 1]*0.7;
e = [1 1 1]*0.6;
f = [0.1 1 0.1]*0.6;

figure(1)
subplot(1,2,2)
plot(-100,-100,'o','markersize',10,'color',c)
hold on
plot(-100,-100,'.','markersize',18,'color',d)
plot(-100,-100,'-k')
plot(-100,-100,'o','markersize',18,'color',f)
%legend('first protein','second protein','matched alpha carbons','unmatched alpha carbons')

[bestMatchedDistance,bestTotalDistance,s1,s2,code] = rectangularMatching(sphere1,sphere2,100,0,0.035);

subplot(1,2,1)
title('The matching algorithm applied to two fictitious protein structures')
xlim([-0.9 0.9])
ylim([-0.9 0.9])
zlim([-1 1])

subplot(1,2,2)
xlim([-0.9 0.9])
ylim([-0.9 0.9])
zlim([-1 1])

export_fig('../../figures/exampleMatching.pdf')

figure(2)
plot3(sphere1(:,1),sphere1(:,2),sphere1(:,3),'o','markersize',10,'color',c);
hold on
plot3(sphere2(:,1),sphere2(:,2),sphere2(:,3),'.','markersize',18,'color',d);
view(19.5,14)
box on
title('Two fictitious protein structures')
legend('a figure ''8'' ','a noisy figure ''8'' with two outliers','location','northwest')
xlabel('x')
ylabel('y')
zlabel('z')

export_fig('../../figures/exampleMatchingPoints.pdf')
