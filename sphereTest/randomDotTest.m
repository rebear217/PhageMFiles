sphere1 = rand(200,3);
sphere2 = [sphere1(1:50,:) ; sphere1(70:end,:)];

displacement = 0.5*(-0.5 + rand(size(sphere2)));
smallDisplacement = 0.02*(-0.5 + rand(size(sphere2)));

%the entire dataset is shifted a little to make a second structure
sphere2 = sphere2 + smallDisplacement;
%some points are shifted a lot in the second structure
sphere2(1:5) = sphere2(1:5) + displacement(1:5);

c = [1 0.1 0]*0.7;
d = [0.1 0.1 1]*0.7;
e = [1 1 1]*0.6;
f = [0.1 1 0.1]*0.6;

%%

figure(1)
subplot(1,2,2)
plot(-100,-100,'o','markersize',10,'color',c)
hold on
plot(-100,-100,'.','markersize',18,'color',d)
plot(-100,-100,'-k')
plot(-100,-100,'o','markersize',18,'color',f)
legend('first protein','second protein','matched alpha carbons','unmatched alpha carbons')

tic
[bestMatchedDistance,bestTotalDistance,s1,s2,code] = rectangularMatching(sphere1,sphere2,10,0,0.035);
toc

subplot(1,2,1)
title('The matching algorithm applied to two fictitious protein structures')
xlim([-0.9 0.9])
ylim([-0.9 0.9])
zlim([-1 1])

subplot(1,2,2)
xlim([-0.9 0.9])
ylim([-0.9 0.9])
zlim([-1 1])

%%

export_fig('dotsMatching.pdf')

figure(2)
plot3(sphere1(:,1),sphere1(:,2),sphere1(:,3),'o','markersize',10,'color',c);
hold on
plot3(sphere2(:,1),sphere2(:,2),sphere2(:,3),'.','markersize',18,'color',d);
view(19.5,14)
box on
title('Two fictitious protein structures')
legend('random dots','same dots with points removed','location','northwest')
xlabel('x')
ylabel('y')
zlabel('z')

export_fig('dotsMatchingPoints.pdf')

