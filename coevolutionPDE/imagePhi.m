n = 30;

x = (0:n-1)/(n-1);
y = x;

sigma = 0.01;

phi = @(y,mesh)(exp(-(y.^2)/(2*sigma)));

%phistar = @(y)(erf(2*y)+1)/2;
%phi = @(y,mesh)max(phistar(mesh))-phistar(y);

[X,Y] = meshgrid(1-x,y);

Phi = phi(X-Y,x);

surf(Phi);
view(2)

colormap(hot);
xlabel('phage');
ylabel('bacteria');

axis tight
box on

colorbar


export_fig('./figures/imagePhi.pdf')

