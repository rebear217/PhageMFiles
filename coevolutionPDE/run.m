function [iNe,iMe,Se,thisFig] = run(Time,iN,iM,inS,timeMin)

    n = length(iN);

    parameters.muN = 1e-8;
    parameters.muM=1e-8;
    parameters.beta = 2;
    parameters.C=1;
    parameters.d=0.1;
    parameters.S0=10;

    % construct the initial condition
    ic = [iN iM inS];
    parameters.x = (0:n-1)/(n-1);
    x = parameters.x;
    parameters.y = (1-n:n-1)/(n-1);

    phageColour = [1 0.2 0.2];
    bacteriaColour = [0.2 0.2 1];

    sigma = 0.001;

    %lock and key (matching alleles) recognition system:
    parameters.phi = exp(-(parameters.y.^2)/(2*sigma));
    parameters.phistar = parameters.phi;
    %graded (gene-for-gene) recognition system:
    %parameters.phistar = (erf(2*parameters.y)+1)/2;
    %parameters.phi = max(parameters.phistar)-parameters.phistar;
    
    %Because convolution has been implemented the wrong way round, phistar is
    %the infection probability function, not phi. phi is the adjoint of phistar.

    %solve the system with the given initial conditions
    %options = odeset('RelTol',1e-5,'AbsTol',1e-4*ones(2*n+1,1));
    options = odeset('RelTol',1e-4,'AbsTol',1e-7,'NormControl','on');

    [T,V] = ode113(@(t,u)RHS(t,u,parameters),[timeMin timeMin+Time],ic);

    % reverse exponential coordinate transformation
    N = exp(V(:,1:n));
    M = exp(V(:,n+1:2*n));
    S = V(:,2*n+1);

    %maxN = findmax(N);
    %maxM = findmax(M);
    %figure(6);
    %subplot(2,1,1);spy(maxN',12);
    %subplot(2,1,2);spy(maxM',12);

    %construct the output of run.m matlab function
    tN = length(T);
    EndPt = V(tN,:);
    iNe = EndPt(1:n);
    iMe = EndPt(n+1:2*n);
    Se = EndPt(2*n+1);

    % plot time vs. resource concentration
    %figure(4);
    %plot(T,S);
    %title('Resource Concentration');
    %xlabel('t');
    %ylabel('S');

    thisFig = figure;
    set(thisFig,'pos',[742         532        1455         813])

    for i = 1:20:tN,
        Nd = N(i,:)/sum(N(i,:));
        Md = M(i,:)/sum(M(i,:));
        mu = dot(x,Nd);
        c = norm(Md-Nd,2);
        nu = sqrt(dot((x-mu).^2,Nd));

        figure(thisFig)

        p1 = plot3(T(i)*ones(size(x)),x,N(i,:),'color',bacteriaColour,'linewidth',1);
        hold on
        p2 = plot3(T(i)*ones(size(x)),x,M(i,:),'color',phageColour','linewidth',1);
        axis tight
        box on
        view(-73.5,58);
        xlabel('time')
        ylabel('recognition phenotype');
        zlabel('densities');
        legend([p1,p2],{'bacteria','phage'});
        legend('boxoff')
        drawnow
    end

end

function F = RHS(t,u,parameters)

    x = parameters.x;
    phi = parameters.phi;
    phistar = parameters.phistar;
    %rho = parameters.rho;
    beta = parameters.beta;
    C = parameters.C;
    muN = parameters.muN;
    muM = parameters.muM;
    d = parameters.d;
    S0 = parameters.S0;

    m = length(u);
    n = (m-1)/2;

    fx=ftoff(x);
    gx=gtoff(x);

    W = ones(1,n);
    W(1) = 1/2;
    W(end) = 1/2;

    N = u(1:n);
    M = u(n+1:2*n);
    S = u(2*n+1);

    Ndot = muN*D(N) + (C*x.*gx*S)./(fx+S) - d - IntOp(phi,exp(M))';
    Mdot = muM*D(M) + beta*IntOp(phistar,exp(N))' - d;
    Sdot = d*(S0-S)-S*sum(W.*x.*exp(N)'./(fx+S))/n;

    F = [Ndot Mdot Sdot]';

end

function y = gtoff(x)
    y=1-x;
    %for no trade-off, use this:
    %y = ones(size(x));
end

function y = ftoff(x)
    y = x;
	%for no trade-off, use this:
    %y = ones(size(x));
end

function O = IntOp(phi,f)
    n = length(f);
    O = Convolve(phi,f)'/(n-1);
end

function Cv = Convolve(K,f)
    %weights for trapezium rule
    N2 = length(f);
    W = ones(N2,1);W(1) = 1/2;W(end) = 1/2;
    f = W.*f;
    F = [zeros(N2-1,1); f]';
    Cu = real(ifft(fft(K).*fft(F)));
    Cv = [Cu(2*N2-1) Cu(1:N2-1)];

    %this is the wrong way round - k-i should reverse roles.
    %the solution to this has been to reverse the roles
    %of phi and phistar

    %for k = 1:N,
    %    Cvv(k) = 0;
    %    for i = 1:N;
    %        Cvv(k) = Cvv(k) + K(k-i+N)*f(i);
    %    end
    %end
    %norm(Cv-Cvv,2)
    %pause

end

function DD = D(Z)

    LP = length(Z);
    %Neuman BCs
    Q(1) = 2*(Z(2)-Z(1));
    Q(LP) = 2*(Z(LP-1)-Z(LP));
    Q(2:LP-1) = Z(1:LP-2)-2*Z(2:LP-1)+Z(3:LP);

    D = [0 diff(Z(1:LP-1))' 0];
    DD = (LP*LP*Q + (LP*D/2).^2);

end
