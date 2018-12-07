function [R,t] = matchPoints3d(p,q,w)
    [n,m] = size(p);
    if not(m == 3)
        error('Data must be n x 3 in size for this to work')
    end
    if nargin < 3
        w = ones(n,1);
    end
    W = repmat(w,1,3);
    pBar = sum(W.*p,1)/sum(w);
    qBar = sum(W.*q,1)/sum(w);
    
    x = p - repmat(pBar,n,1);
    y = q - repmat(qBar,n,1);
    
    X = x';
    Y = y';
    
    S = X*diag(w)*Y';
    [U,~,V] = svd(S);
    R = V*U';
    t = qBar' - R*pBar';
end