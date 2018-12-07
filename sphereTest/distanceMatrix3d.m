function dm = distanceMatrix3d(p,q)
    [n,m] = size(p);
    [N,M] = size(q);
    if not(m == 3) && (M == 3)
        error('Data must be N x 3 in size for this to work')
    end
    Onn = ones(n,1);
    OnN = ones(N,1);
    dx = p(:,1)*OnN' - Onn*q(:,1)';
    dy = p(:,2)*OnN' - Onn*q(:,2)';
    dz = p(:,3)*OnN' - Onn*q(:,3)';
    
    dm = sqrt(dx.^2 + dy.^2 + dz.^2);
end