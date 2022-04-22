function distance = totalEuclidean(p,q)

    separation = abs(p-q).^2;
    distance = sqrt(sum(separation,2));
    
end