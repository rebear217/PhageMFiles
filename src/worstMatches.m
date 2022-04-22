function [indexes,distances] = worstMatches(p,q,k)

    separation = abs(p-q).^2;
    distances = sqrt(sum(separation,2));
    [~,I] = sort(distances,'descend');
    indexes = I(1:k);
    
end