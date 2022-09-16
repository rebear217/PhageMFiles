function S = mySmooth(S,n)
    if nargin < 2
        n = 1;
    end
	for j = 1:n
        S = (S(1:end-2) + S(3:end))/2;
        S = [2*S(1) - S(2); S; 2*S(end) - S(end-1)];
    end
end