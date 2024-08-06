function [f,p,r] = compute_f(H,T)

N = length(T);
numT = 0;
numH = 0;
numI = 0;
for n=1:N
    Tn = (T(n+1:end))==T(n);
    Hn = (H(n+1:end))==H(n);
    numT = numT + sum(Tn);
    numH = numH + sum(Hn);
    numI = numI + sum(Tn .* Hn);
end
p = 1;
r = 1;
if numH > 0
    p = numI / numH;
end
if numT > 0
    r = numI / numT;
end
if (p+r) == 0
    f = 0;
else
    f = 2 * p * r / (p + r);
end
