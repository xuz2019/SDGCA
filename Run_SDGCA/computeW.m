function Wstar = computeW(S,D,W)
    S = min(max(S, 0),1);
    S = (S + S') / 2;
    D = min(max(D, 0),1);
    D = (D + D') / 2;
    W1 = 1 - (1 - S + D) .* (1 - W);
    W2 = (1 + S - D) .* W;
    flag = S - D;
    W1(flag < 0) = 0;
    W2(flag >= 0) = 0;
    Wstar = W1 + W2;
end

