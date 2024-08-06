function [S,D] = OptimizeSDGCA(L,ML,CL)
    [N,~] = size(L);
    maxIter = 300;
    I = eye(N);
    mu1 = 1;
    mu2 = 1;

    S = zeros(N,N);
    D = S;
    Y1 = S;
    Y2 = S;
    F1 = S;
    F2 = S;
    
    for i=1:maxIter
        St = S;
        S = (2*mu1*F1 - D' - Y1) / (2*L + 2*mu1*I);
            
        F1 = Y1 / (2*mu1) + S;
        F1(ML > 0) = 0;
        F1 = min(max(F1, 0),1);
        F1 = F1 + ML;
        F1 = (F1 + F1') / 2;

        Dt = D;
        D = (2*mu2*F2 - S' - Y2) / (2*L + 2*mu2*I);

        F2 = Y2 / (2*mu2) + D;
        F2(CL > 0) = 0;
        F2 = min(max(F2, 0), 1);
        F2 = F2 + CL;
        F2 = (F2 + F2') / 2;

        Y1 = Y1 + mu1*(S-F1);
        Y2 = Y2 + mu2*(D-F2);
        
        rho1 = 1.1;
        rho2 = 1.1;
        mu1 = min(mu1*rho1, 1e6);
        mu2 = min(mu2*rho2, 1e6);

        e = [norm(S-St,'fro'), norm(D-Dt,'fro'), norm(S-F1,'fro'), norm(D-F2,'fro')];

        if max(e) < 1e-3
            break;
        end
    end
end

