function D = computeD(bcs, baseClsSegs)
    [N,M] = size(bcs);
    tau = 0.8;
    D = zeros(N,N);
    simOfCluster = full(simxjac(baseClsSegs));
    RWofCluster = RandomWalkofCluster(simOfCluster);
    RWofCluster = RWofCluster-diag(diag(RWofCluster))+diag(ones(size(RWofCluster,1),1));
    disOfCluster = 1 - RWofCluster;
    for m = 1:M
        D = D + disOfCluster(bcs(:,m),bcs(:,m));
    end
    D = D / M;
    D(D < tau) = 0;
end

function s = simxjac(a,b)

    if ~exist('b', 'var')
      b = a;
    end
    
    n = size(a,1);
    m = size(b,1);
    d = size(a,2);
    if (d~=size(b,2))
      disp('simxjac: data dimensions do not match');
      return;
    end
    
    temp = a * 1 * b';
    asquared = sum((a.^2),2);
    bsquared = sum((b.^2),2);
    s = temp ./ ((asquared * ones(1,m)) + (ones(n,1) * bsquared') - temp);
end

function R = RandomWalkofCluster(W)

    N = size(W,1);
    k = 20;
    beta = 1;

    W = W - diag(diag(W));
    D = diag(1 ./ sum(W, 1));
    D(isinf(D)) = 0;
    W_tilde = D * W;
    
    tmpO = W_tilde;
    O_tilde = W_tilde*W_tilde';
    
    for i = 1:k-1
        tmpO = tmpO*W_tilde;
        O_tilde = O_tilde + beta * (tmpO*tmpO');
    end
    O_i = repmat(diag(O_tilde), 1, N);
    R = O_tilde./sqrt(O_i.*O_i');
    
    sum_W = sum(W_tilde, 2);
    isolatedIdx = find(sum_W<10e-10);
    if numel(isolatedIdx)>0
        R(isolatedIdx,:) = 0;
        R(:,isolatedIdx) = 0;
    end
end