function NWCA=computeNWCA(baseClsSegs,NECI,M)
    baseClsSegs = baseClsSegs';
    N = size(baseClsSegs,1);
    
    NWCA = (bsxfun(@times, baseClsSegs, NECI')) * baseClsSegs' / M;
    NWCA = NWCA / max(max(NWCA));
    NWCA = NWCA-diag(diag(NWCA))+eye(N);
end