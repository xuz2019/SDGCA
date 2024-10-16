function [bcs, baseClsSegs] = getAllSegs(baseCls)
    % n: the number of data points.
    % M: the number of base clusterings.
    % nCls: the number of clusters (in all base clusterings).
    [N,M] = size(baseCls);
    bcs = baseCls;
    nClsOrig = max(bcs,[],1);
    C = cumsum(nClsOrig);
    bcs = bsxfun(@plus, bcs,[0 C(1:end-1)]);
    % nCls = nClsOrig(end)+C(end-1);
    nCls = C(end);
    baseClsSegs=sparse(bcs(:),repmat((1:N)',1,M), 1,nCls,N);
    baseClsSegs = full(baseClsSegs);
end