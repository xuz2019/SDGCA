function NECI = computeNECI(bcs, baseClsSegs, para_theta)
    M = size(bcs,2);
    ETs = getAllClsEntropy(bcs, baseClsSegs);
    NECI = exp(-ETs./para_theta./M);
end

function Es = getAllClsEntropy(bcs, baseClsSegs)
    baseClsSegs = baseClsSegs';
    [~, nCls] = size(baseClsSegs);
    Es = zeros(nCls,1);
    normK = computeNormK(bcs);
    for i = 1:nCls
        partBcs = bcs(baseClsSegs(:,i)~=0,:);
        Es(i) = getOneClsEntropy(partBcs, normK(i));
    end
end

function E = getOneClsEntropy(partBcs, normK_i)
    E = 0;
    for i = 1:size(partBcs,2)
        tmp = sort(partBcs(:,i));
        uTmp = unique(tmp);
        
        if numel(uTmp) <= 1
            continue;
        end
        cnts = zeros(size(uTmp));
        for j = 1:numel(uTmp)
            cnts(j)=sum(sum(tmp==uTmp(j)));
        end
        
        cnts = cnts./sum(cnts(:));
        E = E-sum(cnts.*log2(cnts));
    end
    E = E ./ log2(normK_i);
end

function normK = computeNormK(bcs)
    clust = max(bcs, [],1);
    normK = zeros(1, clust(end));
    normK(1:clust(1)) = repmat(clust(1), 1,clust(1));

    for i=1: length(clust)-1
        normK(clust(i)+1:clust(i+1)) = repmat(clust(i+1)-clust(i), 1, clust(i+1)-clust(i));
    end
end