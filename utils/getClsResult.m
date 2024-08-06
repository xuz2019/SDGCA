function result = getClsResult(CA,clsNum)
    CA(CA < 0) = 0;
    CA(CA > 1) = 1;
    CA = max(CA,CA');
    s = squareform(CA - diag(diag(CA)),'tovector');
    d = 1 - s;
    result = cluster(linkage(d,'average'),'maxclust',clsNum);
end


