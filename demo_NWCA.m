clc;clear;
close all;
addpath(genpath(pwd));

%%
dataName = 'Ecoli';
dataPath = ['./datasets/', dataName, '.mat'];
load(dataPath);

%%
rng(19);
M = 20;
cntTimes = 20;
poolSize = size(members, 2);

%%
bcIdx = zeros(cntTimes, M);
for i = 1:cntTimes
    tmp = randperm(poolSize);
    bcIdx(i,:) = tmp(1:M);
end

%%
if (min(gt) == 0)
    gt = gt + 1;
end
clsNums = length(unique(gt));
[N, ~] = size(members);

%%
NMI = zeros(1, cntTimes);
ARI = NMI;
F = NMI;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NWCApara =  struct('Ecoli', 0.09);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
for runIdx = 1:cntTimes
% parfor runIdx = 1:cntTimes
    baseCls = members(:,bcIdx(runIdx,:));
    [bcs, baseClsSegs] = getAllSegs(baseCls);
    CA = baseClsSegs' * baseClsSegs / M;

    %%
    NWCA = computeNWCA(baseClsSegs, computeNECI(bcs, baseClsSegs, NWCApara.(dataName)), M);
    result = getClsResult(NWCA,clsNums);
    if (min(result) == 0)
        result = result + 1;
    end
    NMI(runIdx) = compute_nmi(result, gt);
    ARI(runIdx) = RandIndex(result, gt);
    F(runIdx) = compute_f(result, gt);

end

%%
disp('           mean    variance')
disp(['NMI       ', num2str(mean(NMI), '%.3f'), '     ', num2str(std(NMI), '%.3f')]);
disp(['ARI       ', num2str(mean(ARI), '%.3f'), '     ', num2str(std(ARI), '%.3f')]);
disp(['F-score   ', num2str(mean(F), '%.3f'), '     ', num2str(std(F), '%.3f')]);