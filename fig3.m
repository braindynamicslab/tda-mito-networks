%WGGGG
%
%
%
%
% aim - generate figures for figure 3
% author - saggar@stanford.edu
% date - 2.19.2023

%% fig 3a
clear; clc; close all;

rawdata = readtable('../manuscript/revision/fig3g_data_matrix1.xlsx','Sheet','raw_data_1');
[r,p]=corr(table2array(rawdata(:,6:end)),'rows','pairwise');
fdr=mafdr(p(:),'BHFDR','True');
r_fdr = r(:);
r_fdr(fdr>0.05) = NaN;

figure, imagesc(reshape(r_fdr, [102 102]))
xticks(1:102)
yticks(1:102)
xticklabels(rawdata.Properties.VariableNames(6:end))
xtickangle(90)
yticklabels(rawdata.Properties.VariableNames(6:end))
colormap("hot");

%% fig 3g

% compute network level values
net1 = {'motor','CPu','OFC','mPFC','Visual','NAc'};
net2 = {'cereb','VN','VTA','thal','CA3','DGv','DGd'};
net3 = {'SN','PAG','Hypo','amyg'};

net_means1_reg = [];

net_means1 = [];
for n = 1:1:length(net1)
    tmp = find(contains(rawdata.Properties.VariableNames,net1{n}));
    tmp2 = table2array(rawdata(:,tmp));
    tmp2 = (tmp2 - nanmean(tmp2))./nanstd(tmp2);
    net_means1 = [net_means1 nanmean(tmp2,2)];
end
net_means1_reg = net_means1;
net_means1 = nanmean(net_means1,2);

net_means2 = [];
net_means2_reg = [];
for n = 1:1:length(net2)
    tmp = find(contains(rawdata.Properties.VariableNames,net2{n}));
    tmp2 = table2array(rawdata(:,tmp));
    tmp2 = (tmp2 - nanmean(tmp2))./nanstd(tmp2);
    net_means2 = [net_means2 nanmean(tmp2,2)];
end
net_means2_reg = net_means2;
net_means2 = nanmean(net_means2,2);

net_means3 = [];
net_means3_reg = [];
for n = 1:1:length(net3)
    tmp = find(contains(rawdata.Properties.VariableNames,net3{n}));
    tmp2 = table2array(rawdata(:,tmp));
    tmp2 = (tmp2 - nanmean(tmp2))./nanstd(tmp2);
    net_means3 = [net_means3 nanmean(tmp2,2)];
end
net_means3_reg = net_means3;
net_means3 = nanmean(net_means3,2);

%% correlations
r=[];
p=[];
[r(1),p(1)]=corr(net_means1,rawdata.zScoreEPM,'rows','pairwise')
[r(2),p(2)]=corr(net_means1,rawdata.OFTZScored,'rows','pairwise')
[r(3),p(3)]=corr(net_means1,rawdata.NSF,'rows','pairwise')
[r(4),p(4)]=corr(net_means1,rawdata.SocialAvoidance,'rows','pairwise')

[r(5),p(5)]=corr(net_means2,rawdata.zScoreEPM,'rows','pairwise')
[r(6),p(6)]=corr(net_means2,rawdata.OFTZScored,'rows','pairwise')
[r(7),p(7)]=corr(net_means2,rawdata.NSF,'rows','pairwise')
[r(8),p(8)]=corr(net_means2,rawdata.SocialAvoidance,'rows','pairwise')

[r(9),p(9)]=corr(net_means3,rawdata.zScoreEPM,'rows','pairwise')
[r(10),p(10)]=corr(net_means3,rawdata.OFTZScored,'rows','pairwise')
[r(11),p(11)]=corr(net_means3,rawdata.NSF,'rows','pairwise')
[r(12),p(12)]=corr(net_means3,rawdata.SocialAvoidance,'rows','pairwise')

%% regional correlations
for n = 1:1:length(net1)
    [r1_regEPM(n), p1_regEPM(n)] = corr(net_means1_reg(:,n), rawdata.zScoreEPM,'rows','pairwise');
    [r1_regOFT(n), p1_regOFT(n)] = corr(net_means1_reg(:,n), rawdata.OFTZScored,'rows','pairwise');
    [r1_regNSF(n), p1_regNSF(n)] = corr(net_means1_reg(:,n), rawdata.NSF,'rows','pairwise');        
    [r1_regSA(n), p1_regSA(n)] = corr(net_means1_reg(:,n), rawdata.SocialAvoidance,'rows','pairwise');    
end

for n = 1:1:length(net2)
    [r2_regEPM(n), p2_regEPM(n)] = corr(net_means2_reg(:,n), rawdata.zScoreEPM,'rows','pairwise');
    [r2_regOFT(n), p2_regOFT(n)] = corr(net_means2_reg(:,n), rawdata.OFTZScored,'rows','pairwise');
    [r2_regNSF(n), p2_regNSF(n)] = corr(net_means2_reg(:,n), rawdata.NSF,'rows','pairwise');        
    [r2_regSA(n), p2_regSA(n)] = corr(net_means2_reg(:,n), rawdata.SocialAvoidance,'rows','pairwise');    
end

for n = 1:1:length(net3)
    [r3_regEPM(n), p3_regEPM(n)] = corr(net_means3_reg(:,n), rawdata.zScoreEPM,'rows','pairwise');
    [r3_regOFT(n), p3_regOFT(n)] = corr(net_means3_reg(:,n), rawdata.OFTZScored,'rows','pairwise');
    [r3_regNSF(n), p3_regNSF(n)] = corr(net_means3_reg(:,n), rawdata.NSF,'rows','pairwise');        
    [r3_regSA(n), p3_regSA(n)] = corr(net_means3_reg(:,n), rawdata.SocialAvoidance,'rows','pairwise');    
end
