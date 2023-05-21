%WGGGG
%
%
% aim - to analyze mitochondrial functioning in mice data using multi-slice
% community detection algorithm
% date - 12.16.2020
% author - saggar@stanford.edu
%


%% multislice
clear; clc;

data_all = readtable('../data/2.5.2021/Mito-behavior_correlations_Ayelet_updated.xlsx','Sheet','Sheet1');
metainfo_tbl = readtable('../data/2.5.2021/data_all_metainfo.xlsx');

data_controls = data_all(1:11,:);
data_stressed_cort = data_all(12:16,:);
data_stressed_sd = data_all(17:22,:);
data_stressed_recovered = data_all(23:27,:);

%within group missing interpolation
data_controls_mat = table2array(data_controls(:,6:102+6-1));
data_controls_mat = fillmissing(data_controls_mat, 'linear', 1); % dim=1, interpolating across subjects

data_stressed_cort_mat = table2array(data_stressed_cort(:,6:102+6-1));
data_stressed_cort_mat = fillmissing(data_stressed_cort_mat, 'linear', 1); % dim=1, interpolating across subjects

data_stressed_sd_mat = table2array(data_stressed_sd(:,6:102+6-1));
data_stressed_sd_mat = fillmissing(data_stressed_sd_mat, 'linear', 1); % dim=1, interpolating across subjects

data_stressed_recovered_mat = table2array(data_stressed_recovered(:,6:102+6-1));
data_stressed_recovered_mat = fillmissing(data_stressed_recovered_mat, 'linear', 1); % dim=1, interpolating across subjects

data_all_mat = [data_controls_mat;data_stressed_cort_mat;data_stressed_sd_mat;data_stressed_recovered_mat];

% create slices for each mito feature
A = {};
for sl = 1:1:6
    A{sl} = suppress_diag(corr(data_all_mat(:,find(metainfo_tbl.mc_func==sl))));       
end

%%
% visualize A
reg_names = {'cereb','amyg','ca3','hypo','cpu','ofc','dgd','sn','mpfc','nac','mot','thal','vta','dgv','pag','vn','vis'};
mf_names = {'CS','CI','CII','CIV','mtDNA','MHI'};
ax=[];
for sl = 1:1:6
   figure;  set(gcf,'color','w'); %imagesc(A{sl});
   ax(sl)= gca;
   s= surfl(A{sl}); axis([0 20 0 20 -100 100])
    rotate(s,[10 0 0], 90)
     axis([0 20 -200 200 -20 20])
      box off
    grid off
end
close all;

[Q, S, B] = categoricalMultiSlice_mf(A, 1, 0.1);
% for each slice compute the allegiance matrix
alleg_mat = zeros(6,17,17);
for sl = 1:1:6
    tmp = S(:,sl);
    tmp1 = unique(tmp);
    tmp2 = zeros(17,17);
    for t = 1:1:length(tmp1)
       tmp2(find(tmp==tmp1(t)), find(tmp==tmp1(t))) = 1;
    end
    alleg_mat(sl,:,:) = tmp2;
end

[ciu,mean_q] = community_louvain_iterative_local(squeeze(mean(alleg_mat)))
[comm, comm2, order] = arrangeComm(ciu,squeeze(mean(alleg_mat)))
figure; imagesc(comm)
set(gcf,'color','w'); 
xticks(1:17); xticklabels(reg_names(order));
yticks(1:17); yticklabels(reg_names(order));
comm_uthr = comm;
comm(comm<0.5) = 0; % keeping only top 30% 
g= graph(comm, 'omitselfloops');
figure; plot(g,'Layout','force','NodeLabel',reg_names(order),'MarkerSize',20,'LineWidth',g.Edges.Weight*2,'NodeCData',ciu(order),'NodeFontSize',20)


%% utility functions
function [ciu,mean_q] = community_louvain_iterative_local(A, obj_fun, gamma)

    if nargin < 2
       obj_fun = 'modularity';
       gamma = 1;
    end
    m=[];
    q=[];
    for iter=1:1:1000
        %[m(:,iter),q(iter)] = community_louvain(A,[],[],'negative_asym');
        [m(:,iter),q(iter)] = community_louvain(A,gamma,[],obj_fun);
    end
    d=agreement_weighted(m,q);
    %tau = prctile(d(:),66);
    ciu=consensus_und(d,.001,100);
    mean_q = mean(q);

end
% from http://netwiki.amath.unc.edu/GenLouvain/GenLouvain
function [Q,S, B] = categoricalMultiSlice_mf(A, gamma, omega)
    N=length(A{1});
    T=length(A);
    B=spalloc(N*T,N*T,(N+T)*N*T);
    twomu=0;
    for s=1:T
        k=sum(A{s});
        twom=sum(k);
        twomu=twomu+twom;
        indx=[1:N]+(s-1)*N;
        B(indx,indx)=A{s}-gamma*k'*k/twom;
    end
    twomu=twomu+T*omega*N*(T-1);
    all2all = N*[(-T+1):-1,1:(T-1)];
    B = B + omega*spdiags(ones(N*T,2*T-2),all2all,N*T,N*T);
    [S,Q] = genlouvain(B,[],[],0);
    %[S,Q] = genlouvain(B,[],[], [], 'moverandw');
    Q = Q/twomu
    S = reshape(S,N,T);

end
