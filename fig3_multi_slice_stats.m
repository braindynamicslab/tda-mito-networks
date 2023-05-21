%WGGGG
%
%
% aim - to analyze mitochondrial functioning in mice data using multi-slice
% community detection algorithm. Examining the degree to which the observed
% communities in the Mito-features are also preserved in Gene co-expression
% data and structural connectome data.
%
% date - 12.16.2020
% author - saggar@stanford.edu
%


%% all conn
clear; clc;
gene_all = readtable('../data/connectome/all_conn.xlsx','Sheet','gene_str');
str_all = readtable('../data/connectome/all_conn.xlsx','Sheet','str_conn');
mf_all = readtable('../data/connectome/all_conn.xlsx','Sheet','mf_conn_str');
gene_all_arr = table2array(gene_all(:,2:end-1));
[~,tmp] = sort(gene_all.order_mf);
gene_all_arr = gene_all_arr(tmp,tmp);

str_all_arr = table2array(str_all(:,2:end-1));
[~,tmp] = sort(str_all.order_mf);
str_all_arr = str_all_arr(tmp,tmp);

mf_all_arr = table2array(mf_all(:,2:end-1));

figure; imagesc(mf_all_arr); title('MF connectivity');
set(gcf,'color','w'); colorbar;
xticks(1:15); xticklabels(mf_all.Var1);
xtickangle(90)
yticks(1:15); yticklabels(mf_all.Var1);

figure; imagesc(log10(str_all_arr)); title('Struc connectivity');
set(gcf,'color','w'); colorbar;
xticks(1:15); xticklabels(mf_all.Var1);
xtickangle(90)
yticks(1:15); yticklabels(mf_all.Var1);

figure; imagesc(gene_all_arr); title('Gene connectivity');
set(gcf,'color','w'); colorbar;
xticks(1:15); xticklabels(mf_all.Var1);
xtickangle(90)
yticks(1:15); yticklabels(mf_all.Var1);

%% Two difference comparisons used: Richardi et al style and Modularity based
clc;
ciu = mf_all.ciu;

str_all_arr = log10(str_all_arr)+10; 
str_within_gene = strength_within(gene_all_arr, ciu, max(ciu));
str_within_mf = strength_within(mf_all_arr, ciu, max(ciu));
str_within_str = strength_within(str_all_arr, ciu, max(ciu));

mod_gene = calMod(gene_all_arr,ciu);
mod_mf = calMod(mf_all_arr,ciu);
mod_str = calMod(str_all_arr,ciu);

str_within_gene_rnd = [];
str_within_mf_rnd = [];
str_within_str_rnd = [];
mod_gene_rnd = [];
mod_mf_rnd = [];
mod_str_rnd = [];

nperm = 10000;
for perm = 1:1:nperm
    nciu = ciu(randperm(length(ciu)));
    str_within_gene_rnd(perm) = strength_within(gene_all_arr, nciu, max(nciu));
    str_within_mf_rnd(perm) = strength_within(mf_all_arr, nciu, max(nciu));
    str_within_str_rnd(perm) = strength_within(str_all_arr, nciu, max(nciu));
    
    mod_gene_rnd(perm) = calMod(gene_all_arr,nciu);
    mod_mf_rnd(perm) = calMod(mf_all_arr,nciu);
    mod_str_rnd(perm) = calMod(str_all_arr,nciu);

end

(sum(str_within_gene_rnd>str_within_gene)+1)/(nperm+1)
(sum(str_within_mf_rnd>str_within_mf)+1)/(nperm+1)
(sum(str_within_str_rnd>str_within_str)+1)/(nperm+1)
(sum(mod_gene_rnd>mod_gene)+1)/(nperm+1)
(sum(mod_mf_rnd>mod_mf)+1)/(nperm+1)
(sum(mod_str_rnd>mod_str)+1)/(nperm+1)

nhist(str_within_mf_rnd)
set(gcf,'color','w');
hold on
plot(str_within_mf,1,'ro','MarkerSize',10,'MarkerFaceColor','r')
axis([0 4 0 700])
clf
nhist(str_within_str_rnd)
hold on
plot(str_within_str,1,'ro','MarkerSize',10,'MarkerFaceColor','r')
clf
nhist(str_within_gene_rnd)
hold on
plot(str_within_gene,1,'ro','MarkerSize',10,'MarkerFaceColor','r')
clf
nhist(mod_mf_rnd)
hold on
plot(mod_mf,1,'ro','MarkerSize',10,'MarkerFaceColor','r')
mod_mf
axis([0 .3 0 700])
clf
nhist(mod_gene_rnd)
hold on
plot(mod_gene,1,'ro','MarkerSize',10,'MarkerFaceColor','r')
clf
nhist(mod_str_rnd)


%% utility functions
function str_within = strength_within(w, ciu, num_comm)
    str_within = [];
    for c = 1:1:num_comm
        str_within = [str_within; mean(w(find(ciu==c),find(ciu==c)),2)./mean(w(find(ciu==c),find(ciu~=c)),2)];
    end
    str_within = mean(str_within);
end
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

