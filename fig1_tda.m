%WGGGG
%
%
% aim - to analyze mitochondrial functioning in mice data using TDA
% date - 12.16.2020
% author - saggar@stanford.edu
%

%% load data
clear; clc; close all;
data_cort = readtable('../data/../data/2.5.2021/Mito_delta_scores_updated 2.5.21.xlsx','Sheet','cort');
data_sd = readtable('../data/../data/2.5.2021/Mito_delta_scores_updated 2.5.21.xlsx','Sheet','sd');

% handle missing data
tmp = table2array(data_cort(:,2:6));
tmp(tmp==0) = NaN;
data_cort_fm = fillmissing(tmp, 'linear', 2); % dim=2, interpolating across subjects

tmp = table2array(data_sd(:,2:7));
tmp(tmp==0) = NaN;
data_sd_fm = fillmissing(tmp, 'linear', 2); % dim=2, interpolating across subjects


%% run Mapper for examining differences/similarity across brainregions within each cohort of mice

metricType = 'euclidean';
runType = 'cort';
metaInfo.mf = data_cort.MC_func;
metaInfo.region = data_cort.BrainReg;
metaInfo.anatomy = data_cort.Anatomical;
metaInfo.exc = data_cort.ExcitationVInhibition;

num_k = 5;
res_val = 8;
gain_val = 70;
[nodeTpMat_cort, nodeBynode_cort, tpMat_cort, matVar_cort, filter_cort] = runBDLMapper_wrapper(data_cort_fm, metaInfo, runType, metricType, num_k, res_val, gain_val);

%% visualize tpMat by regions
[data_cort_regsort, indx] =  sortrows(data_cort, 'BrainReg');
imagesc(tpMat_cort(indx,indx)./norm(tpMat_cort(indx, indx), Inf))
set(gcf,'Color','White');
xticks(1:102)
xtickangle(90)
yticks(1:102)
yticklabels(data_cort_regsort.Var1)
xticklabels(data_cort_regsort.Var1)
title('CORT')
for n=1:1:17
    line([0.5,102.5],[6 6]*n+0.5,'LineWidth',2,'Color','w')
    line([6 6]*n+0.5,[0.5,102.5],'LineWidth',2,'Color','w')
end
%% for SD
runType = 'sd';
metaInfo.mf = (data_sd.MC_func);
metaInfo.region = (data_sd.BrainReg);
metaInfo.anatomy = (data_sd.Anatomical);
metaInfo.exc = (data_sd.ExcitationVInhibition);

num_k = 5;
res_val = 8;
gain_val = 70;

[nodeTpMat_sd, nodeBynode_sd, tpMat_sd, matVar_sd, filter_sd] = runBDLMapper_wrapper(data_sd_fm, metaInfo, runType, metricType, num_k ,res_val, gain_val);

%% visualize tpMat by regions
[data_sd_regsort, indx] =  sortrows(data_sd, 'BrainReg');
imagesc(tpMat_sd(indx,indx)./norm(tpMat_sd(indx, indx), Inf),[0 0.055])
set(gcf,'Color','White');
xticks(1:102)
xtickangle(90)
yticks(1:102)
yticklabels(data_sd_regsort.Var1)
xticklabels(data_sd_regsort.Var1)
title('SD')
for n=1:1:17
    line([0.5,102.5],[6 6]*n+0.5,'LineWidth',2,'Color','w')
    line([6 6]*n+0.5,[0.5,102.5],'LineWidth',2,'Color','w')
end

%% participation coefficient values
p_cort = participation_coef(nodeBynode_cort, matVar_cort.nodeComm_region);
z_cort = module_degree_zscore(nodeBynode_cort, matVar_cort.nodeComm_region);
figure; kde2d(p_cort, z_cort); title('CORT');
set(gcf, 'Color', 'w');

p_sd = participation_coef(nodeBynode_sd, matVar_sd.nodeComm_region);
z_sd = module_degree_zscore(nodeBynode_sd, matVar_sd.nodeComm_region);
figure; kde2d(p_sd, z_sd); title('SD');
set(gcf, 'Color', 'w');


%% running Mapper for examining differences/similarity across mice
metricType = 'euclidean';
runType = 'data_all_mat';
metaInfo.mf = data_all.type3;
metaInfo.region = data_all.type2;
metaInfo.anatomy = ones(1,27);
metaInfo.exc = ones(1,27);

num_k = 6;
res_val = 5;
gain_val = 70;

[nodeTpMat_all, nodeBynode_all, tpMat_all, matVar_all, filter_all] = runBDLMapper_wrapper(data_all_mat, metaInfo, runType, metricType, num_k, res_val, gain_val);





%% utility functions
function [nodeTpMat, nodeBynode, tpMat, matVar, filter] = runBDLMapper_wrapper(parcelData, metaInfo, runType, metricType, num_k, res_val, gain_val)

    % run mapper
    data_z = parcelData;
    nclust = 10;
    num_bin_clusters = nclust;

    json_thrval = 0.5;
    json_saveOrNot = 1;

    [nodeTpMat, nodeBynode, tpMat, filter] = runBDLMapper(data_z, metricType, res_val, gain_val, num_k, num_bin_clusters);
    
    % json generated that can be used by DyNeuSR to visualize in Python
    % toolbox
    json_output = sprintf('../d3-nodepies/data/graph_run_%s_metric_%s_res%d_gain%d_numk%d_3bin_nclust%d.json', runType, metricType, res_val, gain_val, num_k, nclust);
    matVar = createJSON_NodePie_cme_enhanced(nodeTpMat, nodeBynode, metaInfo, json_output, json_thrval, json_saveOrNot);    
    matVar.metaInfo = metaInfo;
    
    
end
function [nodeTpMat, nodeBynode, tpMat, filter] = runBDLMapper(data, metricType, res_val, gain_val, num_k, num_bin_clusters)

    X = data;


    resolution = [res_val res_val];
    gain = gain_val;
    
    fprintf(1,'Estimating distance matrix\n');
    tic
    distMat = estimateDistance(X, metricType);
    toc
    
    
    fprintf(1,'Estimating knn graph\n');
    tic
    % create knn graph, estimate geodesic distances, embed using cmdscale and apply mapper
    [knnGraphTbl, knnGraph_dense_bin, knnGraph_dense_wtd, knnGraph_dense_bin_conn, knnGraph_dense_wtd_conn]= createPKNNG_bdl(distMat, num_k);

    knn_g_wtd = graph(knnGraph_dense_bin_conn);

    % estimate geodesic distances
    dist_geo_wtd = round(distances(knn_g_wtd,'Method','positive'));
    toc
    
    fprintf(1,'Estimating embedding\n');
    tic
    
    % embed using cmdscale
    [y,e] = cmdscale(dist_geo_wtd);

    filter = [y(:,1), y(:,2)];
    toc
    
    fprintf(1,'Running mapper\n');
    tic
    
    
    [adja, adja_pruned, pts_in_vertex, pts_in_vertex_pruned] = mapper2d_bdl_hex_binning(distMat, filter, resolution, gain, num_bin_clusters, 3); %using triangulation
    
    toc
    
    fprintf(1,'Creating final output\n');
    tic
    % creating matrices for d3 visualization
    numNodes = length(pts_in_vertex_pruned);
    numTp = size(X,1);
    nodeTpMat = zeros(numNodes, numTp);
    for node = 1:1:numNodes
        tmp = pts_in_vertex_pruned{node};
        nodeTpMat(node, tmp) = 1;
    end

    nodeBynode = adja_pruned;
    tpMat = getMatTp_wtd(nodeBynode, nodeTpMat);
    fprintf(1,'Done\n');
    toc

end

function distMat = estimateDistance(X, metricType)
    distMat = squareform(pdist(X, metricType));
end
function matVar = createJSON_NodePie_cme_enhanced(nodeTpMat, nodeBynode, metaInfo, outpath, thrval, saveOrNot)

    regionInfo = metaInfo.region;
    mfInfo = metaInfo.mf;
    anaInfo = metaInfo.anatomy;
    excInfo = metaInfo.exc;
    
    numNodes = size(nodeTpMat,1);
    numTRs = size(nodeTpMat,2);
    thr = thrval; % 1 std. dev. above the mean for rsn metaInfor

    d = [];
    d.directed = 0;
    d.multigraph = 0;
    d.graph = {};
    tmp = sum(nodeTpMat,2);
    numLinks = 1;
    
    prop_nodes_region = [];
    prop_nodes_mf = [];
    prop_nodes_anatomy = [];
    prop_nodes_exc = [];

    numMetaGroups_region = length(unique(regionInfo));
    numMetaGroups_mf = length(unique(mfInfo));
    numMetaGroups_anatomy = length(unique(anaInfo));
    numMetaGroups_exc = length(unique(excInfo));
    
    nodeComm_region = [];
    nodeComm_mf = [];
    nodeComm_anatomy = [];
    nodeComm_exc = [];
    
    
    for node = 1:1:numNodes
        
        d.nodes{node}.row_count = tmp(node);
        trs = find(nodeTpMat(node,:));
        d.nodes{node}.trs = sprintf('[%s]',num2str(find(nodeTpMat(node,:))));

        prop_region = zeros(numMetaGroups_region,1); 
        prop_mf = zeros(numMetaGroups_mf,1); 
        prop_anatomy = zeros(numMetaGroups_anatomy,1); 
        prop_exc = zeros(numMetaGroups_exc,1); 
        
        for tr = trs
            prop_region(regionInfo(tr)) = prop_region(regionInfo(tr)) + 1;
            prop_mf(mfInfo(tr)) = prop_mf(mfInfo(tr)) + 1;            
            prop_anatomy(anaInfo(tr)) = prop_anatomy(anaInfo(tr)) + 1;            
            prop_exc(excInfo(tr)) = prop_exc(excInfo(tr)) + 1;            
        end

        [~,nodeComm_region(node)] = max(prop_region);
        prop_nodes_region(node,:) = prop_region;        
        
        [~,nodeComm_mf(node)] = max(prop_mf);
        prop_nodes_mf(node,:) = prop_mf;        
        
        [~,nodeComm_anatomy(node)] = max(prop_anatomy);
        prop_nodes_anatomy(node,:) = prop_anatomy;        
        
        [~,nodeComm_exc(node)] = max(prop_exc);
        prop_nodes_exc(node,:) = prop_exc;    
        
        for grp = 1:1:numMetaGroups_region
            d.nodes{node}.proportions_region{grp}.group = grp;
            d.nodes{node}.proportions_region{grp}.value = prop_region(grp);
            d.nodes{node}.proportions_region{grp}.row_count = tmp(node);  
        end
        
        for grp = 1:1:numMetaGroups_region
            
            if grp <= numMetaGroups_mf
                d.nodes{node}.proportions_mf{grp}.group = grp;
                d.nodes{node}.proportions_mf{grp}.value = prop_mf(grp);
                d.nodes{node}.proportions_mf{grp}.row_count = tmp(node);  
            else
                d.nodes{node}.proportions_mf{grp}.group = grp;
                d.nodes{node}.proportions_mf{grp}.value = 0;
                d.nodes{node}.proportions_mf{grp}.row_count = tmp(node);                 
            end
  
        end
        
        for grp = 1:1:numMetaGroups_region
            
            if grp <= numMetaGroups_anatomy
                d.nodes{node}.proportions_anatomy{grp}.group = grp;
                d.nodes{node}.proportions_anatomy{grp}.value = prop_anatomy(grp);
                d.nodes{node}.proportions_anatomy{grp}.row_count = tmp(node);  
            else
                d.nodes{node}.proportions_anatomy{grp}.group = grp;
                d.nodes{node}.proportions_anatomy{grp}.value = 0;
                d.nodes{node}.proportions_anatomy{grp}.row_count = tmp(node);                 
            end
  
        end        

        for grp = 1:1:numMetaGroups_region
            
            if grp <= numMetaGroups_exc
                d.nodes{node}.proportions_exc{grp}.group = grp;
                d.nodes{node}.proportions_exc{grp}.value = prop_exc(grp);
                d.nodes{node}.proportions_exc{grp}.row_count = tmp(node);  
            else
                d.nodes{node}.proportions_exc{grp}.group = grp;
                d.nodes{node}.proportions_exc{grp}.value = 0;
                d.nodes{node}.proportions_exc{grp}.row_count = tmp(node);                 
            end
  
        end          
        
       
        d.nodes{node}.id = node-1;

        links = find(nodeBynode(node,:));
        if ~isempty(links)        
            for l = 1:1:length(links)
                d.links{numLinks}.source = node - 1;
                d.links{numLinks}.target = links(l) - 1;
                numLinks = numLinks + 1;
            end

        end

    end
    
    if saveOrNot == 1 % only save when asked!
        savejson('',d,outpath);
    end
    
    matVar.nodeComm_region = nodeComm_region;
    matVar.nodeComm_mf = nodeComm_mf;
    matVar.nodeComm_anatomy = nodeComm_anatomy;
    matVar.nodeComm_exc = nodeComm_exc;
    
    
    matVar.prop_nodes_region = prop_nodes_region;
    matVar.prop_nodes_mf = prop_nodes_mf;
    matVar.prop_nodes_anatomy = prop_nodes_anatomy;
    matVar.prop_nodes_exc = prop_nodes_exc;

end