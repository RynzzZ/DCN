function [spkInds, unit_ids] = getGoodSpkInds(session)

% given session name, returns cell array where each entry is vector of
% spike indices (with respect to open ephys clock) for units marked as 
% 'good' in kilosort


% initializations
rootFolder = 'Z:\Qianyun\DCN\';
gitFolder = 'D:\DCN_Project\Github\DCN\';
sessionFolder = fullfile(rootFolder, 'Data', session);
GITDIR = 'D:\Cerebellum_DCN_Project\Github';

ephysFolder = dir(fullfile(sessionFolder, [session(1:4), '-*']));
ephysFolder = fullfile(ephysFolder.folder, ephysFolder.name, 'Record Node 101');
addpath(genpath(fullfile(GITDIR, 'npy-matlab')));
addpath(genpath(fullfile(GITDIR, 'analysis-tools')));
allSpkInds = readNPY(fullfile(ephysFolder, 'spike_times.npy'));
clusters = readNPY(fullfile(ephysFolder, 'spike_clusters.npy'));


% old kilosort format
if exist(fullfile(ephysFolder, 'cluster_groups.csv'), 'file')
    clusterGroups = readtable(fullfile(ephysFolder, 'cluster_groups.csv'));
    unit_ids = clusterGroups.cluster_id(strcmp(clusterGroups.group, 'good'));

% new kilosort format
elseif exist(fullfile(ephysFolder, 'cluster_group.tsv'), 'file')  % new kilosort1 format
    clusterInfo = tdfread(fullfile(ephysFolder, 'cluster_info.tsv'));
    bins = all(clusterInfo.group(:,1:4)=='good',2);  % this is a bit of a hack
    unit_ids = clusterInfo.id(bins);
end


% get spike times for individual units
spkInds = cell(1,length(unit_ids));
for i = 1:length(unit_ids)
    spkInds{i} = allSpkInds(clusters==unit_ids(i));
end





