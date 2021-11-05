function sessionEphysInfo = getSessionEphysInfo(session)

% gets ephys folder name, base name for recording, number of channels,
% sampling frequency, number of samples, bitvols conversion factor

% get name of ephys folder
rootFolder = 'Z:\Qianyun\DCN\';
gitFolder = 'D:\DCN_Project\Github\';
sessionFolder = fullfile(rootFolder, 'Data', session);
sessionFolderFiles = dir(fullfile(rootFolder, 'Data', session, '20*'));
ephysFolder = [sessionFolderFiles([sessionFolderFiles.isdir] == 1).name];
ephysFolder = fullfile(sessionFolder, ephysFolder, 'Record Node 101\');
sessionEphysInfo.ephysFolder = ephysFolder;

% get source name (e.g. 100, 107) and number of channels
contFiles = dir(fullfile(sessionEphysInfo.ephysFolder, '*.continuous'));
contFiles = contFiles(~contains({contFiles.name}, 'AUX')); % remove AUX channels
sessionEphysInfo.channelNum = length(contFiles);
sessionEphysInfo.fileNameBase = contFiles(1).name(1:3);

% get probe mapping file
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')
ephysInfo = readtable(fullfile(rootFolder, 'spreadSheets', 'ephysInfo.xlsx'), 'Sheet', 'ephysInfo');
warning('on', 'MATLAB:table:ModifiedAndSavedVarnames')
sessionEphysInfo.mapFile = ephysInfo.map{strcmp(session, ephysInfo.session)};

% get fs, microvolts conversion factor, and number of samples
addpath(fullfile(gitFolder, 'analysis-tools'))
[~, sessionEphysInfo.timeStamps, info] = load_open_ephys_data_faster(...
    fullfile(sessionEphysInfo.ephysFolder, contFiles(end).name)); % for some reason taking the first contFiles failed on one session...
sessionEphysInfo.fs = info.header.sampleRate;
sessionEphysInfo.bitVolts = info.header.bitVolts;

% get number of samples
file = fullfile(sessionEphysInfo.ephysFolder, [contFiles(1).name(1:end-12), 's.dat']);
temp = dir(file);
sessionEphysInfo.smps = temp.bytes/2/sessionEphysInfo.channelNum; % 2 bytes per sample

% channel used for synchronization
sessionEphysInfo.syncSignal = ephysInfo.syncSignal{strcmp(session, ephysInfo.session)};