function [spike, videoTracking, sessionEphysInfo, getVoltage, data, channelNum_OpenEphys] = loadSessionData(session, varargin)

% initializations
rootFolder = 'Z:\Qianyun\DCN\';
gitFolder = 'D:\DCN_Project\Github\DCN\';
sessionFolder = fullfile(rootFolder, 'Data', session);

% load spikeAnalyzed.mat
if ~exist(fullfile(sessionFolder, 'spikeAnalyzed.mat'), 'file')
    warning('spikeAnalyzed.mat NOT exist! Run analyzeSession.m first!');
    spike = [];
else
    disp('loading spikeAnalyzed.mat...')
    load(fullfile(sessionFolder, 'spikeAnalyzed.mat'), 'spike');
    
    % get sampling frequency for audio sig channel and mic channel in spike
    s.audioFs = round(length(spike.audioSignalTimes)/(spike.audioSignalTimes(end) - spike.audioSignalTimes(1)));
    s.micFs = round(length(spike.micSignalTimes)/(spike.micSignalTimes(end) - spike.micSignalTimes(1)));
end


% load cameraAnalyzed.mat
if ~exist(fullfile(sessionFolder, 'cameraAnalyzed.mat'), 'file')
    warning('cameraAnalyzed.mat NOT exist! Run analyzeSession.m first!');
    videoTracking = [];
else
    disp('loading cameraAnalyzed.mat...')
    load(fullfile(sessionFolder, 'cameraAnalyzed.mat'));
    if ~any(strcmp('jawDistance', videoTracking.Properties.VariableNames))
        error('JawDistance does NOT exist in videoTracking! Run analyzeSession.m first!');
    end
end

% load neural data
if ~exist(fullfile(sessionFolder, 'sessionEphysInfo.mat'), 'file')
    warning('sessioinEphysInfo.mat NOT exist!');
    sessionEphysInfo = [];
    getVoltage = [];
    data = [];
else
    % load sessionEphysInfo
    disp('loading ephys data...')
    load(fullfile(rootFolder, 'Data', session, 'sessionEphysInfo.mat'), 'sessionEphysInfo');
    mapFile = sessionEphysInfo.mapFile;
    load(fullfile('Z:\obstacleData\ephys\channelMaps\kilosort', [mapFile, '.mat']), 'channelNum_OpenEphys');
    
    % function to extract voltage from binary file
    getVoltage = @(data) ...
        double(data)*sessionEphysInfo.bitVolts; % extract voltage from memmapfile, converting to votlage, highpassing, and only return specific channel
    
    % load raw ephys voltage
    contFiles = dir(fullfile(sessionEphysInfo.ephysFolder, '*.continuous'));
    data = memmapfile(fullfile(sessionEphysInfo.ephysFolder, [contFiles(1).name(1:end-12), 's.dat']), ...
        'Format', {'int16', [sessionEphysInfo.channelNum sessionEphysInfo.smps], 'Data'}, 'Writable', false);
end

end