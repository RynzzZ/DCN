function [responses, baselineFR] = processSessionAuditoryData(session, freq, loudness, ephysChannel, responses, varargin)

% settings
s.toneDuration = 0.2; % unit: sec
s.RLFBBNDuration = 0.2; 
s.RLFBBNIntervalDuration = 0.2;
s.RLFBFDuration = 0.2;
s.RLFBFIntervalDuration = 0.2;
s.intervalDuration = 0.4;
s.threshold = -150;
s.thresholdType = 'thresholdDown';
s.plotType = '';
s.suppressFig = true;
s.plotType = []; % choose among 'responseMap', 'RLFBBN' or 'RLFBF'.

if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs

% initializations
rootFolder = 'Z:\Qianyun\DCN\';
gitFolder = 'D:\DCN_Project\Github\DCN\';
sessionFolder = fullfile(rootFolder, 'Data', session);

%-------------------------------load data---------------------------------%
% load spikeAnalyzed.mat
if ~exist(fullfile(sessionFolder, 'spikeAnalyzed.mat'), 'file')
    error('spikeAnalyzed.mat NOT exist! Run analyzeSession.m first!');
else
    disp('loading spikeAnalyzed.mat...')
    load(fullfile(sessionFolder, 'spikeAnalyzed.mat'), 'spike');
    
    % get sampling frequency for audio sig channel and mic channel in spike
    s.audioFs = round(length(spike.audioSignalTimes)/(spike.audioSignalTimes(end) - spike.audioSignalTimes(1)));
    s.micFs = round(length(spike.micSignalTimes)/(spike.micSignalTimes(end) - spike.micSignalTimes(1)));
end

% load neural data
if ~exist(fullfile(sessionFolder, 'sessionEphysInfo.mat'), 'file')
    error('sessioinEphysInfo.mat NOT exist!');
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
    ephysData = getVoltage(data.Data.Data(channelNum_OpenEphys(ephysChannel), :));
    ephysTime = sessionEphysInfo.convertedEphysTimestamps;
end

% get baseline firing rate
startTime = spike.keyboardTimes(find(spike.keyboardTimes>= ephysTime(1), 1, 'first'));
endTime = ephysTime(end)-1;

startInd = find(ephysTime <= startTime, 1, 'last');
endInd = find(ephysTime >= endTime, 1, 'first');

startEphysData = getVoltage(data.Data.Data(channelNum_OpenEphys(ephysChannel), 1:startInd));
endEphysData = getVoltage(data.Data.Data(channelNum_OpenEphys(ephysChannel), endInd:length(ephysTime)));

[~, ~, n1] = crossdet(startEphysData, s.threshold, s.thresholdType);
[~, ~, n2] = crossdet(endEphysData, s.threshold, s.thresholdType);

baselineFR = (n1+n2)/(startTime - ephysTime(1) + 1);

if isempty(responses)
    clear responses
    responses.J = [];
    responses.K = [];
    responses.L = [];
    responses.H = [];
    responses.U = [];
    responses.O = [];
    responses.P = [];
    responses.Y = [];
    responses.E = [];
    responses.B = [];
end

%------------------------plot to determine the threshold------------------%
if ~s.suppressFig
    ind = strfind(convertCharsToStrings(spike.keyboardInput), 'K');
    ind = ind(1);
    startTime = spike.keyboardTimes(ind);
    endTime = startTime + 10;
    
    startInd = find(sessionEphysInfo.convertedEphysTimestamps >= startTime, 1, 'first');
    endInd = find(sessionEphysInfo.convertedEphysTimestamps <= endTime, 1, 'last');
    
    x = sessionEphysInfo.convertedEphysTimestamps(startInd:endInd);
    y = getVoltage(data.Data.Data(channelNum_OpenEphys(ephysChannel), startInd:endInd));
    figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
    plot(x, y);
    axis tight
end
%---------------------------processing the data---------------------------%
disp(session)
disp('start processing...')
resultMatrix = cell(length(spike.keyboardInput), 1);
parfor i = 1:length(spike.keyboardInput)
    if strcmp(s.plotType, 'responseMap')
        if strcmp(spike.keyboardInput(i), 'E') || strcmp(spike.keyboardInput(i), 'B')
            continue
        end
    elseif strcmp(s.plotType, 'RLFBBN')
        if ~strcmp(spike.keyboardInput(i), 'B')
            continue
        end
    elseif strcmp(s.plotType, 'RLFBF')
        if ~strcmp(spike.keyboardInput(i), 'E')
            continue
        end
    end
    
    keyboardTime = spike.keyboardTimes(i);
    if any(strcmp(fieldnames(spike), 'optoTrainTimes'))
        if sum(spike.optoTrainTimes >= keyboardTime & spike.optoTrainTimes <= keyboardTime + 5) > 50
            continue
        end
    end
    
    switch(spike.keyboardInput(i))
        case 'J'
            keyboardTime = spike.keyboardTimes(i);
            toneNumber = length(freq.J);
            [resultMatrix{i, 1}] = processAuditoryData(keyboardTime, spike.audioSignal, spike.audioSignalTimes, ...
                ephysData, ephysTime, toneNumber, s, 'J');
      
            
        case 'K'
            keyboardTime = spike.keyboardTimes(i);
            toneNumber = length(freq.K);
            [resultMatrix{i, 1}] = processAuditoryData(keyboardTime, spike.audioSignal, spike.audioSignalTimes, ...
                ephysData, ephysTime, toneNumber, s, 'K');
            
        case 'L'
            keyboardTime = spike.keyboardTimes(i);
            toneNumber = length(freq.L);
            [resultMatrix{i, 1}] = processAuditoryData(keyboardTime, spike.audioSignal, spike.audioSignalTimes, ...
                ephysData, ephysTime, toneNumber, s, 'L');
            
        case 'H'
            keyboardTime = spike.keyboardTimes(i);
            toneNumber = length(freq.H);
            [resultMatrix{i, 1}] = processAuditoryData(keyboardTime, spike.audioSignal, spike.audioSignalTimes, ...
                ephysData, ephysTime, toneNumber, s, 'H');
            
        case 'U'
            keyboardTime = spike.keyboardTimes(i);
            toneNumber = length(freq.U);
            [resultMatrix{i, 1}] = processAuditoryData(keyboardTime, spike.audioSignal, spike.audioSignalTimes, ...
                ephysData, ephysTime, toneNumber, s, 'U');
            
        case 'O'
            keyboardTime = spike.keyboardTimes(i);
            toneNumber = length(freq.O);
            [resultMatrix{i, 1}] = processAuditoryData(keyboardTime, spike.audioSignal, spike.audioSignalTimes, ...
                ephysData, ephysTime, toneNumber, s, 'O');
            
        case 'P'
            keyboardTime = spike.keyboardTimes(i);
            toneNumber = length(freq.P);
            [resultMatrix{i, 1}] = processAuditoryData(keyboardTime, spike.audioSignal, spike.audioSignalTimes, ...
                ephysData, ephysTime, toneNumber, s, 'P');
            
        case 'Y'
            keyboardTime = spike.keyboardTimes(i);
            toneNumber = length(freq.Y);
            [resultMatrix{i, 1}] = processAuditoryData(keyboardTime, spike.audioSignal, spike.audioSignalTimes, ...
                ephysData, ephysTime, toneNumber, s, 'Y');
            
        case 'B'
            keyboardTime = spike.keyboardTimes(i);
            toneNumber = length(loudness.RLFBBN);
            [resultMatrix{i, 1}] = processAuditoryData(keyboardTime, spike.audioSignal, spike.audioSignalTimes, ...
                ephysData, ephysTime, toneNumber, s, 'B');
            
        case 'E'
            keyboardTime = spike.keyboardTimes(i);
            toneNumber = length(loudness.RLFBF);
            [resultMatrix{i, 1}] = processAuditoryData(keyboardTime, spike.audioSignal, spike.audioSignalTimes, ...
                ephysData, ephysTime, toneNumber, s, 'E');
            
    end
    
    
end

disp('Parfor loop completed! Start formatting the results!')

for i = 1:length(spike.keyboardInput)
    switch(spike.keyboardInput(i))
        case 'J'
            FR = [];
            FR = resultMatrix{i, 1};
            if isempty(responses.J)
                responses.J = FR;
            else
                responses.J = [responses.J; FR];
            end
            
        case 'K'
            FR = [];
            FR = resultMatrix{i, 1};
            if isempty(responses.K)
                responses.K = FR;
            else
                responses.K = [responses.K; FR];
            end
            
        case 'L'
            FR = [];
            FR = resultMatrix{i, 1};
            if isempty(responses.L)
                responses.L = FR;
            else
                responses.L = [responses.L; FR];
            end
            
        case 'H'
            FR = [];
            FR = resultMatrix{i, 1};
            if isempty(responses.H)
                responses.H = FR;
            else
                responses.H = [responses.H; FR];
            end
            
        case 'U'
            FR = [];
            FR = resultMatrix{i, 1};
            if isempty(responses.U)
                responses.U = FR;
            else
                responses.U = [responses.U; FR];
            end
            
        case 'O'
            FR = [];
            FR = resultMatrix{i, 1};
            if isempty(responses.O)
                responses.O = FR;
            else
                responses.O = [responses.O; FR];
            end
            
        case 'P'
            FR = [];
            FR = resultMatrix{i, 1};
            if isempty(responses.P)
                responses.P = FR;
            else
                responses.P = [responses.P; FR];
            end
            
        case 'Y'
            FR = [];
            FR = resultMatrix{i, 1};
            if isempty(responses.Y)
                responses.Y = FR;
            else
                responses.Y = [responses.Y; FR];
            end
            
        case 'E'
            FR = [];
            FR = resultMatrix{i, 1};
            if isempty(responses.E)
                responses.E = FR;
            else
                responses.E = [responses.E; FR];
            end
            
        case 'B'
            FR = [];
            FR = resultMatrix{i, 1};
            if isempty(responses.B)
                responses.B = FR;
            else
                responses.B = [responses.B; FR];
            end         
    end
end



end