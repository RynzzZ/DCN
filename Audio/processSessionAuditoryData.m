function [responses, optoResponses, intervalResponses, optoIntervalResponses] = processSessionAuditoryData(session, freq, loudness, ephysChannel, responses, optoResponses, intervalResponses, optoIntervalResponses, varargin)

% settings
s.RMToneDuration = 0.2; % unit: sec
s.RMIntervalDuration = 0.4;
s.BBNDuration = 0.2; 
s.BBNIntervalDuration = 0.2;
s.BFDuration = 0.2;
s.BFIntervalDuration = 0.2;

s.optoTrainFreq = 50;
s.threshold = -420;
s.thresholdType = 'thresholdDown';
s.plotType = '';
s.suppressFig = true;
s.plotType = []; % choose among 'responseMap', 'RLFBBN' or 'RLFBF'.
s.toneEphysLag = 0.00;

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
    
    % determine target channel
    if contains(sessionEphysInfo.ephysFolder, '101')
        ephysChannel = channelNum_OpenEphys(ephysChannel);
    end
    
    % load raw ephys voltage
    contFiles = dir(fullfile(sessionEphysInfo.ephysFolder, '*.continuous'));
    data = memmapfile(fullfile(sessionEphysInfo.ephysFolder, [contFiles(1).name(1:end-12), 's.dat']), ...
        'Format', {'int16', [sessionEphysInfo.channelNum sessionEphysInfo.smps], 'Data'}, 'Writable', false);
    ephysData = getVoltage(data.Data.Data(ephysChannel, :));
    ephysTime = sessionEphysInfo.convertedEphysTimestamps;
    
    disp('finish loading ephys data!');
end

resultMatrix = cell(length(spike.keyboardInput), 2);
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
    responses.baselineFR = [];
end

if isempty(intervalResponses)
    clear intervalResponses
    intervalResponses.J = [];
    intervalResponses.K = [];
    intervalResponses.L = [];
    intervalResponses.H = [];
    intervalResponses.U = [];
    intervalResponses.O = [];
    intervalResponses.P = [];
    intervalResponses.Y = [];
    intervalResponses.E = [];
    intervalResponses.B = [];
    intervalResponses.baselineFR = [];
end

if any(strcmp(fieldnames(spike), 'optoTrainTimes'))
    optoResultMatrix = cell(length(spike.keyboardInput), 2);
    if isempty(optoResponses)
        clear optoResponses
        optoResponses.J = [];
        optoResponses.K = [];
        optoResponses.L = [];
        optoResponses.H = [];
        optoResponses.U = [];
        optoResponses.O = [];
        optoResponses.P = [];
        optoResponses.Y = [];
        optoResponses.E = [];
        optoResponses.B = [];
    end
    
    if isempty(optoIntervalResponses)
        clear optoIntervalResponses
        optoIntervalResponses.J = [];
        optoIntervalResponses.K = [];
        optoIntervalResponses.L = [];
        optoIntervalResponses.H = [];
        optoIntervalResponses.U = [];
        optoIntervalResponses.O = [];
        optoIntervalResponses.P = [];
        optoIntervalResponses.Y = [];
        optoIntervalResponses.E = [];
        optoIntervalResponses.B = [];
    end
    
end


%------------------------plot to determine the threshold------------------%
if ~s.suppressFig
    ind = strfind(convertCharsToStrings(spike.keyboardInput), 'E');
    ind = ind(3);
    startTime = spike.keyboardTimes(ind);
    endTime = startTime + 20;
    
    startInd = find(sessionEphysInfo.convertedEphysTimestamps >= startTime, 1, 'first');
    endInd = find(sessionEphysInfo.convertedEphysTimestamps <= endTime, 1, 'last');
    
    x = sessionEphysInfo.convertedEphysTimestamps(startInd:endInd);
    y = getVoltage(data.Data.Data(ephysChannel, startInd:endInd));
    figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
    plot(x, y); hold on;
    plot([x(1), x(end)], [s.threshold, s.threshold], '-r');
    axis tight  
end


% get baseline firing rate
TimeDuration = 5;
ind = strfind(convertCharsToStrings(spike.keyboardInput), 'B');
ind = ind(1);
startTime = spike.keyboardTimes(ind) - TimeDuration;
endTime = ephysTime(end) - TimeDuration;

startInd = find(ephysTime >= startTime, 1, 'first');
startInd2 = find(ephysTime <= startTime + TimeDuration, 1, 'last');
endInd = find(ephysTime >= endTime, 1, 'first');

startEphysData = getVoltage(data.Data.Data(ephysChannel, startInd:startInd2));
endEphysData = getVoltage(data.Data.Data(ephysChannel, endInd:length(ephysTime)));

[~, ~, n1] = crossdet(startEphysData, s.threshold, s.thresholdType);
[~, ~, n2] = crossdet(endEphysData, s.threshold, s.thresholdType);

baselineFR = (n1+n2)/( TimeDuration*2 );
responses.baselineFR = baselineFR;
optoResponses.baselineFR = baselineFR;

%---------------------------processing the data---------------------------%
disp(session)
disp('start processing...')

for i = 1:length(spike.keyboardInput)
    s.optoTrainTimes = [];
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
    
    if any(strcmp(fieldnames(spike), 'optoTrainTimes'))
        s.optoTrainTimes = spike.optoTrainTimes;
    end
    
    switch(spike.keyboardInput(i))
        case 'J'
            s.toneDuration = s.RMToneDuration;
            s.intervalDuration = s.RMIntervalDuration;
            keyboardTime = spike.keyboardTimes(i);
            toneNumber = length(freq.J);
  
            [resultMatrix{i, 1}, optoResultMatrix{i, 1}, resultMatrix{i, 2}, optoResultMatrix{i, 2}] = processAuditoryData(keyboardTime, spike.audioSignal, ...
                spike.audioSignalTimes, ephysData, ephysTime, toneNumber, s, 'J');
            
            
        case 'K'
            s.toneDuration = s.RMToneDuration;
            s.intervalDuration = s.RMIntervalDuration;
            keyboardTime = spike.keyboardTimes(i);
            toneNumber = length(freq.K);
  
            [resultMatrix{i, 1}, optoResultMatrix{i, 1}, resultMatrix{i, 2}, optoResultMatrix{i, 2}] = processAuditoryData(keyboardTime, spike.audioSignal, ...
                spike.audioSignalTimes, ephysData, ephysTime, toneNumber, s, 'K');
            
        case 'L'
            s.toneDuration = s.RMToneDuration;
            s.intervalDuration = s.RMIntervalDuration;
            keyboardTime = spike.keyboardTimes(i);
            toneNumber = length(freq.L);
  
            [resultMatrix{i, 1}, optoResultMatrix{i, 1}, resultMatrix{i, 2}, optoResultMatrix{i, 2}] = processAuditoryData(keyboardTime, spike.audioSignal, ...
                spike.audioSignalTimes, ephysData, ephysTime, toneNumber, s, 'L');
            
        case 'H'
            s.toneDuration = s.RMToneDuration;
            s.intervalDuration = s.RMIntervalDuration;
            keyboardTime = spike.keyboardTimes(i);
            toneNumber = length(freq.H);
  
            [resultMatrix{i, 1}, optoResultMatrix{i, 1}, resultMatrix{i, 2}, optoResultMatrix{i, 2}] = processAuditoryData(keyboardTime, spike.audioSignal, ...
                spike.audioSignalTimes, ephysData, ephysTime, toneNumber, s, 'H');
            
        case 'U'
            s.toneDuration = s.RMToneDuration;
            s.intervalDuration = s.RMIntervalDuration;
            keyboardTime = spike.keyboardTimes(i);
            toneNumber = length(freq.U);
  
            [resultMatrix{i, 1}, optoResultMatrix{i, 1}, resultMatrix{i, 2}, optoResultMatrix{i, 2}] = processAuditoryData(keyboardTime, spike.audioSignal, ...
                spike.audioSignalTimes, ephysData, ephysTime, toneNumber, s, 'U');
            
        case 'O'
            s.toneDuration = s.RMToneDuration;
            s.intervalDuration = s.RMIntervalDuration;
            keyboardTime = spike.keyboardTimes(i);
            toneNumber = length(freq.O);
  
            [resultMatrix{i, 1}, optoResultMatrix{i, 1}, resultMatrix{i, 2}, optoResultMatrix{i, 2}] = processAuditoryData(keyboardTime, spike.audioSignal, ...
                spike.audioSignalTimes, ephysData, ephysTime, toneNumber, s, 'O');
            
        case 'P'
            s.toneDuration = s.RMToneDuration;
            s.intervalDuration = s.RMIntervalDuration;
            keyboardTime = spike.keyboardTimes(i);
            toneNumber = length(freq.P);
  
            [resultMatrix{i, 1}, optoResultMatrix{i, 1}, resultMatrix{i, 2}, optoResultMatrix{i, 2}] = processAuditoryData(keyboardTime, spike.audioSignal, ...
                spike.audioSignalTimes, ephysData, ephysTime, toneNumber, s, 'P');
            
        case 'Y'
            s.toneDuration = s.RMToneDuration;
            s.intervalDuration = s.RMIntervalDuration;
            keyboardTime = spike.keyboardTimes(i);
            toneNumber = length(freq.Y);
  
            [resultMatrix{i, 1}, optoResultMatrix{i, 1}, resultMatrix{i, 2}, optoResultMatrix{i, 2}] = processAuditoryData(keyboardTime, spike.audioSignal, ...
                spike.audioSignalTimes, ephysData, ephysTime, toneNumber, s, 'Y');
            
        case 'B'
            s.toneDuration = s.BBNDuration;
            s.intervalDuration = s.BBNIntervalDuration;
            keyboardTime = spike.keyboardTimes(i);
            toneNumber = length(loudness.RLFBBN);
  
            [resultMatrix{i, 1}, optoResultMatrix{i, 1}, resultMatrix{i, 2}, optoResultMatrix{i, 2}] = processAuditoryData(keyboardTime, spike.audioSignal, ...
                spike.audioSignalTimes, ephysData, ephysTime, toneNumber, s, 'B');
            
        case 'E'
            s.toneDuration = s.BFDuration;
            s.intervalDuration = s.BFIntervalDuration;
            keyboardTime = spike.keyboardTimes(i);
            toneNumber = length(loudness.RLFBF);
  
            [resultMatrix{i, 1}, optoResultMatrix{i, 1}, resultMatrix{i, 2}, optoResultMatrix{i, 2}] = processAuditoryData(keyboardTime, spike.audioSignal, ...
                spike.audioSignalTimes, ephysData, ephysTime, toneNumber, s, 'E');          
    end
end

disp('Parfor loop completed! Start formatting the results!')

for i = 1:length(spike.keyboardInput)
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
    key = spike.keyboardInput(i);
    if contains('EBJKLHUOPY', key)
        str = ['FR = []; FR = resultMatrix{i, 1}; if isempty(responses.', key, '); responses.', key, ' = FR; else; responses.', key, ' = [responses.', key, '; FR]; end;' ];
        eval(str);
        str = ['FR = []; FR = resultMatrix{i, 2}; if isempty(intervalResponses.', key, '); intervalResponses.', key, ' = FR; else; intervalResponses.', key, ' = [intervalResponses.', key, '; FR]; end;' ];
        eval(str);
        str = ['FR = []; FR = optoResultMatrix{i, 1}; if isempty(optoResponses.', key, '); optoResponses.', key, ' = FR; else; optoResponses.', key, ' = [optoResponses.', key, '; FR]; end;' ];
        eval(str);
        str = ['FR = []; FR = optoResultMatrix{i, 2}; if isempty(optoIntervalResponses.', key, '); optoIntervalResponses.', key, ' = FR; else; optoIntervalResponses.', key, ' = [optoIntervalResponses.', key, '; FR]; end;' ];
        eval(str);
    end    
end

responses.ephysThresh = s.threshold;
optoResponses.ephysThresh = s.threshold;

disp('Done');

end