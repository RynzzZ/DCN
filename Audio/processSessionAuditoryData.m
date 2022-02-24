function [responses, baselineFR] = processSessionAuditoryData(session, freq, loudness, ephysChannel, responses, varargin)

% settings
s.toneDuration = 0.2; % unit: sec
s.RLFBBNDuration = 0.2; 
s.RLFBBNIntervalDuration = 0.2;
s.RLFBFDuration = 0.2;
s.RLFBFIntervalDuration = 0.2;
s.intervalDuration = 0.05;
s.threshold = -300;
s.thresholdType = 'thresholdDown';
s.plotType = '';

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
end

% get baseline firing rate
ephysTime = sessionEphysInfo.convertedEphysTimestamps;
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

ind = strfind(convertCharsToStrings(spike.keyboardInput), 'J');
ind = ind(1);
startTime = spike.keyboardTimes(ind);
endTime = startTime + 10;

startInd = find(sessionEphysInfo.convertedEphysTimestamps >= startTime, 1, 'first');
endInd = find(sessionEphysInfo.convertedEphysTimestamps <= endTime, 1, 'last');

x = sessionEphysInfo.convertedEphysTimestamps(startInd:endInd);
y = getVoltage(data.Data.Data(channelNum_OpenEphys(ephysChannel), startInd:endInd));
figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
plot(x, y);

%---------------------------processing the data---------------------------%

disp('start processing...')
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

    disp(spike.keyboardInput(i));
    
    switch(spike.keyboardInput(i))
        case 'J'
            keyboardTime = spike.keyboardTimes(i);
            tempInd = int32(find(spike.audioSignalTimes >= keyboardTime, 1, 'first'));
            
            spikeStartInds = nan(length(freq.J), 1);
            ephysStartInds = nan(length(freq.J), 1);
            spikeStopInds = nan(length(freq.J), 1);
            ephysStopInds = nan(length(freq.J), 1);
            FR = nan(length(freq.J), 1);
            
            for j = 1:length(freq.J)
                if mod(j, 10) == 1
                    fprintf('%f ', j/length(freq.J))
                end
                
                spikeStartInds(j) = int32(find(spike.audioSignal(tempInd:end) >= 0.003, 1, 'first') + tempInd);
                startTime = spike.audioSignalTimes(spikeStartInds(j));
                ephysStartInds(j) = find(sessionEphysInfo.convertedEphysTimestamps >= startTime, 1, 'first');
                
                spikeStopInds(j) = int32(spikeStartInds(j) + s.toneDuration*s.audioFs);
                stopTime = spike.audioSignalTimes(spikeStopInds(j));
                ephysStopInds(j) = find(sessionEphysInfo.convertedEphysTimestamps <= stopTime, 1, 'last');
                
                
                chunkEphysData = getVoltage(data.Data.Data(channelNum_OpenEphys(ephysChannel), ephysStartInds(j):ephysStopInds(j)));
                [~, locs, ncrs] = crossdet(chunkEphysData, s.threshold, s.thresholdType);
                % spkNum = ncrs - sum(diff(locs)/sessionEphysInfo.fs <= 0.006);
                disp(ncrs);
                FR(j) = ncrs/(stopTime - startTime);
              
                tempInd = int32(spikeStopInds(j) + s.intervalDuration/2*s.audioFs);
                
                % quality check
                if mod(j, 8) == 0
                    figure;
                    x = linspace(startTime, stopTime, ephysStopInds(j) - ephysStartInds(j) + 1);
                    plot(x, chunkEphysData); axis tight; box off; hold on;
                    y = ones(size(locs))*s.threshold;
                    scatter(x(locs), y, '.', 'r');
                    title(['frequency = ', num2str(freq.J(j)/1000), 'kHz']);
                end
            end
                       
            if size(FR, 2) == 1; FR = FR'; end
            if isempty(responses.J)
                responses.J = FR;
            else
                responses.J = [responses.J; FR];
            end
            
            figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
            plot(freq.(spike.keyboardInput(i)), FR);
            title(['Keyboard Input ', spike.keyboardInput(i)]);
            
            
        case 'K'
            keyboardTime = spike.keyboardTimes(i);
            tempInd = int32(find(spike.audioSignalTimes >= keyboardTime, 1, 'first'));
            
            spikeStartInds = nan(length(freq.K), 1);
            ephysStartInds = nan(length(freq.K), 1);
            spikeStopInds = nan(length(freq.K), 1);
            ephysStopInds = nan(length(freq.K), 1);
            FR = nan(length(freq.K), 1);
            
            for j = 1:length(freq.K)
                if mod(j, 10) == 1
                    fprintf('%f ', j/length(freq.K))
                end
                spikeStartInds(j) = int32(find(spike.audioSignal(tempInd:end) >= 0.003, 1, 'first') + tempInd);
                startTime = spike.audioSignalTimes(spikeStartInds(j));
                ephysStartInds(j) = find(sessionEphysInfo.convertedEphysTimestamps >= startTime, 1, 'first');
                
                spikeStopInds(j) = int32(spikeStartInds(j) + s.toneDuration*s.audioFs);
                stopTime = spike.audioSignalTimes(spikeStopInds(j));
                ephysStopInds(j) = find(sessionEphysInfo.convertedEphysTimestamps <= stopTime, 1, 'last');
                
                
                chunkEphysData = getVoltage(data.Data.Data(channelNum_OpenEphys(ephysChannel), ephysStartInds(j):ephysStopInds(j)));
                [~, locs, ncrs] = crossdet(chunkEphysData, s.threshold, s.thresholdType);
                % spkNum = ncrs - sum(diff(locs)/sessionEphysInfo.fs <= 0.006);
                disp(ncrs);
                FR(j) = ncrs/(stopTime - startTime);
              
                tempInd = int32(spikeStopInds(j) + s.intervalDuration/2*s.audioFs);
            end
            
            if size(FR, 2) == 1; FR = FR'; end
            if isempty(responses.K)
                responses.K = FR;
            else
                responses.K = [responses.K; FR];
            end
            
            figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
            plot(freq.(spike.keyboardInput(i)), FR);
            title(['Keyboard Input ', spike.keyboardInput(i)]);
            
        case 'L'

            keyboardTime = spike.keyboardTimes(i);
            tempInd = int32(find(spike.audioSignalTimes >= keyboardTime, 1, 'first'));
            
            spikeStartInds = nan(length(freq.L), 1);
            ephysStartInds = nan(length(freq.L), 1);
            spikeStopInds = nan(length(freq.L), 1);
            ephysStopInds = nan(length(freq.L), 1);
            FR = nan(length(freq.L), 1);
            
            for j = 1:length(freq.L)
                if mod(j, 10) == 1
                    fprintf('%f ', j/length(freq.L))
                end
                spikeStartInds(j) = int32(find(spike.audioSignal(tempInd:end) >= 0.003, 1, 'first')) + tempInd;
                startTime = spike.audioSignalTimes(spikeStartInds(j));
                ephysStartInds(j) = find(sessionEphysInfo.convertedEphysTimestamps >= startTime, 1, 'first');
                
                spikeStopInds(j) = int32(spikeStartInds(j) + s.toneDuration*s.audioFs);
                stopTime = spike.audioSignalTimes(spikeStopInds(j));
                ephysStopInds(j) = find(sessionEphysInfo.convertedEphysTimestamps <= stopTime, 1, 'last');
                
                
                chunkEphysData = getVoltage(data.Data.Data(channelNum_OpenEphys(ephysChannel), ephysStartInds(j):ephysStopInds(j)));
                [~, locs, ncrs] = crossdet(chunkEphysData, s.threshold, s.thresholdType);
                % spkNum = ncrs - sum(diff(locs)/sessionEphysInfo.fs <= 0.006);
                disp(ncrs);
                FR(j) = ncrs/(stopTime - startTime);
              
                tempInd = int32(spikeStopInds(j) + s.intervalDuration/2*s.audioFs);
            end
            
            if size(FR, 2) == 1; FR = FR'; end
            if isempty(responses.L)
                responses.L = FR;
            else
                responses.L = [responses.L; FR];
            end
            
            figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
            plot(freq.(spike.keyboardInput(i)), FR);
            title(['Keyboard Input ', spike.keyboardInput(i)]);
            
        case 'H'

            keyboardTime = spike.keyboardTimes(i);
            tempInd = int32(find(spike.audioSignalTimes >= keyboardTime, 1, 'first'));
            
            spikeStartInds = nan(length(freq.H), 1);
            ephysStartInds = nan(length(freq.H), 1);
            spikeStopInds = nan(length(freq.H), 1);
            ephysStopInds = nan(length(freq.H), 1);
            FR = nan(length(freq.H), 1);
            
            for j = 1:length(freq.H)
                if mod(j, 10) == 1
                    fprintf('%f ', j/length(freq.H))
                end

                spikeStartInds(j) = int32(find(spike.audioSignal(tempInd:end) >= 0.003, 1, 'first') + int32(tempInd));
                startTime = spike.audioSignalTimes(spikeStartInds(j));
                ephysStartInds(j) = find(sessionEphysInfo.convertedEphysTimestamps >= startTime, 1, 'first');
                
                spikeStopInds(j) = int32(spikeStartInds(j) + s.toneDuration*s.audioFs);
                stopTime = spike.audioSignalTimes(spikeStopInds(j));
                ephysStopInds(j) = find(sessionEphysInfo.convertedEphysTimestamps <= stopTime, 1, 'last');
                
                
                chunkEphysData = getVoltage(data.Data.Data(channelNum_OpenEphys(ephysChannel), ephysStartInds(j):ephysStopInds(j)));
                [~, locs, ncrs] = crossdet(chunkEphysData, s.threshold, s.thresholdType);
                % spkNum = ncrs - sum(diff(locs)/sessionEphysInfo.fs <= 0.006);
                disp(ncrs);
                FR(j) = ncrs/(stopTime - startTime);
              
                tempInd = int32(spikeStopInds(j) + s.intervalDuration/2*s.audioFs);
            end
            

            if size(FR, 2) == 1; FR = FR'; end
            if isempty(responses.H)
                responses.H = FR;
            else
                responses.H = [responses.H; FR];
            end
            
            figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
            plot(freq.(spike.keyboardInput(i)), FR);
            title(['Keyboard Input ', spike.keyboardInput(i)]);
            
        case 'U'

            keyboardTime = spike.keyboardTimes(i);
            tempInd = int32(find(spike.audioSignalTimes >= keyboardTime, 1, 'first'));
            
            spikeStartInds = nan(length(freq.U), 1);
            ephysStartInds = nan(length(freq.U), 1);
            spikeStopInds = nan(length(freq.U), 1);
            ephysStopInds = nan(length(freq.U), 1);
            FR = nan(length(freq.U), 1);
            
            for j = 1:length(freq.U)
                if mod(j, 10) == 1
                    fprintf('%f ', j/length(freq.U))
                end
                spikeStartInds(j) = int32(find(spike.audioSignal(tempInd:end) >= 0.003, 1, 'first') + tempInd);
                startTime = spike.audioSignalTimes(spikeStartInds(j));
                ephysStartInds(j) = find(sessionEphysInfo.convertedEphysTimestamps >= startTime, 1, 'first');
                
                spikeStopInds(j) = int32(spikeStartInds(j) + s.toneDuration*s.audioFs);
                stopTime = spike.audioSignalTimes(spikeStopInds(j));
                ephysStopInds(j) = find(sessionEphysInfo.convertedEphysTimestamps <= stopTime, 1, 'last');
                
                
                chunkEphysData = getVoltage(data.Data.Data(channelNum_OpenEphys(ephysChannel), ephysStartInds(j):ephysStopInds(j)));
                [~, locs, ncrs] = crossdet(chunkEphysData, s.threshold, s.thresholdType);
                % spkNum = ncrs - sum(diff(locs)/sessionEphysInfo.fs <= 0.006);
                disp(ncrs);
                FR(j) = ncrs/(stopTime - startTime);
              
                tempInd = int32(spikeStopInds(j) + s.intervalDuration/2*s.audioFs);
            end
            
            if size(FR, 2) == 1; FR = FR'; end
            if isempty(responses.U)
                responses.U = FR;
            else
                responses.U = [responses.U; FR];
            end
            
            figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
            plot(freq.(spike.keyboardInput(i)), FR);
            title(['Keyboard Input ', spike.keyboardInput(i)]);
            
        case 'O'

            keyboardTime = spike.keyboardTimes(i);
            tempInd = int32(find(spike.audioSignalTimes >= keyboardTime, 1, 'first'));
            
            spikeStartInds = nan(length(freq.O), 1);
            ephysStartInds = nan(length(freq.O), 1);
            spikeStopInds = nan(length(freq.O), 1);
            ephysStopInds = nan(length(freq.O), 1);
            FR = nan(length(freq.O), 1);
            
            for j = 1:length(freq.O)
                if mod(j, 10) == 1
                    fprintf('%f ', j/length(freq.O))
                end
                spikeStartInds(j) = int32(find(spike.audioSignal(tempInd:end) >= 0.003, 1, 'first') + tempInd);
                startTime = spike.audioSignalTimes(spikeStartInds(j));
                ephysStartInds(j) = find(sessionEphysInfo.convertedEphysTimestamps >= startTime, 1, 'first');
                
                spikeStopInds(j) = int32(spikeStartInds(j) + s.toneDuration*s.audioFs);
                stopTime = spike.audioSignalTimes(spikeStopInds(j));
                ephysStopInds(j) = find(sessionEphysInfo.convertedEphysTimestamps <= stopTime, 1, 'last');
                
                
                chunkEphysData = getVoltage(data.Data.Data(channelNum_OpenEphys(ephysChannel), ephysStartInds(j):ephysStopInds(j)));
                [~, locs, ncrs] = crossdet(chunkEphysData, s.threshold, s.thresholdType);
                % spkNum = ncrs - sum(diff(locs)/sessionEphysInfo.fs <= 0.006);
                disp(ncrs);
                FR(j) = ncrs/(stopTime - startTime);
              
                tempInd = int32(spikeStopInds(j) + s.intervalDuration/2*s.audioFs);
            end
            
            if size(FR, 2) == 1; FR = FR'; end
            if isempty(responses.O)
                responses.O = FR;
            else
                responses.O = [responses.O; FR];
            end
            
            figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
            plot(freq.(spike.keyboardInput(i)), FR);
            title(['Keyboard Input ', spike.keyboardInput(i)]);
            
        case 'P'

            keyboardTime = spike.keyboardTimes(i);
            tempInd = int32(find(spike.audioSignalTimes >= keyboardTime, 1, 'first'));
            
            spikeStartInds = nan(length(freq.P), 1);
            ephysStartInds = nan(length(freq.P), 1);
            spikeStopInds = nan(length(freq.P), 1);
            ephysStopInds = nan(length(freq.P), 1);
            FR = nan(length(freq.P), 1);
            
            for j = 1:length(freq.P)
                if mod(j, 10) == 1
                    fprintf('%f ', j/length(freq.P))
                end
                spikeStartInds(j) = int32(find(spike.audioSignal(tempInd:end) >= 0.003, 1, 'first') + tempInd);
                startTime = spike.audioSignalTimes(spikeStartInds(j));
                ephysStartInds(j) = find(sessionEphysInfo.convertedEphysTimestamps >= startTime, 1, 'first');
                
                spikeStopInds(j) = int32(spikeStartInds(j) + s.toneDuration*s.audioFs);
                stopTime = spike.audioSignalTimes(spikeStopInds(j));
                ephysStopInds(j) = find(sessionEphysInfo.convertedEphysTimestamps <= stopTime, 1, 'last');
                
                
                chunkEphysData = getVoltage(data.Data.Data(channelNum_OpenEphys(ephysChannel), ephysStartInds(j):ephysStopInds(j)));
                [~, locs, ncrs] = crossdet(chunkEphysData, s.threshold, s.thresholdType);
                % spkNum = ncrs - sum(diff(locs)/sessionEphysInfo.fs <= 0.006);
                disp(ncrs);
                FR(j) = ncrs/(stopTime - startTime);
              
                tempInd = int32(spikeStopInds(j) + s.intervalDuration/2*s.audioFs);
            end
            
            if size(FR, 2) == 1; FR = FR'; end
            if isempty(responses.P)
                responses.P = FR;
            else
                responses.P = [responses.P; FR];
            end
            
            figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
            plot(freq.(spike.keyboardInput(i)), FR);
            title(['Keyboard Input ', spike.keyboardInput(i)]);
            
        case 'Y'

            keyboardTime = spike.keyboardTimes(i);
            tempInd = int32(find(spike.audioSignalTimes >= keyboardTime, 1, 'first'));
            
            spikeStartInds = nan(length(freq.Y), 1);
            ephysStartInds = nan(length(freq.Y), 1);
            spikeStopInds = nan(length(freq.Y), 1);
            ephysStopInds = nan(length(freq.Y), 1);
            FR = nan(length(freq.Y), 1);
            
            for j = 1:length(freq.Y)
                if mod(j, 10) == 1
                    fprintf('%f ', j/length(freq.Y))
                end
                spikeStartInds(j) = int32(find(spike.audioSignal(tempInd:end) >= 0.003, 1, 'first') + tempInd);
                startTime = spike.audioSignalTimes(spikeStartInds(j));
                ephysStartInds(j) = find(sessionEphysInfo.convertedEphysTimestamps >= startTime, 1, 'first');
                
                spikeStopInds(j) = int32(spikeStartInds(j) + s.toneDuration*s.audioFs);
                stopTime = spike.audioSignalTimes(spikeStopInds(j));
                ephysStopInds(j) = find(sessionEphysInfo.convertedEphysTimestamps <= stopTime, 1, 'last');
                
                
                chunkEphysData = getVoltage(data.Data.Data(channelNum_OpenEphys(ephysChannel), ephysStartInds(j):ephysStopInds(j)));
                [~, locs, ncrs] = crossdet(chunkEphysData, s.threshold, s.thresholdType);
                % spkNum = ncrs - sum(diff(locs)/sessionEphysInfo.fs <= 0.006);
                disp(ncrs);
                FR(j) = ncrs/(stopTime - startTime);
              
                tempInd = int32(spikeStopInds(j) + s.intervalDuration/2*s.audioFs);
            end
            
            if size(FR, 2) == 1; FR = FR'; end
            if isempty(responses.Y)
                responses.Y = FR;
            else
                responses.Y = [responses.Y; FR];
            end
            
            figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
            plot(freq.(spike.keyboardInput(i)), FR);
            title(['Keyboard Input ', spike.keyboardInput(i)]);
            
        case 'B'
            keyboardTime = spike.keyboardTimes(i);
            tempInd = int32(find(spike.audioSignalTimes >= keyboardTime, 1, 'first'));
            
            for j = 1:length(loudness.RLFBBN)
                if mod(j, 10) == 1
                    fprintf('%f ', j/length(loudness.RLFBBN))
                end
                spikeStartInds(j) = int32(find(spike.audioSignal(tempInd:end) >= 0.002, 1, 'first') + tempInd);
                startTime = spike.audioSignalTimes(spikeStartInds(j));
                ephysStartInds(j) = find(sessionEphysInfo.convertedEphysTimestamps >= startTime, 1, 'first');
                
                spikeStopInds(j) = int32(spikeStartInds(j) + s.RLFBBNDuration*s.audioFs);
                stopTime = spike.audioSignalTimes(spikeStopInds(j));
                ephysStopInds(j) = find(sessionEphysInfo.convertedEphysTimestamps <= stopTime, 1, 'last');
                
                
                chunkEphysData = getVoltage(data.Data.Data(channelNum_OpenEphys(ephysChannel), ephysStartInds(j):ephysStopInds(j)));
                [~, locs, ncrs] = crossdet(chunkEphysData, s.threshold, s.thresholdType);
                % spkNum = ncrs - sum(diff(locs)/sessionEphysInfo.fs <= 0.006);
                disp(ncrs);
                FR(j) = ncrs/(stopTime - startTime);
              
                tempInd = int32(spikeStopInds(j) + s.RLFBBNIntervalDuration/2*s.audioFs);
            end
            
            if size(FR, 2) == 1; FR = FR'; end
            if isempty(responses.B)
                responses.B = FR;
            else
                responses.B = [responses.B; FR];
            end
            
            figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
            plot(loudness.RLFBBN, FR);
            title(['Keyboard Input ', spike.keyboardInput(i)]);
            
            
        case 'E'
            keyboardTime = spike.keyboardTimes(i);
            tempInd = int32(find(spike.audioSignalTimes >= keyboardTime, 1, 'first'));
            
            for j = 1:length(loudness.RLFBF)
                if mod(j, 10) == 1
                    fprintf('%f ', j/length(loudness.RLFBF))
                end
                spikeStartInds(j) = int32(find(spike.audioSignal(tempInd:end) >= 0.002, 1, 'first') + tempInd);
                startTime = spike.audioSignalTimes(spikeStartInds(j));
                ephysStartInds(j) = find(sessionEphysInfo.convertedEphysTimestamps >= startTime, 1, 'first');
                
                spikeStopInds(j) = int32(spikeStartInds(j) + s.RLFBFDuration*s.audioFs);
                stopTime = spike.audioSignalTimes(spikeStopInds(j));
                ephysStopInds(j) = find(sessionEphysInfo.convertedEphysTimestamps <= stopTime, 1, 'last');
                
                
                chunkEphysData = getVoltage(data.Data.Data(channelNum_OpenEphys(ephysChannel), ephysStartInds(j):ephysStopInds(j)));
                [~, locs, ncrs] = crossdet(chunkEphysData, s.threshold, s.thresholdType);
                % spkNum = ncrs - sum(diff(locs)/sessionEphysInfo.fs <= 0.006);
                disp(ncrs);
                FR(j) = ncrs/(stopTime - startTime);
                
                tempInd = int32(spikeStopInds(j) + s.RLFBFIntervalDuration/2*s.audioFs);
            end
            
            if size(FR, 2) == 1; FR = FR'; end
            if isempty(responses.E)
                responses.E = FR;
            else
                responses.E = [responses.E; FR];
            end
            
            figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
            plot(loudness.RLFBF, FR);
            title(['Keyboard Input ', spike.keyboardInput(i)]);
            
    end
end


%
%             startInd = find(sessionEphysInfo.convertedEphysTimestamps >= 7.5, 1, 'first');
%             endInd = find(sessionEphysInfo.convertedEphysTimestamps <= 68.5, 1, 'last');
%             x = sessionEphysInfo.convertedEphysTimestamps(startInd:endInd);
%             y = data.Data.Data(channelNum_OpenEphys(ephysChannel), startInd:endInd);
%             figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
%             plot(x, y); hold on;
%             for i = 1:length(ephysStartInds)
%                 x1 = sessionEphysInfo.convertedEphysTimestamps(ephysStartInds(i));
%                 plot([x1, x1], [min(y), max(y)], 'r-');
%                 x2 = sessionEphysInfo.convertedEphysTimestamps(ephysStopInds(i));
%                 plot([x2, x2], [min(y), max(y)], 'c-');
%             end



end