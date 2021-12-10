function makeRateLevelFunction(session, outputCellChannel, bestFreq, varargin)

% settings
s.toneDuration = 0.2; % unit: sec
s.intervalDuration = 0.2; % unit: sec
s.BFLoudnessStartEnd = [10, 70]; % unit: dB SPL
s.loudnessStepLength = 1; % unit: dB SPL
s.BFKey = 'E';
s.BBNKey = 'B';

if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs

% initializations
rootFolder = 'Z:\Qianyun\DCN\';
gitFolder = 'D:\DCN_Project\Github\DCN\';
sessionFolder = fullfile(rootFolder, 'Data', session);

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

% load neuralData.mat & sessionEphysInfo.mat
if ~exist(fullfile(sessionFolder, 'neuralData.mat'), 'file') || ~exist(fullfile(sessionFolder, 'sessionEphysInfo.mat'), 'file')
    error('neuralData.mat or sessioinEphysInfo.mat NOT exist!');
else
    disp('loading neuralData.mat...')
    load(fullfile(rootFolder, 'Data', session, 'sessionEphysInfo.mat'), 'sessionEphysInfo');
    mapFile = sessionEphysInfo.mapFile;
    load(fullfile('Z:\obstacleData\ephys\channelMaps\kilosort', [mapFile, '.mat']), 'channelNum_OpenEphys');
     
    neuralData = load(fullfile(sessionFolder, 'neuralData.mat'));
    ind = find(neuralData.bestChannels == outputCellChannel);
    unitID = neuralData.unit_ids(ind);
    spkTimes = neuralData.spkTimes{ind};
    bestChannel = neuralData.bestChannels(ind, :);    
end


% process rate level function for best freq and BBN
configColNames = cell({s.BFKey, s.BBNKey});
rows = 6;
config = cell(rows, length(configColNames));

configRowName = cell(rows, 1);
configRowName{1, 1} = 'keyTriggerTime';
configRowName{2, 1} = 'startInds';
configRowName{3, 1} = 'stopInds';
configRowName{4, 1} = 'loudnessSeq';
configRowName{5, 1} = 'meanFR';
configRowName{6, 1} = 'tone';

loudnessLevel = s.BFLoudnessStartEnd(1):s.loudnessStepLength:s.BFLoudnessStartEnd(2);
for k = 1:length(configColNames)
    key = configColNames{k};
    keyTime = spike.keyboardTimes(find(spike.keyboardInput == key));
    keyTime = keyTime(1);
    tempInd = find(spike.audioSignalTimes >= keyTime, 1, 'first');
    startInds = nan(length(loudnessLevel), 1);
    stopInds = nan(length(loudnessLevel), 1);
    FR = nan(length(loudnessLevel), 1);
    
    for i = 1:length(loudnessLevel)
        tempInd = int32(tempInd);
        
        startInds(i) = int32(find(spike.audioSignal(tempInd:end) >= 0.003, 1, 'first') + tempInd);
        startTime = spike.audioSignalTimes(startInds(i));
        stopInds(i) = int32(startInds(i) + s.toneDuration*s.audioFs);
        stopTime = spike.audioSignalTimes(stopInds(i));
        
        spkNum = sum(spkTimes >= startTime & spkTimes <= stopTime);
        FR(i) = spkNum/(stopTime - startTime);
        
        tempInd = stopInds(i) + s.intervalDuration/2*s.audioFs;
    end

    config{1, k} = keyTime;
    config{2, k} = startInds;
    config{3, k} = stopInds;
    config{4, k} = loudnessLevel;
    config{5, k} = FR;
    if strcmp(configColNames{k}, s.BFKey)
        config{6, k} = bestFreq;
    elseif strcmp(configColNames{k}, s.BBNKey)
        config{6, k} = 'BBN';
    end
    
    % sanity check
    figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
    plot(spike.audioSignalTimes(startInds(1):stopInds(end)), spike.audioSignal(startInds(1):stopInds(end)));
    hold on;
    for i = 1:length(loudnessLevel)
        plot([spike.audioSignalTimes(startInds(i)), spike.audioSignalTimes(startInds(i))], ...
            [-5, 5], '-r', 'LineWidth', 1.5);
        plot([spike.audioSignalTimes(stopInds(i)), spike.audioSignalTimes(stopInds(i))], ...
            [-5, 5], '-c', 'LineWidth', 1.5);
    end
    
    title(['crude tone seq - key ', configColNames{k}]);
end


% calculate baseline firing rate
baselineStartTimes(1, 1) = sessionEphysInfo.convertedEphysTimestamps(1);
baselineEndTimes(1, 1) = spike.keyboardTimes(2); 
baselineStartTimes(1, 2) = spike.keyboardTimes(end-1);
baselineEndTimes(1, 2) = sessionEphysInfo.convertedEphysTimestamps(end);
baselineFR = nan(1, 2);
for i = 1:2
    spkNum = sum(spkTimes >= baselineStartTimes(i) & spkTimes <= baselineEndTimes(i));
    baselineFR(i) = spkNum/(baselineEndTimes(i) - baselineStartTimes(i));    
end

meanBaselineFR = sum(baselineFR)/length(baselineFR);


% make rate level function plot!!
figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
colorMatrix = {'#0072BD', '#A2142F'};

for i = 1:size(config, 2)
    
    x = [config{strcmp(configRowName, 'loudnessSeq'), i}]; % db spl level
    y = [config{strcmp(configRowName, 'meanFR'), i}]; % firing rate

    plot(x, y, '-', 'LineWidth', 2.5, 'Color', colorMatrix{i}); hold on;
    box off; 
end
plot(x, repmat(meanBaselineFR, size(x)), '--k', 'LineWidth', 0.5);
xlabel('Loudness (dB SPL)');
ylabel('Firing Rate (spikes/s)');
ylim([0, 160]);
legend([num2str(bestFreq/1000), ' kHz'], 'BBN', 'Baseline');
legend boxoff
title(['Rate Level Function ', session, ' Unit ', num2str(unitID), ' Ch', num2str(outputCellChannel)],...
            'Interpreter', 'none');

saveas(gcf, fullfile(sessionFolder, 'trialFigs', ['Rate Level Function Unit', num2str(unitID), ' Ch' num2str(bestChannel), '.png']));











end