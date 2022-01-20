function makeResponseMap(session, outputCellChannel, varargin)


% settings
s.toneDuration = 0.2; % unit: sec
s.intervalDuration = 0.05; % unit: sec
s.crudeToneSeqStartEnd = [4500, 72000];
s.crudeToneSeqNum = 40;
s.crudeToneLoudnessLevel = [70, 50, 40, 30]; % unit: db spl
s.fineToneSeqStartEnd = [5000, 10000];
s.fineToneSeqNum = 40;
s.fineToneLoudnessLevel = [50, 40, 30, 20]; % unit: db spl
s.baselineFRTimeWindow = 5; % unit: sec

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
    load(fullfile(rootFolder, 'Data', session, 'sessionEphysInfo.mat'), 'sessionEphysInfo');
    mapFile = sessionEphysInfo.mapFile;
    load(fullfile('Z:\obstacleData\ephys\channelMaps\kilosort', [mapFile, '.mat']), 'channelNum_OpenEphys');
     
    neuralData = load(fullfile(sessionFolder, 'neuralData.mat'));
    ind = find(neuralData.bestChannels == outputCellChannel);
    unitID = neuralData.unit_ids(ind);
    spkTimes = neuralData.spkTimes{ind};
    bestChannel = neuralData.bestChannels(ind, :);    
end

% calculate crude and fine tone freq sequence
s.crudeToneSeq = logspace(log10(s.crudeToneSeqStartEnd(1)), log10(s.crudeToneSeqStartEnd(2)), s.crudeToneSeqNum);
s.fineToneSeq = logspace(log10(s.fineToneSeqStartEnd(1)), log10(s.fineToneSeqStartEnd(2)), s.fineToneSeqNum);

% process the crude tone seq
configColNames = cell({'J', 'K', 'L', 'H'});
rows = 7;
config = cell(rows, length(configColNames));

configRowName = cell(rows, 1);
configRowName{1, 1} = 'keyTriggerTime';
configRowName{2, 1} = 'startInds';
configRowName{3, 1} = 'stopInds';
configRowName{4, 1} = 'toneFreqSeq';
configRowName{5, 1} = 'meanFR';
configRowName{6, 1} = 'loudness';
configRowName{7, 1} = 'isCrudeToneSeq';

for k = 1:length(s.crudeToneLoudnessLevel)
    key = configColNames{k};
    keyTime = spike.keyboardTimes(find(spike.keyboardInput == key));
    
    tempInd = int32(find(spike.audioSignalTimes >= keyTime, 1, 'first'));
    startInds = nan(s.crudeToneSeqNum, 1);
    stopInds = nan(s.crudeToneSeqNum, 1);
    FR = nan(s.crudeToneSeqNum, 1);
    for i = 1:s.crudeToneSeqNum
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
    config{4, k} = s.crudeToneSeq;
    config{5, k} = FR;
    config{6, k} = s.crudeToneLoudnessLevel(k);
    config{7, k} = true;
    
    % sanity check
    figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
    plot(spike.audioSignalTimes(startInds(1):stopInds(end)), spike.audioSignal(startInds(1):stopInds(end)));
    hold on;
    for i = 1:s.crudeToneSeqNum
        plot([spike.audioSignalTimes(startInds(i)), spike.audioSignalTimes(startInds(i))], ...
            [-5, 5], '-r', 'LineWidth', 1.5);
        plot([spike.audioSignalTimes(stopInds(i)), spike.audioSignalTimes(stopInds(i))], ...
            [-5, 5], '-c', 'LineWidth', 1.5);
    end
    
    title(['crude tone seq - key ', configColNames{k}]);
end


% process the fine tone seq    
configColNames = [configColNames, {'U', 'O', 'P', 'Y'}];
tempColCount = size(config, 2);
for k = 1:length(s.fineToneLoudnessLevel)
    key = configColNames{tempColCount+k};
    keyTime = spike.keyboardTimes(find(spike.keyboardInput == key));
    
    tempInd = find(spike.audioSignalTimes >= keyTime, 1, 'first');
    startInds = nan(s.fineToneSeqNum, 1);
    stopInds = nan(s.fineToneSeqNum, 1);
    FR = nan(s.fineToneSeqNum, 1);
    for i = 1:s.fineToneSeqNum
        tempInd = int32(tempInd);
        
        startInds(i) = int32(find(spike.audioSignal(tempInd:end) >= 0.003, 1, 'first') + tempInd);
        startTime = spike.audioSignalTimes(startInds(i));
        stopInds(i) = int32(startInds(i) + s.toneDuration*s.audioFs);
        stopTime = spike.audioSignalTimes(stopInds(i));
        
        spkNum = sum(spkTimes >= startTime & spkTimes <= stopTime);
        FR(i) = spkNum/(stopTime - startTime);
        
        tempInd = stopInds(i) + s.intervalDuration/2*s.audioFs;
    end
    
    config{1, tempColCount+k} = keyTime;
    config{2, tempColCount+k} = startInds;
    config{3, tempColCount+k} = stopInds;
    config{4, tempColCount+k} = s.fineToneSeq;
    config{5, tempColCount+k} = FR;
    config{6, tempColCount+k} = s.fineToneLoudnessLevel(k);
    config{7, tempColCount+k} = false;
    
    % sanity check
    figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
    plot(spike.audioSignalTimes(startInds(1):stopInds(end)), spike.audioSignal(startInds(1):stopInds(end)));
    hold on;
    for i = 1:s.crudeToneSeqNum
        plot([spike.audioSignalTimes(startInds(i)), spike.audioSignalTimes(startInds(i))], ...
            [-5, 5], '-r', 'LineWidth', 1.5);
        plot([spike.audioSignalTimes(stopInds(i)), spike.audioSignalTimes(stopInds(i))], ...
            [-5, 5], '-c', 'LineWidth', 1.5);
    end
    
    title(['fine tone seq - key ', configColNames{tempColCount+k}]);
end    

% combine crude tone seq and fine tone seq together!
loudnessRow = find(strcmp(configRowName, 'loudness'));
loudnessLevel = unique([config{loudnessRow, :}]);

responseMap = cell(3, length(loudnessLevel));
for i = 1:length(loudnessLevel)
    responseMap{1, i} = loudnessLevel(i);
    if sum([config{loudnessRow, :}] == loudnessLevel(i)) == 1
        ind = find([config{loudnessRow, :}] == loudnessLevel(i));
        responseMap{2, i} = config{strcmp(configRowName, 'toneFreqSeq'), ind};
        responseMap{3, i} = config{strcmp(configRowName, 'meanFR'), ind};
    else
        inds = find([config{loudnessRow, :}] == loudnessLevel(i));
        toneSeq = [];
        FR = [];
        for j = 1:length(inds)
            toneSeq = [toneSeq, [config{strcmp(configRowName, 'toneFreqSeq'), inds(j)}]];
            FR = [FR, [config{strcmp(configRowName, 'meanFR'), inds(j)}]'];
        end
        [sortedToneSeq, index] = sort(toneSeq);
        sortedFR = FR(index);
        responseMap{2, i} = sortedToneSeq;
        responseMap{3, i} = sortedFR;
    end   
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


% plot response map!!
figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
colorMatrix = {'#0072BD', 	'#EDB120', 	'#77AC30', 	'#7E2F8E', 	'#A2142F'};
rows = length(loudnessLevel);
cols = 1;
plotInd = length(loudnessLevel);

for i = 1:length(loudnessLevel)
    subplot(rows, cols, plotInd);
    
    x = [responseMap{2, i}]/1000;
    y = [responseMap{3, i}];

    semilogx(x, y, '-', 'LineWidth', 1.5, 'Color', colorMatrix{i}); hold on;
    plot(x, repmat(meanBaselineFR, size(x)), '--', 'LineWidth', 0.5, 'Color', colorMatrix{i});
    box off; 
    legend([num2str(responseMap{1, i}) ' dB SPL']);
    legend boxoff;
    
    if i ~= 1
        axis tight;
        h = gca;
        h.YAxis.Visible = 'off';
        h.XAxis.Visible = 'off';
        
    else
        xticks([5 10 20 50]);
        xlim([s.crudeToneSeqStartEnd(1)/1000, s.crudeToneSeqStartEnd(2)/1000]);
        xlabel('frequency (kHz)');
        ylabel('spikes/s');
    end
    
        
    if i == length(loudnessLevel)
        title(['Response Map ', session, ' Unit ', num2str(unitID), ' Ch', num2str(outputCellChannel)],...
            'Interpreter', 'none');
    end
    
    plotInd = plotInd - 1;
    
end
saveas(gcf, fullfile(sessionFolder, 'trialFigs', ['Reponse Map Unit', num2str(unitID), ' Ch' num2str(bestChannel), '.png']));

end