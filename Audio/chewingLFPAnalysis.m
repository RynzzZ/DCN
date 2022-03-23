%% settings
session = '20220302_000';
ephysChannel = 24;
audioFs = 150000;
ephysFs = 30000;
suppressFig = 0;

% audio settings
load('C:\Users\Qianyun Zhang\OneDrive\AudioFiles\Mimic_BBN_Noise\ChewMimics_2msDuration_5to70DBSPL_1DBSPLStep_RandomLevels.mat');
randomLevel1 = randomLevel;
load('C:\Users\Qianyun Zhang\OneDrive\AudioFiles\Mimic_BBN_Noise\ChewMimics_10msDuration_5to70DBSPL_1DBSPLStep_RandomLevels.mat');
randomLevel2 = randomLevel;
intervalDuration = 0.2;

% initializations
rootFolder = 'Z:\Qianyun\DCN\';
gitFolder = 'D:\DCN_Project\Github\DCN\';
sessionFolder = fullfile(rootFolder, 'Data', session);

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%loading data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% load cameraAnalyzed.mat
if ~exist(fullfile(sessionFolder, 'cameraAnalyzed.mat'), 'file')
    error('cameraAnalyzed.mat NOT exist! Run analyzeSession.m first!');
else
    disp('loading cameraAnalyzed.mat...')
    load(fullfile(sessionFolder, 'cameraAnalyzed.mat'));
    if ~any(strcmp('jawDistance', videoTracking.Properties.VariableNames))
        warning('JawDistance does NOT exist in videoTracking! Will NOT analyze chewing LFP!');
    end
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
    ephysTime = sessionEphysInfo.convertedEphysTimestamps;
    
    LFPVoltage = getVoltage(data.Data.Data(channelNum_OpenEphys(ephysChannel), :));
    rmsv = sqrt(movmean(LFPVoltage.^2, 100));  
end



%% process the data

key = '1';
keyInds = strfind(convertCharsToStrings(spike.keyboardInput), key);
keyToUse = [0, 1, 1, 0, 0, 0,  1, 1, 1, 1];
% keyToUse = [0, 1, 1, 1];

if length(keyInds) ~= length(keyToUse)
    error('Length Not Match: KeyInds & keyToUse');
end

if strcmp(key, '1')
    loudnessLevel = randomLevel1;
    BBNDuration = 0.002;
elseif strcmp(key, '2')
    loudnessLevel = randomLevel2;
    BBNDuration = 0.01;
end

[sortedLoudness, sortInds] = sort(loudnessLevel);

LFPMAX = nan(length(keyInds), length(loudnessLevel));
maxInds = nan(length(keyInds), length(loudnessLevel));
LFPMIN = nan(length(keyInds), 1);
sortedLFPMAX = nan(length(keyInds), length(loudnessLevel));
sortedLFPMIN = nan(length(keyInds), length(loudnessLevel));

for i = 1:length(keyInds)
    
    
    if keyToUse(i) == 0
        continue;
    end
    
    ind = keyInds(i);
    keyboardTime = spike.keyboardTimes(ind);
    tempInd = int32(find(spike.audioSignalTimes >= keyboardTime, 1, 'first'));
    
    startTime = nan(length(loudnessLevel), 1);
    stopTime = nan(length(loudnessLevel), 1);
    spikeStartInd = nan(length(loudnessLevel), 1);
    ephysStartInd = nan(length(loudnessLevel), 1);
    spikeStopInd = nan(length(loudnessLevel), 1);
    ephysStopInd = nan(length(loudnessLevel), 1);
    
    % process the big crunch
    for j = 1:length(loudnessLevel)
        if mod(j, 10) == 1
            fprintf('%f \n', j/length(loudnessLevel))
        end
        spikeStartInd(j) = int32(find(spike.audioSignal(tempInd:end) >= 0.0015, 1, 'first') + tempInd - BBNDuration*1*audioFs);
        startTime(j) = spike.audioSignalTimes(spikeStartInd(j));
        ephysStartInd(j) = find(ephysTime >= startTime(j), 1, 'first');
        
        spikeStopInd(j) = int32(spikeStartInd(j) + BBNDuration*3*audioFs);
        stopTime(j) = spike.audioSignalTimes(spikeStopInd(j));
        ephysStopInd(j) = find(ephysTime <= stopTime(j), 1, 'last');
        
        chunkLFPData = rmsv(ephysStartInd(j):ephysStopInd(j));
        [LFPMAX(i, j), maxInd] = max(chunkLFPData);
        maxInds(i, j) = maxInd + ephysStartInd(j);
              
        tempInd = int32(spikeStopInd(j) + intervalDuration/2*audioFs);
    end
    
    ephysTempInd = find(ephysTime >= keyboardTime, 1, 'first');
    tempInd = int32(find(spike.audioSignalTimes >= keyboardTime, 1, 'first'));
    temp = LFPMAX(i, :);
    sortedLFPMAX(i, :) = temp(sortInds);
    LFPMIN(i, 1) = mean(rmsv(ephysTempInd : ephysTempInd + ephysFs*18));
    
    % plot the big crunch sound for every trial
    if ~suppressFig
        figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
        subplot(2, 1, 1);
        plot(ephysTime(ephysTempInd : ephysTempInd + ephysFs*18), rmsv(ephysTempInd : ephysTempInd + ephysFs*18)); 
        axis tight; box off; hold on;
        scatter(ephysTime(maxInds(i, :)), LFPMAX(i, :), '.r');
%         for j = 1:length(spikeStartInd)
%             x1 = ephysTime(ephysStartInd(j)); hold on
%             plot([x1, x1], [20, 150], '-r');
%             x2 = ephysTime(ephysStopInd(j));
%             plot([x2, x2], [20, 150], '-c');
%         end
        subplot(2, 1, 2);
        plot(spike.audioSignalTimes(tempInd:tempInd + audioFs*18), spike.audioSignal(tempInd:tempInd + audioFs*18)); 
        axis tight;  box off;
        hold on;
        for j = 1:length(spikeStartInd)
            x1 = spike.audioSignalTimes(spikeStartInd(j)); hold on
            plot([x1, x1], [-0.5, 0.5], '-r');
            x2 = spike.audioSignalTimes(spikeStopInd(j));
            plot([x2, x2], [-0.5, 0.5], '-c');
        end
    end
    

    
end

%% make the chewing mimic RLF figure!

x = sortedLoudness;
figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
shadedErrorBar(x, nanmean(sortedLFPMAX), nanstd(sortedLFPMAX), ...
    'lineprops', {'linewidth', 3, 'color', '#0072BD'}, 'patchSaturation', .1);
axis tight;
box off;
xlabel('loudness (dB SPL)');
ylabel('LFP (rmsv)');
title('2ms duration BBN RLF (LFP)');
ylim([0, 200]);



figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
subplot(1, 3, 1:2)
x = sortedLoudness;
shadedErrorBar(x, nanmean(sortedLFPMAX), nanstd(sortedLFPMAX), ...
    'lineprops', {'linewidth', 3, 'color', '#0072BD'}, 'patchSaturation', .1);
axis tight;
box off;
xlabel('loudness (dB SPL)');
ylabel('LFP (rmsv)');
title('2ms duration BBN RLF (LFP)');
ylim([0, 300]);

subplot(1, 3, 3)
boxplot(crunchLFPMAX);
hold on;
ylim([0, 300]); box off;
scatter(ones(size(crunchLFPMAX)).*(1+(rand(size(crunchLFPMAX))-0.5)/10), crunchLFPMAX,'r','filled')





%% processing the actual chewing LFP!!
crunchSearchTimeWindow = 3; % sec
crunchTimeWindow = 0.02; % sec
chewingSearchStartTimeShift = 4; % sec
chewingSearchTimeWindow = 12; % sec
chewingTimeWindow = 4; % sec
chewingSearchStepLength = 0.1; % sec

trialTotalCount = spike.totalFoodNum;
trialStartTimes = nan(trialTotalCount, 1);
crunchLFPMAX = nan(trialTotalCount, 1);
crunchMaxInds = nan(trialTotalCount, 1);
chewingWindow = [];
crunchMicRMSV = cell(1, 1);
chewingJawDistance = cell(1, 1);
chewingJawTime = cell(1, 1);
chewingJawInds = cell(1, 1);
selectedTrialNum = 1;
for i = 1:trialTotalCount

    fprintf('trial %d/%d \n', i, trialTotalCount);
    
    %%%%%%%%%%%%%%%%%%%%% Processing crunch %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % locate big crunch chunk
    crunchSearchStartTime = spike.foodTriggerTimes(i) + 1;
    crunchSearchEndTime = crunchSearchStartTime + crunchSearchTimeWindow;
    trialStartTimes(i) = crunchSearchStartTime;
    
    % crunch chunk microphone
    crunchChunkMicStartInd = find(spike.micSignalTimes >= crunchSearchStartTime, 1, 'first');
    crunchChunkMicEndInd = find(spike.micSignalTimes <= crunchSearchEndTime, 1, 'last');
    crunchChunkMic = spike.micSignal(crunchChunkMicStartInd : crunchChunkMicEndInd);
    crunchChunkMic = highpass(crunchChunkMic, 100, 150000);
    crunchChunkMicRMSV = sqrt(movmean(crunchChunkMic.^2, 100));
    
    % crunch chunk ephys(LFP)
    crunchChunkEphysStartInd = find(ephysTime >= crunchSearchStartTime, 1, 'first');
    crunchChunkEphysEndInd = find(ephysTime <= crunchSearchEndTime, 1, 'last');
    crunchChunkLFP = rmsv(crunchChunkEphysStartInd:crunchChunkEphysEndInd);
    smoothedCrunchChunkLFP = smooth(crunchChunkLFP, 0.002);
    
    % locate the big crunch
    [crunchLFPMAX(i), maxInd] = max(crunchChunkLFP);
    crunchMaxInds(i) = maxInd + crunchChunkEphysStartInd;
    
    
    figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
    subplot(2, 1, 1);
    
    % process the mic signal to make spectrogram
    x = crunchChunkMic;
    x = detrend(x);
    x = x/std(x);
    
    winlen = 1024;
    win = blackman(winlen, 'periodic');
    hop = round(winlen/4);
    nfft = round(2*winlen);
    [~, F, T, STPS] = spectrogram(x, win, winlen-hop, nfft, audioFs, 'power');
    STPS = 10*log10(STPS);
    
    % plot the spectrogram
    surf(T, F, STPS);
    caxis([-20, 0])
    colormap hot
    shading interp
    axis tight
    box off
    view(0, 90)
    h = gca; h.XAxis.Visible = 'off';
    ylabel('Frequency, Hz')
    
    subplot(2, 1, 2)
    plot(ephysTime(crunchChunkEphysStartInd:crunchChunkEphysEndInd),crunchChunkLFP); hold on;
    plot(ephysTime(crunchMaxInds(i)), crunchLFPMAX(i), '.r', 'MarkerSize', 10);
    axis tight; box off;
    xlabel('time(sec)')

    %%%%%%%%%%%%%%%%%%%%% Processing chewing %%%%%%%%%%%%%%%%%%%%%%%%%%%
    trialStartTime = spike.foodTriggerTimes(i) + 1;
    chewingEndTime = trialStartTime + chewingSearchStartTimeShift + chewingSearchTimeWindow;
    
    % locate chewing period
    chewingChunkStartTime = trialStartTime + chewingSearchStartTimeShift;
    chewingChunkEndTime = chewingChunkStartTime + chewingTimeWindow;
    fitResults = [];
    
    while(chewingChunkEndTime <= chewingEndTime)
        chewingChunkVideoStartInd = find(videoTracking.frameTimestamps >= chewingChunkStartTime, 1, 'first');
        chewingChunkVideoEndInd = find(videoTracking.frameTimestamps <= chewingChunkEndTime, 1, 'last');
        
        chewingChunkJaw = videoTracking.jawDistance(chewingChunkVideoStartInd:chewingChunkVideoEndInd);
        fitResult = sineFit(1:length(chewingChunkJaw), chewingChunkJaw', 0);
        fitResults = [fitResults; fitResult];
        
        chewingChunkStartTime = chewingChunkStartTime + chewingSearchStepLength;
        chewingChunkEndTime = chewingChunkEndTime + chewingSearchStepLength;
    end
    
    fitMSE = fitResults(:, 5);
    tempInd = find(fitMSE == min(fitMSE(fitMSE>0)));
    minFitMSE = min(fitMSE(fitMSE>0));
    disp(['minFitMSE = ', num2str(minFitMSE)]);
    
    if fitResults(tempInd, 2) > 1.5 && fitResults(tempInd, 3) > 0.03
        disp(['trial ' num2str(i) ' selected'])
        chewingWindow(selectedTrialNum, 1) = trialStartTime + chewingSearchStartTimeShift + chewingSearchStepLength*(tempInd - 1);
        chewingWindow(selectedTrialNum, 2) = chewingWindow(selectedTrialNum, 1) + chewingTimeWindow;
        chewingChunkVideoStartInd = find(videoTracking.frameTimestamps >= chewingWindow(selectedTrialNum, 1), 1, 'first');
        chewingChunkVideoEndInd = find(videoTracking.frameTimestamps <= chewingWindow(selectedTrialNum, 2), 1, 'last');
        chewingChunkJaw = videoTracking.jawDistance(chewingChunkVideoStartInd:chewingChunkVideoEndInd);
        % sanity check
        sineFit(1:length(chewingChunkJaw), chewingChunkJaw'); % sanity check
        
        chewingJawDistance{selectedTrialNum, 1} = chewingChunkJaw;
        chewingJawTime{selectedTrialNum, 1} = videoTracking.frameTimestamps(chewingChunkVideoStartInd:chewingChunkVideoEndInd);
        chewingJawInds{selectedTrialNum, 1} = [chewingChunkVideoStartInd:chewingChunkVideoEndInd]';
        selectedTrialNum = selectedTrialNum + 1;
    end
   
end





