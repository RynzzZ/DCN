function analyzeTrial(session, varargin)

% This function breaks session data into each 'food' trial.
% Save data into 'trialAnalyzed.mat'. Each row is one 'food' trial.


% settings
s.hasLFP = true; % whether this session's recording contains a good LFP channel
s.hasMic = true; % whether this session contains good microphone recordings
s.crunchSearchTimeWindow = 4; % sec


if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs

% initializations
rootFolder = 'Z:\Qianyun\DCN\';
gitFolder = 'D:\DCN_Project\Github\DCN\';
sessionFolder = fullfile(rootFolder, 'Data', session);

% load spikeAnalyzed.mat
if ~exist(fullfile(sessionFolder, 'spikeAnalyzed.mat'), 'file')
    error('spikeAnalyzed.mat NOT exist! Run analyzeSession.m first!');
else
    load(fullfile(sessionFolder, 'spikeAnalyzed.mat'));
end

% get LFP channel data
ephysInfo = readtable(fullfile(rootFolder, 'Spreadsheets', 'ephysInfo.xlsx'));
LFPChannel = ephysInfo.LFP(strcmp(ephysInfo.session, session));
if isnan(LFPChannel) || LFPChannel == 0
    warning('No LFP channel in this session!');
    s.hasLFP = false;
else
    load(fullfile(rootFolder, 'Data', session, 'sessionEphysInfo.mat'));
    mapFile = sessionEphysInfo.mapFile;
    load(fullfile('Z:\obstacleData\ephys\channelMaps\kilosort', [mapFile, '.mat']), 'channelNum_OpenEphys');
    
    % load data
    getVoltage = @(data) double(data)*sessionEphysInfo.bitVolts; % extract voltage from memmapfile, converting to votlage, highpassing, and only return specific channel
    contFiles = dir(fullfile(sessionEphysInfo.ephysFolder, '*.continuous'));
    data = memmapfile(fullfile(sessionEphysInfo.ephysFolder, [contFiles(1).name(1:end-12), 's.dat']), ...
        'Format', {'int16', [sessionEphysInfo.channelNum sessionEphysInfo.smps], 'Data'}, 'Writable', false);
    LFPVoltage = getVoltage(data.Data.Data(channelNum_OpenEphys(LFPChannel), :));
end

% format data into trial structure
if s.hasLFP && s.hasMic % only proceed if this session has both good quality LFP and microphone recordngs
    trialTotalCount = spike.totalFoodNum;
    trial = table();
    for i = 1:trialTotalCount
        
        fprintf('trial %d/%d \n', i, trialTotalCount);
        
        % locate big crunch 
        trial.foodTriggerTime(i) = spike.foodTriggerTimes(i);
        crunchStartTime = spike.foodTriggerTimes(i) - 0.1;
        crunchEndTime = crunchStartTime + s.crunchSearchTimeWindow;
        
        crunchChunkMicStartInd = find(spike.micSignalTimes >= crunchStartTime, 1, 'first');
        crunchChunkMicEndInd = find(spike.micSignalTimes <= crunchEndTime, 1, 'last');
        crunchChunkMic = spike.micSignal(crunchChunkMicStartInd : crunchChunkMicEndInd);
        crunchChunkMic = highpass(crunchChunkMic, 100, 150000);
        
        crunchMicThreshold = max(crunchChunkMic)*0.8;
        bigCrunchTime = spike.micSignalTimes(find(crunchChunkMic >= crunchMicThreshold, 1, 'first') + crunchChunkMicStartInd);
        trial.CrunchTime(i) = bigCrunchTime;
        
        % get LFP around the big crunch
        crunchChunkEphysStartInd = find(sessionEphysInfo.convertedEphysTimestamps >= crunchStartTime, 1, 'first');
        crunchChunkEphysEndInd = find(sessionEphysInfo.convertedEphysTimestamps <= crunchEndTime, 1, 'last');
        crunchChunkLFP = LFPVoltage(crunchChunkEphysStartInd:crunchChunkEphysEndInd);
        
        crunchEphysStartTime = bigCrunchTime;
        crunchEphysEndTime = bigCrunchTime + 0.01;
        crunchEphysStartInd = find(sessionEphysInfo.convertedEphysTimestamps >= crunchEphysStartTime, 1, 'first');
        crunchEphysEndInd = find(sessionEphysInfo.convertedEphysTimestamps <= crunchEphysEndTime, 1, 'last');
        
        crunchLFP = LFPVoltage(crunchEphysStartInd:crunchEphysEndInd);
        crunchLFP(abs(crunchLFP) < 100) = 0;
        
        % calculate LFP peak to peak value for the big crunch
        [maxpks, maxlocs] = findpeaks(crunchLFP);
        [minpks, minlocs] = findpeaks(-crunchLFP);
        
        maxPeak = max(maxpks);
        minPeak = -max(minpks);
        crunchLFP_Peak2Peak = maxPeak - minPeak; % microvolt
        maxPeakLoc = maxlocs(maxpks == max(maxpks)) + crunchEphysStartInd;
        minPeakLoc = minlocs(minpks == max(minpks)) + crunchEphysStartInd;
        
        trial.crunchLFP_Peak2Peak(i) = crunchLFP_Peak2Peak;
        trial.crunchLFP_MaxPeak(i) = maxPeak;
        trial.crunchLFP_MinPeak(i) = minPeak;
        trial.crunchLFP_MaxPeakLoc(i) = maxPeakLoc;
        trial.crunchLFP_MaxPeakTime(i) = sessionEphysInfo.convertedEphysTimestamps(maxPeakLoc);
        trial.crunchLFP_MinPeakLoc(i) = minPeakLoc;
        trial.crunchLFP_MinPeakTime(i) = sessionEphysInfo.convertedEphysTimestamps(minPeakLoc);
        
        % plot big crunch, mic signal + LFP, with LFP Peak to Peak points
        if mod(i, 4) == 1
            figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
            plotInd = 1; 
            plotMatrix = [1, 3; 2, 4; 5, 7; 6, 8];
        end
        
        % mic signal
        subplot(4, 2, plotMatrix(plotInd, 1));
        plot(spike.micSignalTimes(crunchChunkMicStartInd:crunchChunkMicEndInd), crunchChunkMic); hold on
        plot([trial.foodTriggerTime(i), trial.foodTriggerTime(i)], [min(crunchChunkMic), max(crunchChunkMic)], 'k-');
        ylabel('volt');
        xlabel('time (s)');
        title(['session ', session, ' trial ', num2str(i)]);
        legend('mic signal', 'food trigger');
        axis tight;

        % LFP voltage
        subplot(4, 2, plotMatrix(plotInd, 2));
        plot(sessionEphysInfo.convertedEphysTimestamps(crunchChunkEphysStartInd:crunchChunkEphysEndInd),...
            LFPVoltage(crunchChunkEphysStartInd:crunchChunkEphysEndInd));
        hold on; axis tight;
        plot([trial.foodTriggerTime(i), trial.foodTriggerTime(i)], [min(crunchChunkLFP), max(crunchChunkLFP)], 'k-');
        plot(sessionEphysInfo.convertedEphysTimestamps(maxPeakLoc), maxPeak, 'r.', 'MarkerSize', 13);
        plot(sessionEphysInfo.convertedEphysTimestamps(minPeakLoc), minPeak, 'y.', 'MarkerSize', 13);
        ylabel('microVolt');
        xlabel('time (s)');
        title('LFP Signal');
        axis tight;
        plotInd = plotInd + 1;
    end   
    save(fullfile(sessionFolder, 'trialAnalyzed.mat'), 'trial');
end



end