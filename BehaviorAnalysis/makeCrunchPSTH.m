function makeCrunchPSTH(session, ephysChannel, varargin)

% settings
s.supressFigure = false; % whether to only process the data but supress plotting figures
s.crunchSearchTimeWindow = 4; % sec
s.crunchTimeWindow = [-0.005, 0.04]; % sec
s.ephysThresh = -300;
s.crunchRefPeriod = 0.02;
s.micThresh = 0.03;
s.fs = 150000; % sampling rate for mic signal
s.ephysFs = 30000; % sampling rate for ephys
s.spkRateFs = 1000; % sampling rate for unit spkRate

if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs

% initializations
rootFolder = 'Z:\Qianyun\DCN\';
gitFolder = 'D:\DCN_Project\Github\DCN\';
sessionFolder = fullfile(rootFolder, 'Data', session);
saveFolder = fullfile(rootFolder, 'Data', session, 'trialFigs', 'crunch');
figureFolder = 'Z:\Qianyun\DCN\Figures\PSTH\crunch\';


% load spikeAnalyzed.mat
if ~exist(fullfile(sessionFolder, 'spikeAnalyzed.mat'), 'file')
    error('spikeAnalyzed.mat NOT exist! Run analyzeSession.m first!');
else
    disp('loading spikeAnalyzed.mat...')
    load(fullfile(sessionFolder, 'spikeAnalyzed.mat'), 'spike');
end

% load cameraAnalyzed.mat
if ~exist(fullfile(sessionFolder, 'cameraAnalyzed.mat'), 'file')
    error('cameraAnalyzed.mat NOT exist! Run analyzeSession.m first!');
else
    disp('loading cameraAnalyzed.mat...')
    load(fullfile(sessionFolder, 'cameraAnalyzed.mat'), 'videoTracking');
    if ~any(strcmp('jawDistance', videoTracking.Properties.VariableNames))
        warning('JawDistance does NOT exist in videoTracking! Will NOT analyze chewing LFP!');
    end
end


% get good unit channel data
load(fullfile(rootFolder, 'Data', session, 'sessionEphysInfo.mat'), 'sessionEphysInfo');
ephysTime = sessionEphysInfo.convertedEphysTimestamps;
mapFile = sessionEphysInfo.mapFile;
load(fullfile('Z:\obstacleData\ephys\channelMaps\kilosort', [mapFile, '.mat']), 'channelNum_OpenEphys');

getVoltage = @(data) double(data)*sessionEphysInfo.bitVolts; % extract voltage from memmapfile, converting to votlage, highpassing, and only return specific channel
contFiles = dir(fullfile(sessionEphysInfo.ephysFolder, '*.continuous'));
data = memmapfile(fullfile(sessionEphysInfo.ephysFolder, [contFiles(1).name(1:end-12), 's.dat']), ...
    'Format', {'int16', [sessionEphysInfo.channelNum sessionEphysInfo.smps], 'Data'}, 'Writable', false);
% unitVoltage = getVoltage(data.Data.Data(channelNum_OpenEphys(bestChannel), :));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%processing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trialTotalCount = spike.totalFoodNum;
crunchTimes = [];
crunchNum = 1;
spkRate = [];
for i = 1:trialTotalCount
    
    fprintf('trial %d/%d \n', i, trialTotalCount);
    
    %%%%%%%%%%%%%%%%%%%%% Processing crunch %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get big crunch chunk
    crunchStartTime = spike.foodTriggerTimes(i);
    crunchEndTime = crunchStartTime + s.crunchSearchTimeWindow;
    
    crunchChunkBins = spike.micSignalTimes >= crunchStartTime & spike.micSignalTimes <= crunchEndTime;
    crunchChunkTime = spike.micSignalTimes(crunchChunkBins);
    crunchChunkMic = spike.micSignal(crunchChunkBins);
    crunchChunkMic = highpass(crunchChunkMic, 100, s.fs);
    crunchChunkMicRMSV = sqrt(movmean(crunchChunkMic.^2, 100));
    
    
    % locate crunch within crunch chunk
    [~, locs, n] = crossdet(crunchChunkMicRMSV, s.micThresh, 'thresholdUp');
    inds = find(diff(locs)<s.crunchRefPeriod*s.ephysFs);
    locs(inds+1) = [];
    n = n - length(inds);
    
    figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
    tiledlayout('flow')
    nexttile
    plot(crunchChunkTime, crunchChunkMicRMSV); hold on;
    plot(crunchChunkTime(locs), repmat(s.micThresh, size(crunchChunkTime(locs))), 'r.', 'MarkerSize', 10);

    for j = 1:n
        crunchTime = crunchChunkTime(locs(j));
        crunchTimes = [crunchTimes; crunchChunkTime(locs(j))];
        crunchStartTime = crunchTime + s.crunchTimeWindow(1);
        crunchEndTime = crunchTime + s.crunchTimeWindow(2);
        crunchBins = spike.micSignalTimes >= crunchStartTime & spike.micSignalTimes <= crunchEndTime;
        crunchMic(crunchNum, :) = spike.micSignal(crunchBins);
        
        % get crunch ephys data
        crunchEphysStartInd = find(ephysTime >= crunchStartTime, 1, 'first');
        crunchEphysStopInd = find(ephysTime <= crunchEndTime, 1, 'last');
        % crunchTimeWindow = ephysTime(crunchEphysStartInd:crunchEphysStopInd);
        clear crunchEphysData
        crunchEphysData = getVoltage(data.Data.Data(channelNum_OpenEphys(ephysChannel), crunchEphysStartInd:crunchEphysStopInd));
        spkTimeWindow = linspace(s.crunchTimeWindow(1), s.crunchTimeWindow(2), length(crunchEphysData));
        
        nexttile
        plot(spkTimeWindow, crunchEphysData);
        axis tight; hold on;
        
        % find spikes (by threshold crossing)
        [~, spikeInds, ~] = crossdet(crunchEphysData, s.ephysThresh, 'thresholdDown');
        spikeTimes = spkTimeWindow(spikeInds);
        [spkRate(crunchNum, :), spikeRateTimes] = getFiringRate(spikeTimes, 'tLims', [s.crunchTimeWindow(1), s.crunchTimeWindow(2)]);
        plot(spikeRateTimes, spkRate(crunchNum, :), 'r-'); 
        crunchNum = crunchNum + 1;
    end
     % save figs
     saveas(gcf, fullfile(saveFolder, ['trial_', num2str(i), '_crunchMicEphys.png']))
end


figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
yyaxis left
shadedErrorBar(linspace(s.crunchTimeWindow(1)*1000, s.crunchTimeWindow(2)*1000, size(crunchMic, 2)),...
    mean(crunchMic), std(crunchMic), 'lineProps', {'-', 'Color', '#0072BD', 'LineWidth', 2});
ylabel('mic rmsv');
ylim([-0.3, 0.3]);

yyaxis right
shadedErrorBar(linspace(s.crunchTimeWindow(1)*1000, s.crunchTimeWindow(2)*1000, size(spkRate, 2)),...
    mean(spkRate), std(spkRate), 'lineProps', {'-', 'Color', '#A2142F', 'LineWidth', 2});
ylim([-160, 160]);
xlabel('time (ms)');
ylabel('spk rate (spk/s)');


% save figs
saveas(gcf, fullfile(saveFolder, ['crunchPSTH.png']));
saveas(gcf, fullfile(figureFolder, [session, '_CH', num2str(ephysChannel), '_crunchPSTH.png']));



end