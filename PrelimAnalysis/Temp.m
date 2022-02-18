%% temp script 

session = '20211108_000';
trialID = 2;

% settings
s.hasMic = true; % whether this session contains good microphone recordings
s.hasJawTrace = true; % whether this session contains video tracking for jaw trace
s.supressFigure = false; % whether to only process the data but supress plotting figures
s.unitID = 128;
s.crunchSearchTimeWindow = 4; % sec
s.crunchTimeWindow = 0.02; % sec
s.chewingSearchTimeWindow = 10; % sec
s.chewingTimeWindow = 3; % sec
s.chewingSearchStepLength = 0.1; % sec
s.fs = 150000; % sampling rate for mic signal
s.ephysFs = 30000; % sampling rate for ephys
s.spkRateFs = 1000; % sampling rate for unit spkRate

if exist('varargin', 'var'); for trialID = 1:2:length(varargin); s.(varargin{trialID}) = varargin{trialID+1}; end; end  % parse name-value pairs

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
end

% load cameraAnalyzed.mat
if ~exist(fullfile(sessionFolder, 'cameraAnalyzed.mat'), 'file')
    error('cameraAnalyzed.mat NOT exist! Run analyzeSession.m first!');
else
    disp('loading cameraAnalyzed.mat...')
    load(fullfile(sessionFolder, 'cameraAnalyzed.mat'), 'videoTracking');
    if ~any(strcmp('jawDistance', videoTracking.Properties.VariableNames))
        warning('JawDistance does NOT exist in videoTracking! Will NOT analyze chewing LFP!');
        s.hasJawTrace = false;
    end
end


% get good unit channel data
load(fullfile(rootFolder, 'Data', session, 'sessionEphysInfo.mat'), 'sessionEphysInfo');
mapFile = sessionEphysInfo.mapFile;
load(fullfile('Z:\obstacleData\ephys\channelMaps\kilosort', [mapFile, '.mat']), 'channelNum_OpenEphys');

neuralData = load(fullfile(sessionFolder, 'neuralData.mat'));
ind = find(neuralData.unit_ids == s.unitID);
spkRates = neuralData.spkRates(ind, :);
spkRateTimes = neuralData.timeStamps;
bestChannel = neuralData.bestChannels(ind, :);

getVoltage = @(data) double(data)*sessionEphysInfo.bitVolts; % extract voltage from memmapfile, converting to votlage, highpassing, and only return specific channel
contFiles = dir(fullfile(sessionEphysInfo.ephysFolder, '*.continuous'));
data = memmapfile(fullfile(sessionEphysInfo.ephysFolder, [contFiles(1).name(1:end-12), 's.dat']), ...
    'Format', {'int16', [sessionEphysInfo.channelNum sessionEphysInfo.smps], 'Data'}, 'Writable', false);
unitVoltage = getVoltage(data.Data.Data(channelNum_OpenEphys(bestChannel), :));


% load data
ephysInfo = readtable(fullfile(rootFolder, 'Spreadsheets', 'ephysInfo.xlsx'));
LFPChannel = ephysInfo.LFP(strcmp(ephysInfo.session, session));

LFPVoltage = getVoltage(data.Data.Data(channelNum_OpenEphys(LFPChannel), :));
rmsv = sqrt(movmean(LFPVoltage.^2, 100));


fprintf('trial %d \n', trialID);
%%%%%%%%%%%%%%%%%%%%% Processing crunch %%%%%%%%%%%%%%%%%%%%%%%%%%%
% get big crunch chunk
crunchStartTime = spike.foodTriggerTimes(trialID) - 0.1;
crunchEndTime = crunchStartTime + s.crunchSearchTimeWindow;

% get Mic crunch chunk data
crunchChunkMicStartInd = find(spike.micSignalTimes >= crunchStartTime, 1, 'first');
crunchChunkMicEndInd = find(spike.micSignalTimes <= crunchEndTime, 1, 'last');
crunchChunkMic = spike.micSignal(crunchChunkMicStartInd : crunchChunkMicEndInd);
crunchChunkMic = highpass(crunchChunkMic, 100, 150000);
crunchChunkMicRMSV = sqrt(movmean(crunchChunkMic.^2, 100));

% get Unit crunch chunk data
crunchChunkSpkRateStartInd = find(spkRateTimes >= crunchStartTime, 1, 'first');
crunchChunkSpkRateEndInd = find(spkRateTimes <= crunchEndTime, 1, 'last');
crunchChunkSpkRateTimes = spkRateTimes(crunchChunkSpkRateStartInd:crunchChunkSpkRateEndInd);
crunchChunkSpkRate = spkRates(crunchChunkSpkRateStartInd:crunchChunkSpkRateEndInd); % unit spkRate

% get LFP crunch chunk data (span = 0.002)
crunchChunkLFPStartInd = find(sessionEphysInfo.convertedEphysTimestamps >= crunchStartTime, 1, 'first');
crunchChunkLFPEndInd = find(sessionEphysInfo.convertedEphysTimestamps <= crunchEndTime, 1, 'last');
crunchChunkLFPTimes = sessionEphysInfo.convertedEphysTimestamps(crunchChunkLFPStartInd:crunchChunkLFPEndInd);
crunchChunkLFP = rmsv(crunchChunkLFPStartInd:crunchChunkLFPEndInd);
smoothedCrunchChunkLFP = smooth(crunchChunkLFP, 0.002);

%%%%%%%%%%%%%%%%%%%%% Processing chewing %%%%%%%%%%%%%%%%%%%%%%%%%%%
trialStartTime = spike.foodTriggerTimes(trialID) - 0.1;
chewingEndTime = trialStartTime + s.chewingSearchTimeWindow;

% locate chewing period
chewingChunkStartTime = trialStartTime + 4;
chewingChunkEndTime = chewingChunkStartTime + s.chewingTimeWindow;
fitMSE = [];
while(chewingChunkEndTime <= chewingEndTime)
    chewingChunkVideoStartInd = find(videoTracking.frameTimestamps >= chewingChunkStartTime, 1, 'first');
    chewingChunkVideoEndInd = find(videoTracking.frameTimestamps <= chewingChunkEndTime, 1, 'last');
    
    chewingChunkJaw = videoTracking.jawDistance(chewingChunkVideoStartInd:chewingChunkVideoEndInd);
    fitResult = sineFit(1:length(chewingChunkJaw), chewingChunkJaw', 0);
    fitMSE = [fitMSE, fitResult(end)];
    
    chewingChunkStartTime = chewingChunkStartTime + s.chewingSearchStepLength;
    chewingChunkEndTime = chewingChunkEndTime + s.chewingSearchStepLength;
end

tempInd = find(fitMSE == min(fitMSE(fitMSE>0)));
disp(['minFitMSE = ', num2str(min(fitMSE(fitMSE>0)))]);
chewingChunkStartTime = trialStartTime + 4 + s.chewingSearchStepLength*(tempInd - 1);
chewingChunkEndTime = chewingChunkStartTime + s.chewingTimeWindow;
chewingChunkVideoStartInd = find(videoTracking.frameTimestamps >= chewingChunkStartTime, 1, 'first');
chewingChunkVideoEndInd = find(videoTracking.frameTimestamps <= chewingChunkEndTime, 1, 'last');

chewingChunkJaw = videoTracking.jawDistance(chewingChunkVideoStartInd:chewingChunkVideoEndInd);
sineFit(1:length(chewingChunkJaw), chewingChunkJaw'); % sanity check

% get Unit chewing peirod data
chewingChunkSpkRateStartInd = find(spkRateTimes >= chewingChunkStartTime, 1, 'first');
chewingChunkSpkRateEndInd = find(spkRateTimes <= chewingChunkEndTime, 1, 'last');
chewingChunkSpkRateTimes = spkRateTimes(chewingChunkSpkRateStartInd:chewingChunkSpkRateEndInd);
chewingChunkSpkRate = spkRates(chewingChunkSpkRateStartInd:chewingChunkSpkRateEndInd); % unit spkRate

% get LFP chewing period data (span = 0.009)
chewingChunkLFPStartInd = find(sessionEphysInfo.convertedEphysTimestamps >= chewingChunkStartTime, 1, 'first');
chewingChunkLFPEndInd = find(sessionEphysInfo.convertedEphysTimestamps <= chewingChunkEndTime, 1, 'last');
chewingChunkLFPTimes = sessionEphysInfo.convertedEphysTimestamps(chewingChunkLFPStartInd:chewingChunkLFPEndInd);
chewingChunkLFP = rmsv(chewingChunkLFPStartInd:chewingChunkLFPEndInd);
smoothedChewingChunkLFP = smooth(chewingChunkLFP, 0.009);

%% process chunk data (crunch + chewing)
chunkLength = 8.86; %sec

% get chunk start and stop time
chunkStartTime = spike.foodTriggerTimes(trialID);
chunkEndTime = chunkStartTime + chunkLength;

% get Mic chunk data
chunkMicStartInd = find(spike.micSignalTimes >= chunkStartTime, 1, 'first');
chunkMicEndInd = find(spike.micSignalTimes <= chunkEndTime, 1, 'last');
chunkMic = spike.micSignal(chunkMicStartInd : chunkMicEndInd);
chunkMic = highpass(chunkMic, 100, 150000);
chunkMicRMSV = sqrt(movmean(chunkMic.^2, 100));

% get LFP chunk data (span = 0.002)
chunkLFPStartInd = find(sessionEphysInfo.convertedEphysTimestamps >= chunkStartTime, 1, 'first');
chunkLFPEndInd = find(sessionEphysInfo.convertedEphysTimestamps <= chunkEndTime, 1, 'last');
chunkLFPTimes = sessionEphysInfo.convertedEphysTimestamps(chunkLFPStartInd:chunkLFPEndInd);
chunkLFP = LFPVoltage(chunkLFPStartInd:chunkLFPEndInd);
chunkLFPrmsv = rmsv(chunkLFPStartInd:chunkLFPEndInd);
smoothedChunkLFPrmsv = smooth(chunkLFPrmsv, 0.002);

% get tongue and jaw data
chunkVideoStartInd = find(videoTracking.frameTimestamps >= chunkStartTime, 1, 'first');
chunkVideoEndInd = find(videoTracking.frameTimestamps <= chunkEndTime, 1, 'last');
jawTrace = videoTracking.jawDistance(chunkVideoStartInd:chunkVideoEndInd);
tongueConfidence = videoTracking.tongue_confidence(chunkVideoStartInd:chunkVideoEndInd);
tongueExist = videoTracking.tongueExist(chunkVideoStartInd:chunkVideoEndInd);

%% plot
figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
colorMatrix = {'#EDB120', 	'#77AC30', 	'#7E2F8E', 	'#0072BD'};
fontSize = 14;

% plot tongue traces
x = linspace(0, chunkLength, length(jawTrace));
y = (tongueConfidence/max(tongueConfidence)) + 3.3;
plot(x, y, '-', 'Color', colorMatrix{1}, 'LineWidth', 1.2); hold on; box off;

% plot jaw traces
y = (jawTrace/max(jawTrace)) + 2;
plot(x, y, '-', 'Color', colorMatrix{2}, 'LineWidth', 1.2);

% plot LFP
% x = linspace(0, chunkLength, length(chunkLFP));
% y = chunkLFP/(max(chunkLFP) - min(chunkLFP)) + 1.5;
% plot(x, y, '-', 'Color', colorMatrix{3}, 'LineWidth', 1.2);

% plot LFP
x = linspace(0, chunkLength, length(smoothedChunkLFPrmsv));
y = smoothedChunkLFPrmsv/max(smoothedChunkLFPrmsv) + 0.8;
plot(x, y, '-', 'Color', colorMatrix{3}, 'LineWidth', 1.2);

% mic rmsv
x = linspace(0, chunkLength, length(chunkMicRMSV));
y = chunkMicRMSV/max(chunkMicRMSV)-0.2;
plot(x, y, '-', 'Color', colorMatrix{4}, 'LineWidth', 1.2);

h = gca; h.XAxis.Visible = 'off'; 
h.YAxis.Visible = 'off';


% pellet dispensed
plot([0, 0], [0, 4.5], 'k--');
text(-0.5, 4.6, 'pellet dispensed', 'FontWeight', 'bold', 'FontSize', fontSize);

% pellet arrival
pelletArrivalTime = videoTracking.frameTimestamps(5044);
pelletArrivalTime = pelletArrivalTime - chunkStartTime;
plot([pelletArrivalTime, pelletArrivalTime], [0, 4.5], 'k--'); 
text(pelletArrivalTime-0.35, 4.6, 'pellet arrives', 'FontWeight', 'bold', 'FontSize', fontSize);

% crunch line
crunchTime = 57.36 - chunkStartTime;
plot([crunchTime, crunchTime], [2, 4.5], 'k--' );
text(crunchTime-0.17, 4.6, 'crunch', 'FontWeight', 'bold', 'FontSize', fontSize);

% ryhthmic chewing line
plot([5, chunkLength], [4.5, 4.5], 'k--');
text(6.5, 4.6, 'ryhthmic chewing', 'FontWeight', 'bold', 'FontSize', fontSize);

% text
text(-0.95, 3.33, 'tongue reach', 'Color', colorMatrix{1}, 'FontWeight', 'bold', 'FontSize', fontSize);
text(-0.9, 2.3, 'jaw distance', 'Color', colorMatrix{2}, 'FontWeight', 'bold', 'FontSize', fontSize);
text(-0.35, 1.13, 'LFP', 'Color', colorMatrix{3}, 'FontWeight', 'bold', 'FontSize', fontSize);
text(-0.33, 0.07, 'mic', 'Color', colorMatrix{4}, 'FontWeight', 'bold', 'FontSize', fontSize)

box off; axis tight;
ylim([-0.1, 4.6])




%% Plot!!

figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
rows = 6; cols = 2;

% series of plots for big crunch
plotMatrix = [1, 3, 5, 7, 9, 11];
colorMatrix = {' ', '#0072BD', 	'#EDB120', 	'#77AC30', 	'#7E2F8E', 	'#A2142F'};

% plot mic signal spectrogram
plotInd = 1;
subplot(rows, cols, plotMatrix(plotInd));

% process the mic signal to make spectrogram
x = crunchChunkMic;
x = detrend(x);
x = x/std(x);

winlen = 1024;
win = blackman(winlen, 'periodic');
hop = round(winlen/4);
nfft = round(2*winlen);
[~, F, T, STPS] = spectrogram(x, win, winlen-hop, nfft, s.fs, 'power');
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
title([session(1:8), '_', session(10:end) ' Trial ', num2str(trialID), ' first big crunch'], 'Interpreter','none');

% plot mic signal RMSV
plotInd = plotInd + 1;
subplot(rows, cols, plotMatrix(plotInd));

plot(spike.audioSignalTimes(crunchChunkMicStartInd:crunchChunkMicEndInd), crunchChunkMicRMSV, ...
    '-', 'Color', colorMatrix{plotInd});
hold on; box off; axis tight;
plot([spike.foodTriggerTimes(trialID), spike.foodTriggerTimes(trialID)], [min(crunchChunkMicRMSV), max(crunchChunkMicRMSV)], 'k-');
ylabel('RMS Value');
legend('Mic RMSV', 'Food Trigger', 'Location', 'northwest'); legend boxoff;
h = gca; h.XAxis.Visible = 'off';


% plot tongue and jaw traces
crunchChunkVideoStartInd = find(videoTracking.frameTimestamps >= crunchStartTime, 1, 'first');
crunchChunkVideoEndInd = find(videoTracking.frameTimestamps <= crunchEndTime, 1, 'last');
jawTrace = videoTracking.jawDistance(crunchChunkVideoStartInd:crunchChunkVideoEndInd);

plotInd = plotInd + 1;
subplot(rows, cols, plotMatrix(plotInd));
plot(videoTracking.frameTimestamps(crunchChunkVideoStartInd:crunchChunkVideoEndInd), ...
    videoTracking.tongue_confidence(crunchChunkVideoStartInd:crunchChunkVideoEndInd), ...
    '-', 'Color', colorMatrix{plotInd});
hold on; box off; axis tight;
plot([spike.foodTriggerTimes(trialID), spike.foodTriggerTimes(trialID)], [0, 1], 'k-');
ylabel('confidence');
ylim([0, 1]);
legend('Tongue Existence','Location', 'northwest'); legend boxoff;
h = gca; h.XAxis.Visible = 'off';

plotInd = plotInd + 1;
subplot(rows, cols, plotMatrix(plotInd));
plot(videoTracking.frameTimestamps(crunchChunkVideoStartInd:crunchChunkVideoEndInd), jawTrace,...
    '-', 'Color', colorMatrix{plotInd});
hold on; box off; axis tight;
plot([spike.foodTriggerTimes(trialID), spike.foodTriggerTimes(trialID)], [min(jawTrace), max(jawTrace)], 'k-');
ylabel('pixel');
legend('Jaw Distance', 'Location', 'northwest'); legend boxoff;
h = gca; h.XAxis.Visible = 'off';

% crunch chunk LFP (span = 0.002)
plotInd = plotInd + 1;
subplot(rows, cols, plotMatrix(plotInd));

plot(crunchChunkLFPTimes, smoothedCrunchChunkLFP, '-', 'Color', colorMatrix{plotInd});
hold on; box off; axis tight;
plot([spike.foodTriggerTimes(trialID), spike.foodTriggerTimes(trialID)], [min(smoothedCrunchChunkLFP), max(smoothedCrunchChunkLFP)], 'k-');
ylabel('rms value');
legend('LFP RMSV, span = 0.002', 'Location', 'northwest'); legend boxoff;
h = gca; h.XAxis.Visible = 'off';


% Unit spkRate
plotInd = plotInd + 1;
subplot(rows, cols, plotMatrix(plotInd));

plot(crunchChunkSpkRateTimes, crunchChunkSpkRate, '-', 'Color', colorMatrix{plotInd});
hold on; box off; axis tight;
plot([spike.foodTriggerTimes(trialID), spike.foodTriggerTimes(trialID)], [min(crunchChunkSpkRate), max(crunchChunkSpkRate)], 'k-');
ylabel('spikes/s');
legend('spike rate', 'Location', 'northwest'); legend boxoff;
xlabel('time (s)');

%%%%%%%%%%% series of plots for chewing period %%%%%%%%%%%
plotMatrix = [2, 4, 6, 8, 10, 12];
colorMatrix = {' ', '#0072BD', 	'#EDB120', 	'#77AC30', '#7E2F8E', '#A2142F'};

% plot mic signal spectrogram
plotInd = 1;
subplot(rows, cols, plotMatrix(plotInd));

% process the mic signal to make spectrogram
chewingChunkMicStartInd = find(spike.micSignalTimes >= chewingChunkStartTime, 1, 'first');
chewingChunkMicEndInd = find(spike.micSignalTimes <= chewingChunkEndTime, 1, 'last');
chewingChunkMic = spike.micSignal(chewingChunkMicStartInd : chewingChunkMicEndInd);
chewingChunkMic = highpass(chewingChunkMic, 100, 150000);
chewingChunkMicRMSV = sqrt(movmean(chewingChunkMic.^2, 100));

% plot mic signal spectrogram
x = chewingChunkMic;
x = detrend(x);
x = x/std(x);

winlen = 1024;
win = blackman(winlen, 'periodic');
hop = round(winlen/4);
nfft = round(2*winlen);
[~, F, T, STPS] = spectrogram(x, win, winlen-hop, nfft, s.fs, 'power');
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
title([session(1:8), '_', session(10:end) ' Trial ', num2str(trialID), ' rhythmic chewing'], 'Interpreter','none');

% plot mic signal RMSV
plotInd = plotInd + 1;
subplot(rows, cols, plotMatrix(plotInd));

plot(spike.audioSignalTimes(chewingChunkMicStartInd:chewingChunkMicEndInd), chewingChunkMicRMSV, ...
    '-', 'Color', colorMatrix{plotInd});
box off; axis tight;
ylabel('RMS Value');
legend('Mic RMSV', 'Location', 'northwest'); legend boxoff;
h = gca; h.XAxis.Visible = 'off';


% plot tongue traces
plotInd = plotInd + 1;
subplot(rows, cols, plotMatrix(plotInd));
plot(videoTracking.frameTimestamps(chewingChunkVideoStartInd:chewingChunkVideoEndInd), ...
    videoTracking.tongue_confidence(chewingChunkVideoStartInd:chewingChunkVideoEndInd), ...
    '-', 'Color', colorMatrix{plotInd});
hold on; box off; axis tight;
ylabel('confidence');
ylim([0, 1]);
legend('Tongue Existence','Location', 'northwest'); legend boxoff;
h = gca; h.XAxis.Visible = 'off';

% plot jaw traces
plotInd = plotInd + 1;
subplot(rows, cols, plotMatrix(plotInd));
plot(videoTracking.frameTimestamps(chewingChunkVideoStartInd:chewingChunkVideoEndInd), ...
    chewingChunkJaw, '-', 'Color', colorMatrix{plotInd});
box off; axis tight;
ylabel('pixel');
legend('Jaw Distance', 'Location', 'northwest'); legend boxoff;
h = gca; h.XAxis.Visible = 'off';

% LFP rmsv (span = 0.009)
plotInd = plotInd + 1;
subplot(rows, cols, plotMatrix(plotInd));

plot(chewingChunkLFPTimes, smoothedChewingChunkLFP, '-', 'Color', colorMatrix{plotInd});
hold on; box off; axis tight;
ylabel('microvolt');
legend('LFP RMSV, span = 0.009', 'Location', 'northwest'); legend boxoff;
h = gca; h.XAxis.Visible = 'off';

% Unit spkRates
plotInd = plotInd + 1;
subplot(rows, cols, plotMatrix(plotInd));

plot(chewingChunkSpkRateTimes, chewingChunkSpkRate, '-', 'Color', colorMatrix{plotInd});
box off; axis tight;
ylabel('spikes/s');
legend('spike rate', 'Location', 'northwest'); legend boxoff;
xlabel('time (s)');




%% Cross Correlation!

% calculate the cross corr b/w jaw distance and unit spkRate
maxLag = 0.2*s.spkRateFs;
smoothedChewingChunkSpkRate = smooth(chewingChunkSpkRate, 0.03);
cameraTimestamps = videoTracking.frameTimestamps(chewingChunkVideoStartInd:chewingChunkVideoEndInd);
chewingChunkJaw_interp = interp1(cameraTimestamps, chewingChunkJaw, chewingChunkSpkRateTimes);
chewingChunkSpkRate_xcorr = smoothedChewingChunkSpkRate(~isnan(chewingChunkJaw_interp));
chewingChunkJaw_xcorr = chewingChunkJaw_interp(~isnan(chewingChunkJaw_interp));

[c, lags] = xcorr(chewingChunkJaw_xcorr, chewingChunkSpkRate_xcorr, maxLag, 'coeff');
chewingCrossCorr_JawUnit = c;
[maxCorrJawUnit, maxCorrLocJawUnit] = max(c);

% calculate the cross corr b/w jaw distance and LFP
maxLag = 0.2*s.ephysFs;
cameraTimestamps = videoTracking.frameTimestamps(chewingChunkVideoStartInd:chewingChunkVideoEndInd);
chewingChunkJaw_interp = interp1(cameraTimestamps, chewingChunkJaw, chewingChunkLFPTimes);
chewingChunkLFP_xcorr = smoothedChewingChunkLFP(~isnan(chewingChunkJaw_interp));
chewingChunkJaw_xcorr = chewingChunkJaw_interp(~isnan(chewingChunkJaw_interp));

[c, lags] = xcorr(chewingChunkJaw_xcorr, chewingChunkLFP_xcorr, maxLag, 'coeff');
chewingCrossCorr_JawLFP = c;
[maxCorrJawLFP, maxCorrLocJawLFP] = max(c);


% calculate the cross corr b/w LFP and SpkRate
maxLag = 0.2*s.ephysFs;
chewingChunkSpkRate_interp = interp1(chewingChunkSpkRateTimes, smoothedChewingChunkSpkRate, chewingChunkLFPTimes);
chewingChunkLFP_xcorr = smoothedChewingChunkLFP(~isnan(chewingChunkSpkRate_interp));
chewingChunkSpkRate_xcorr = chewingChunkSpkRate_interp(~isnan(chewingChunkSpkRate_interp));

[c, lags] = xcorr(chewingChunkSpkRate_xcorr, chewingChunkLFP_xcorr, maxLag, 'coeff');
chewingCrossCorr_SpkRateLFP= c;
[maxCorrSpkRateLFP, maxCorrLocSpkRateLFP] = max(c);

%% plot cross corrs

% cross corr: Jaw + Unit SpkRate
figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
tiledlayout('flow');
nexttile

maxLag = 0.2*s.spkRateFs;
x = [-maxLag:maxLag]/s.spkRateFs;
y = chewingCrossCorr_JawUnit;
plot(x, y, '-', 'LineWidth', 2);

hold on; axis tight; box off;
xlabel('time (s)');
title(['trial ', num2str(trialID), ' cross corr b/w chewing jaw trace & Unit SpkRate'], 'Interpreter','none');

% plot the line for max cross corr and find its corresponding time point
maxCorrLoc = maxCorrLocJawUnit;
maxCorrLoc = (maxCorrLoc - maxLag)/s.spkRateFs;
yLimits = get(gca,'YLim');
plot([maxCorrLoc, maxCorrLoc], [min(yLimits), max(yLimits)], 'r-');
legend('chewing & unit spkRate cross corr', ['max corr loc = ', num2str(maxCorrLoc), ' s']);
legend boxoff;

% cross corr: Jaw + LFP
nexttile

maxLag = 0.2*s.ephysFs;
x = [-maxLag:maxLag]./s.ephysFs;
y = chewingCrossCorr_JawLFP;
plot(x, y, '-', 'LineWidth', 2);

hold on; axis tight; box off;
xlabel('time (s)');
title(['trial ', num2str(trialID), ' cross corr b/w chewing jaw trace & LFP'], 'Interpreter','none');

% plot the line for max cross corr and find its corresponding time point
maxCorrLoc = maxCorrLocJawLFP;
maxCorrLoc = (maxCorrLoc - maxLag)/s.ephysFs;
yLimits = get(gca,'YLim');
plot([maxCorrLoc, maxCorrLoc], [min(yLimits), max(yLimits)], 'r-');
legend('chewing & LFP cross corr', ['max corr loc = ', num2str(maxCorrLoc), ' s']);
legend boxoff;


% cross corr: SpkRate + LFP
nexttile

maxLag = 0.2*s.ephysFs;
x = [-maxLag:maxLag]./s.ephysFs;
y = chewingCrossCorr_SpkRateLFP;
plot(x, y, '-', 'LineWidth', 2);

hold on; axis tight; box off;
xlabel('time (s)');
title(['trial ', num2str(trialID), ' cross corr b/w chewing spkRate & LFP'], 'Interpreter','none');

% plot the line for max cross corr and find its corresponding time point
maxCorrLoc = maxCorrLocSpkRateLFP;
maxCorrLoc = (maxCorrLoc - maxLag)/s.ephysFs;
yLimits = get(gca,'YLim');
plot([maxCorrLoc, maxCorrLoc], [min(yLimits), max(yLimits)], 'r-');
legend('unit spkRate & LFP cross corr', ['max corr loc = ', num2str(maxCorrLoc), ' s']);
legend boxoff;




%% test

figure;
subplot(3, 1, 1); plot(smoothedChewingChunkSpkRate); axis tight; 

subplot(3, 1, 2); plot(smoothedChewingChunkLFP); axis tight; 

subplot(3, 1, 3); plot(chewingChunkJaw_xcorr); axis tight;



%% jaw phase & LFP

ephysChannel = 30;

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
    
    LFPVoltage = getVoltage(data.Data.Data(channelNum_OpenEphys(ephysChannel), :));
    rmsv = sqrt(movmean(LFPVoltage.^2, 100));
end
ephysTime = sessionEphysInfo.convertedEphysTimestamps;

ephysStartInds = nan(length(goodChewsStartStopTime), 1);
ephysStopInds = nan(length(goodChewsStartStopTime), 1);
LFPChunks = cell(length(goodChewsStartStopTime), 1);
for i = 1:length(goodChewsStartStopTime)
    ephysStartInds(i) = find(ephysTime > goodChewsStartStopTime(i, 1), 1, 'first');
    ephysStopInds(i) = find(ephysTime < goodChewsStartStopTime(i, 2), 1, 'last');
    
    % LFPChunks{i, 1} = getVoltage(data.Data.Data(channelNum_OpenEphys(ephysChannel), ephysStartInds(i, 1):ephysStopInds(i, 1))); 
    LFPChunks{i, 1} = rmsv(ephysStartInds(i):ephysStopInds(i));
end

figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
rows = 2;
cols = 5;
plotInd = 1;

selectedChews = [];
for i = 1:length(LFPChunks)
    subplot(rows, cols, plotInd);
    
    x1 = linspace(goodChewsStartStopTime(i, 1), goodChewsStartStopTime(i, 2), length(goodChewsJawDistance{i, 1}));
    y1 = goodChewsJawDistance{i, 1};
    % y1 = lowpass(y1, 5, 1/(x1(2) - x1(1)));
    
    x2 = linspace(ephysTime(ephysStartInds(i)), ephysTime(ephysStopInds(i)), length(LFPChunks{i, 1}));
    y2 = LFPChunks{i, 1};
    y2 = smooth(LFPChunks{i, 1}, 0.002);
    if any((y2-100)>-20)
        selectedChews = [selectedChews, i];
        plot(x1, y1*10, 'LineWidth', 2);
        hold on; box off;
        plot(x2, y2-100, 'LineWidth', 2);
        
        plotInd = plotInd + 1;
    end
    
    if plotInd > rows*cols
        figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
        rows = 2;
        cols = 5;
        plotInd = 1;    
    end
    
    
end

load('Z:\Qianyun\DCN\Data\Common\goodChews.mat');
jawDistanceMaxLength = 30;
LFPMaxLength = 7500;

% goodChewsJawDistanceMatrix = nan(length(selectedChews), jawDistanceMaxLength);
% goodChewsLFPMatrix = nan(length(selectedChews), LFPMaxLength);
temp = length(goodChewsJawDistanceMatrix);
for i = 1:length(selectedChews)
    ind = selectedChews(i);
    goodChewsJawDistanceMatrix(temp+i, :) = interp1(1:length(goodChewsJawDistance{ind, 1}), goodChewsJawDistance{ind, 1}, ...
        linspace(1, length(goodChewsJawDistance{ind, 1}), jawDistanceMaxLength), 'spline');
    
    goodChewsLFPMatrix(temp+i, :) = interp1(1:length(LFPChunks{ind, 1}), LFPChunks{ind, 1}, ...
        linspace(1, length(LFPChunks{ind, 1}), LFPMaxLength), 'spline');
    
end


% plot: overlay jaw movement and LFP rmsv
figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
yyaxis left
shadedErrorBar(linspace(0, 1, size(goodChewsJawDistanceMatrix, 2)), mean(goodChewsJawDistanceMatrix), std(goodChewsJawDistanceMatrix), ...
    'lineProps', {'-', 'Color', '#77AC30', 'LineWidth', 2});
box off;
ylabel('Jaw Distance(pixel)');


yyaxis right
% shadedErrorBar(linspace(0, 1, LFPMaxLength), mean(goodChewsLFPMatrix), std(goodChewsLFPMatrix), ...
%      'lineProps', {'-', 'Color', '#A2142F'});
shadedErrorBar(linspace(0, 1, LFPMaxLength), mean(goodChewsLFPMatrix), std(goodChewsLFPMatrix), ...
     'lineProps', {'c-', 'LineWidth', 2});

% for i = 1:length(goodChewsLFPMatrix)
%     plot(linspace(0, 1, LFPMaxLength), goodChewsLFPMatrix(i, :), '-'); hold on;
% end
box off;
ylabel('RMS Value');
xlabel('Jaw Phase');

legend('Jaw Movement', 'LFP');

legend('Jaw Movement', 'CWC', 'LFP');
ylabel('LFP RMSV or SpkRate (spk/s)')

save('Z:\Qianyun\DCN\Data\Common\goodChews.mat', 'goodChewsJawDistanceMatrix', 'goodChewsLFPMatrix');


%% Individual Chews - LFP Rasters overlay Jaw Distance


load('Z:\Qianyun\DCN\Data\Common\goodChews.mat');

figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
rows = 2;
cols = 5;
plotInd = 1;

jawDistanceMaxLength = 30;
LFPMaxLength = 7500;
peakInds = cell(length(goodChewsJawDistanceMatrix), 1);
for i = 1:length(goodChewsJawDistanceMatrix)
    subplot(rows, cols, plotInd);
    
    yyaxis left
    x1 = linspace(0, 1, jawDistanceMaxLength);
    y1 = goodChewsJawDistanceMatrix(i, :);
    plot(x1, y1, 'LineWidth', 2);
    
    yyaxis right
    x2 = linspace(0, 1, LFPMaxLength);
    y2 = goodChewsLFPMatrix(i, :);
    plot(x2, y2, 'LineWidth', 2); hold on;
    
    [~, locs] = findpeaks(y2, 'MinPeakHeight', max(y2)*0.8, 'MinPeakDistance', 150);
    peakInds{i, 1} = locs;
    
    plot(x2(locs), goodChewsLFPMatrix(i, locs), 'c.', 'MarkerSize', 15);
    
    plotInd = plotInd + 1;
    
    if plotInd > rows*cols
        figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
        rows = 2;
        cols = 5;
        plotInd = 1;    
    end   
end

temp = nan(length(peakInds), 1);
for i = 1:length(peakInds)
    temp(i) = peakInds{i, 1}(1);
end

[sorted, sortedIndex] = sort(temp);

% plot: 3 subplot figure - jaw movement, LFPCWC raster and JawLFPCWC curve
figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
rows = 2;
cols = 2;
plotInd = 1;

subplot(rows, cols, plotInd);
shadedErrorBar(linspace(0, 1, size(goodChewsJawDistanceMatrix, 2)), mean(goodChewsJawDistanceMatrix), std(goodChewsJawDistanceMatrix), ...
    'lineProps', {'-', 'Color', '#77AC30', 'LineWidth', 2});
box off;
ylabel('Jaw Distance(pixel)');
legend('Jaw Movement')
ax = gca;
ax.XAxis.Visible = 'off';

plotInd = 3;
subplot(rows, cols, plotInd);

for i = 1:length(spkInds)
    x1 = linspace(0, 1, LFPMaxLength);
    inds = spkInds{i, 1};
    x2 = x1(inds);
    if i == 1
        plot(x2, repmat(length(spkInds)+1-i, 1, length(x2)), '.', 'MarkerSize', 20, 'color', '#A2142F');
    else
        plot(x2, repmat(length(spkInds)+1-i, 1, length(x2)), '.', 'MarkerSize', 20, 'color', '#A2142F', 'HandleVisibility', 'off');
    end
    hold on; box off;
end

for i = 1:length(peakInds)
    x1 = linspace(0, 1, LFPMaxLength);
    inds = peakInds{i, 1};
    x2 = x1(inds);
    if i == 1
        plot(x2, repmat(length(peakInds)+1-i, 1, length(x2)), 'c.', 'MarkerSize', 20 );
    else
        plot(x2, repmat(length(peakInds)+1-i, 1, length(x2)), 'c.', 'MarkerSize', 20, 'HandleVisibility', 'off' );
    end
    hold on; box off;
end

legend('CWC', 'LFP')

% histogram(CWCRaster/7500); hold on;
% histogram(LFPRaster/7500); box off

xlabel('Jaw Phase');
ylabel('Individual Chew Number');
xlim([0, 1]);

% temp = [];
% L = 0;
% for i = 1:length(peakInds)
%     temp(L+1:L+length(peakInds{i, 1})) = peakInds{i, 1};
%     L = length(temp);
% end
% LFPRaster = temp;
% 
% temp = [];
% L = 0;
% for i = 1:length(spkInds)
%     temp(L+1:L+length(spkInds{i, 1})) = spkInds{i, 1};
%     L = length(temp);
% end
% CWCRaster = temp;

subplot(2, 2, [2, 4]);
yyaxis left
shadedErrorBar(linspace(0, 1, size(goodChewsJawDistanceMatrix, 2)), mean(goodChewsJawDistanceMatrix), std(goodChewsJawDistanceMatrix), ...
    'lineProps', {'-', 'Color', '#77AC30', 'LineWidth', 2});
box off;
ylabel('Jaw Distance(pixel)');


yyaxis right
shadedErrorBar(linspace(0, 1, LFPMaxLength), mean(goodChewsLFPMatrix), std(goodChewsLFPMatrix), ...
     'lineProps', {'c-', 'LineWidth', 2}); hold on;
shadedErrorBar(linspace(0, 1, size(spkRate, 2)), mean(spkRate), std(spkRate), ...
     'lineProps', {'-', 'Color', '#A2142F', 'LineWidth', 2}); hold on;
box off;
ylabel('LFP RMS Value or Spk/s');
xlabel('Jaw Phase');

legend('Jaw Movement', 'LFP', 'CWC');



%% Individual Chews - Process CWC 
ephysChannel = 23;
ephysThresh = -200;

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
    
    ephysVoltage = getVoltage(data.Data.Data(channelNum_OpenEphys(ephysChannel), :));
end
ephysTime = sessionEphysInfo.convertedEphysTimestamps;


goodChewsStartStopTime = trialInfo.goodChewsStartStopTime;
ephysStartInds = nan(length(goodChewsStartStopTime), 1);
ephysStopInds = nan(length(goodChewsStartStopTime), 1);
ephysChunks = cell(length(goodChewsStartStopTime), 1);
for i = 1:length(goodChewsStartStopTime)
    ephysStartInds(i) = find(ephysTime > goodChewsStartStopTime(i, 1), 1, 'first');
    ephysStopInds(i) = find(ephysTime < goodChewsStartStopTime(i, 2), 1, 'last');
    
    ephysChunks{i, 1} = getVoltage(data.Data.Data(channelNum_OpenEphys(ephysChannel), ephysStartInds(i, 1):ephysStopInds(i, 1))); 
    % ephysChunks{i, 1} = rmsv(ephysStartInds(i):ephysStopInds(i));
end

figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
rows = 2;
cols = 5;
plotInd = 1;

goodChewsJawDistance = trialInfo.goodChewsJawDistance;
for i = 1:length(ephysChunks)
    subplot(rows, cols, plotInd);
    
    x1 = linspace(goodChewsStartStopTime(i, 1), goodChewsStartStopTime(i, 2), length(goodChewsJawDistance{i, 1}));
    y1 = goodChewsJawDistance{i, 1};
    
    x2 = linspace(ephysTime(ephysStartInds(i)), ephysTime(ephysStopInds(i)), length(ephysChunks{i, 1}));
    y2 = ephysChunks{i, 1};
    
    plot(x1, y1*15, 'LineWidth', 2);
    hold on; box off;
    plot(x2, y2, 'LineWidth', 1);
    axis tight;
    
    plotInd = plotInd + 1;

    
    if plotInd > rows*cols
        figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
        rows = 2;
        cols = 5;
        plotInd = 1;    
    end
    
    
end

selectedChews = [1:18, 20:28, 30:32, 34:111, 113:156, 159:178];
figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
rows = 2;
cols = 5;
plotInd = 1;

for i = 1:length(selectedChews)
    subplot(rows, cols, plotInd);
    ind = selectedChews(i);
    x1 = linspace(goodChewsStartStopTime(ind, 1), goodChewsStartStopTime(ind, 2), length(goodChewsJawDistance{ind, 1}));
    y1 = goodChewsJawDistance{ind, 1};
    
    x2 = linspace(ephysTime(ephysStartInds(ind)), ephysTime(ephysStopInds(ind)), length(ephysChunks{ind, 1}));
    y2 = ephysChunks{ind, 1};
    
    [~, spikeInds, ~] = crossdet(y2, ephysThresh, 'thresholdDown');
    spkInds{i, 1} = spikeInds;
    
    plot(x1, y1*15, 'LineWidth', 2);
    hold on; box off;
    plot(x2, y2, 'LineWidth', 1);
    axis tight;
    plot(x2(spikeInds), repmat(ephysThresh, 1, length(spikeInds)), 'c.', 'MarkerSize', 10);
    
    plotInd = plotInd + 1;
    

    
    if plotInd > rows*cols
        figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
        rows = 2;
        cols = 5;
        plotInd = 1;    
    end
    
    
end

load('Z:\Qianyun\DCN\Data\Common\goodChews.mat');
% temp = length(goodChewsJawDistanceMatrix);
jawDistanceMaxLength = 30;
ephysMaxLength = 7500;

goodChewsJawDistanceMatrixCWC = nan(length(selectedChews), jawDistanceMaxLength);
goodChewsEphysMatrix = nan(length(selectedChews), ephysMaxLength);

spkInds = cell(1, 1);
spkRate = nan(1, 250);
tempInd = 1;
for i = 1:length(selectedChews)
    ind = selectedChews(i);
    goodChewsJawDistanceMatrixCWC(i, :) = interp1(1:length(goodChewsJawDistance{ind, 1}), goodChewsJawDistance{ind, 1}, ...
        linspace(1, length(goodChewsJawDistance{ind, 1}), jawDistanceMaxLength), 'spline');
    
    goodChewsEphysMatrix(i, :) = interp1(1:length(ephysChunks{ind, 1}), ephysChunks{ind, 1}, ...
        linspace(1, length(ephysChunks{ind, 1}), ephysMaxLength), 'spline');
    
    [~, spikeInds, ~] = crossdet(goodChewsEphysMatrix(i, :), ephysThresh, 'thresholdDown');
    
    
    if ~isempty(spkInds)
        spkInds{tempInd, 1} = spikeInds;
        spikeTimes = spikeInds/sessionEphysInfo.fs;
        [spkRate(tempInd, :), spikeRateTimes] = getFiringRate(spikeTimes, 'tLims', [1/sessionEphysInfo.fs, ephysMaxLength/sessionEphysInfo.fs]);
        tempInd = tempInd + 1;
    end
    
end


figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
yyaxis left
shadedErrorBar(linspace(0, 1, size(goodChewsJawDistanceMatrixCWC, 2)), mean(goodChewsJawDistanceMatrixCWC), std(goodChewsJawDistanceMatrixCWC), ...
    'lineProps', {'-', 'Color', '#77AC30', 'LineWidth', 2});
box off;
ylabel('Jaw Distance(pixel)');


yyaxis right
shadedErrorBar(linspace(0, 1, size(spkRate, 2)), mean(spkRate), std(spkRate), ...
     'lineProps', {'-', 'Color', '#A2142F', 'LineWidth', 2});

box off;
ylabel('SpkRate(Spk/s)');
xlabel('Jaw Phase');

legend('Jaw Movement', 'CWC');


save(fullfile(sessionFolder, 'CWCGoodChews.mat'), 'goodChewsJawDistanceMatrixCWC', 'goodChewsEphysMatrix', 'spkInds', 'spkRate');

