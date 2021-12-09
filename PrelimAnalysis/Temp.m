%% temp script 

session = '20211111_002';
trialID = 9;

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


