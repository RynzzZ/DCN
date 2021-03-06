function analyzeTrialUnit(session, varargin)

% This function breaks session data into each 'food' trial.
% (1) Save data into 'trialAnalyzed.mat'. 
%       - trial table: each row is one 'food' trial.
%       - chewingCrossCorr: each row is the cross corr b/w chewing jaw traces and the corresponding LFP for one 'food' trial. Length: [-maxLag:maxLag]
%       - allCrunchLFP: each row is the smoothed LFP centered around the big crunch for one 'food' trial. Length: [-s.crunchTimeWindow : s.crunchTimeWindow]*s.ephysFs
% (2) Plot series of crunch and chewing plots for each trial.
% (3) Plots are saved in the 'trialFigs' folder in the session folder.

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



% format ephys data into trial structure
if s.hasMic && s.hasJawTrace
    trialTotalCount = spike.totalFoodNum;
    trial = table();
    
    maxLag = 0.2*s.spkRateFs;
    chewingJawUnitCrossCorr = nan(trialTotalCount, maxLag*2+1);
    for i = 1:trialTotalCount
        
        fprintf('trial %d/%d \n', i, trialTotalCount);

        %%%%%%%%%%%%%%%%%%%%% Processing crunch %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get big crunch chunk 
        trial.foodTriggerTime(i) = spike.foodTriggerTimes(i);
        crunchStartTime = spike.foodTriggerTimes(i) - 0.1;
        crunchEndTime = crunchStartTime + s.crunchSearchTimeWindow;
        trial.trialStartTime(i) = crunchStartTime;
        
        crunchChunkMicStartInd = find(spike.micSignalTimes >= crunchStartTime, 1, 'first');
        crunchChunkMicEndInd = find(spike.micSignalTimes <= crunchEndTime, 1, 'last');
        crunchChunkMic = spike.micSignal(crunchChunkMicStartInd : crunchChunkMicEndInd);
        crunchChunkMic = highpass(crunchChunkMic, 100, 150000);
        crunchChunkMicRMSV = sqrt(movmean(crunchChunkMic.^2, 100));
        
        % get ephysData around the big crunch        
        crunchChunkEphysStartInd = find(spkRateTimes >= crunchStartTime, 1, 'first');
        crunchChunkEphysEndInd = find(spkRateTimes <= crunchEndTime, 1, 'last');
        crunchChunkEphysTimes1 = spkRateTimes(crunchChunkEphysStartInd:crunchChunkEphysEndInd);
        crunchChunkEphys1 = spkRates(crunchChunkEphysStartInd:crunchChunkEphysEndInd); % unit spkRate
        
        crunchChunkEphysStartInd = find(sessionEphysInfo.convertedEphysTimestamps >= crunchStartTime, 1, 'first');
        crunchChunkEphysEndInd = find(sessionEphysInfo.convertedEphysTimestamps <= crunchEndTime, 1, 'last');
        crunchChunkEphysTimes2 = sessionEphysInfo.convertedEphysTimestamps(crunchChunkEphysStartInd:crunchChunkEphysEndInd);
        crunchChunkEphys2 = unitVoltage(crunchChunkEphysStartInd:crunchChunkEphysEndInd); % unit raw data
        
        %%%%%%%%%%%%%%%%%%%%% Processing chewing %%%%%%%%%%%%%%%%%%%%%%%%%%%
        trialStartTime = spike.foodTriggerTimes(i) - 0.1;
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
        
        % get ephysData of the chewing period
        chewingChunkSpkRateStartInd = find(spkRateTimes >= chewingChunkStartTime, 1, 'first');
        chewingChunkSpkRateEndInd = find(spkRateTimes <= chewingChunkEndTime, 1, 'last');
        chewingChunkEphysTimes1 = spkRateTimes(chewingChunkSpkRateStartInd:chewingChunkSpkRateEndInd);
        chewingChunkEphys1 = spkRates(chewingChunkSpkRateStartInd:chewingChunkSpkRateEndInd); % unit spkRate
        
        chewingChunkEphysStartInd = find(sessionEphysInfo.convertedEphysTimestamps >= chewingChunkStartTime, 1, 'first');
        chewingChunkEphysEndInd = find(sessionEphysInfo.convertedEphysTimestamps <= chewingChunkEndTime, 1, 'last');
        chewingChunkEphysTimes2 = sessionEphysInfo.convertedEphysTimestamps(chewingChunkEphysStartInd:chewingChunkEphysEndInd);
        chewingChunkEphys2 = unitVoltage(chewingChunkEphysStartInd:chewingChunkEphysEndInd); % unit raw data
        
        % calculate the cross corr b/w jaw distance and unit spkRate
        cameraTimestamps = videoTracking.frameTimestamps(chewingChunkVideoStartInd:chewingChunkVideoEndInd);
        ephysTimestamps = spkRateTimes(chewingChunkSpkRateStartInd:chewingChunkSpkRateEndInd);
        chewingChunkJaw_interp = interp1(cameraTimestamps, chewingChunkJaw, ephysTimestamps);
        chewingChunkSpkRate_xcorr = chewingChunkEphys1(~isnan(chewingChunkJaw_interp));
        chewingChunkJaw_xcorr = chewingChunkJaw_interp(~isnan(chewingChunkJaw_interp));
        
        [c, lags] = xcorr(chewingChunkJaw_xcorr, chewingChunkSpkRate_xcorr, maxLag);
        chewingJawUnitCrossCorr(i, :) = c;
        [maxCorr, maxCorrLoc] = max(c);
        
        % save key variables into trial table
        trial.chewingChunkStartTime(i) = chewingChunkStartTime;
        trial.chewingChunkVideoStartInd(i) = chewingChunkVideoStartInd;
        trial.chewingChunkEphysStartInd(i) = chewingChunkEphysStartInd;
        trial.chewingChunkSpkRateStartInd(i) = chewingChunkSpkRateStartInd;
        trial.chewingChunkEndTime(i) = chewingChunkEndTime;
        trial.chewingChunkVideoEndInd(i) = chewingChunkVideoEndInd;
        trial.chewingChunkEphysEndInd(i) = chewingChunkEphysEndInd;
        trial.chewingChunkSpkRateEndInd(i) = chewingChunkSpkRateEndInd;
        trial.chewingChunkJawSpkRateMaxCorr(i) = maxCorr;
        trial.chewingChunkJawSpkRateMaxCorrLag(i) = (maxCorrLoc-maxLag)/s.spkRateFs; % unit: sec        
               
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % PLOT!
        if ~s.supressFigure
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
            title([session(1:8), '_', session(10:end) ' Trial ', num2str(i), ' first big crunch'], 'Interpreter','none');
            
            % plot mic signal RMSV
            plotInd = plotInd + 1;
            subplot(rows, cols, plotMatrix(plotInd));
            
            plot(spike.audioSignalTimes(crunchChunkMicStartInd:crunchChunkMicEndInd), crunchChunkMicRMSV, ...
                '-', 'Color', colorMatrix{plotInd});
            hold on; box off; axis tight;
            plot([trial.foodTriggerTime(i), trial.foodTriggerTime(i)], [min(crunchChunkMicRMSV), max(crunchChunkMicRMSV)], 'k-');
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
            plot([trial.foodTriggerTime(i), trial.foodTriggerTime(i)], [0, 1], 'k-');
            ylabel('confidence');
            ylim([0, 1]);
            legend('Tongue Existence','Location', 'northwest'); legend boxoff;
            h = gca; h.XAxis.Visible = 'off';
            
            plotInd = plotInd + 1;
            subplot(rows, cols, plotMatrix(plotInd));
            plot(videoTracking.frameTimestamps(crunchChunkVideoStartInd:crunchChunkVideoEndInd), jawTrace,...
                '-', 'Color', colorMatrix{plotInd});
            hold on; box off; axis tight;
            plot([trial.foodTriggerTime(i), trial.foodTriggerTime(i)], [min(jawTrace), max(jawTrace)], 'k-');
            ylabel('pixel');
            legend('Jaw Distance', 'Location', 'northwest'); legend boxoff;
            h = gca; h.XAxis.Visible = 'off';           
            
            % Unit spkRate
            plotInd = plotInd + 1;
            subplot(rows, cols, plotMatrix(plotInd));

            plot(crunchChunkEphysTimes1, crunchChunkEphys1, '-', 'Color', colorMatrix{plotInd});
            hold on; box off; axis tight;
            plot([trial.foodTriggerTime(i), trial.foodTriggerTime(i)], [min(crunchChunkEphys1), max(crunchChunkEphys1)], 'k-');
            ylabel('spikes/s');
            legend('spike rate', 'Location', 'northwest'); legend boxoff;
            h = gca; h.XAxis.Visible = 'off';
            
            % Unit Raw Data
            plotInd = plotInd + 1;
            subplot(rows, cols, plotMatrix(plotInd));
            
            plot(crunchChunkEphysTimes2, crunchChunkEphys2, '-', 'Color', colorMatrix{plotInd});
            hold on; box off; axis tight;
            plot([trial.foodTriggerTime(i), trial.foodTriggerTime(i)], [min(crunchChunkEphys2), max(crunchChunkEphys2)], 'k-');
            ylabel('microvolt');
            legend('raw data', 'Location', 'northwest'); legend boxoff;
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
            title([session(1:8), '_', session(10:end) ' Trial ', num2str(i), ' rhythmic chewing'], 'Interpreter','none');
            
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
            
            
            % LFP rmsv / Unit spkRates
            plotInd = plotInd + 1;
            subplot(rows, cols, plotMatrix(plotInd));
            
            plot(chewingChunkEphysTimes1, chewingChunkEphys1, '-', 'Color', colorMatrix{plotInd});
            box off; axis tight;
            ylabel('spikes/s');
            legend('spike rate', 'Location', 'northwest'); legend boxoff;
            h = gca; h.XAxis.Visible = 'off';
            
            % LFP rmsv / Unit Raw Data
            plotInd = plotInd + 1;
            subplot(rows, cols, plotMatrix(plotInd));
            
            plot(chewingChunkEphysTimes2, chewingChunkEphys2, '-', 'Color', colorMatrix{plotInd});
            hold on; box off; axis tight;
            ylabel('microvolt');
            legend('raw data', 'Location', 'northwest'); legend boxoff;
            xlabel('time (s)');
            
            % save the plot!
            saveas(gcf, fullfile(sessionFolder, 'trialFigs', ['trial', num2str(i), '_MicJawUnit_Ch',...
                num2str(bestChannel),'.png']));            
        end
    end   
    settings = s;
    save(fullfile(sessionFolder, 'trialAnalyzedUnit.mat'), 'trial', 'chewingJawUnitCrossCorr', 'settings');
    
    % plot the cross corr curve b/w chewing jaw distance and unit spkRate over all the trials in
    % this session
    figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
    x = [-maxLag:maxLag]/s.spkRateFs;
    y = mean(chewingJawUnitCrossCorr);
    err = std(chewingJawUnitCrossCorr)/sqrt(size(chewingJawUnitCrossCorr, 1));
    
    shadedErrorBar(x, y, err, 'lineProps', {'-', 'Color', '#0072BD'});
    hold on; axis tight;
    xlabel('time (s)');
    title([session, ' cross corr b/w chewing jaw trace & Unit SpkRate'], 'Interpreter','none');
    
    % plot the line for max cross corr and find its corresponding time point
    maxCorrLoc = find(y == max(y));
    maxCorrLoc = (maxCorrLoc - maxLag)/s.spkRateFs;
    yLimits = get(gca,'YLim');
    plot([maxCorrLoc, maxCorrLoc], [min(yLimits), max(yLimits)], 'r-');
    legend('chewing & unit spkRate cross corr', ['max corr loc = ', num2str(maxCorrLoc), ' s']);
    legend boxoff;
    
    % save the cross corr plot
    saveas(gcf, fullfile(sessionFolder, 'trialFigs', ['chewing jaw trace & Unit Ch', num2str(bestChannel), ' cross corr.png']));
    
    disp('Trial Analysis Finished!');
else
    disp('Crunches and chewings are NOT analyzed, because no good JawTrace or no good Mic');
end

end