function analyzeTrialLFP(session, varargin)

% This function breaks session data into each 'food' trial.
% (1) Save data into 'trialAnalyzed.mat'. 
%       - trial table: each row is one 'food' trial.
%       - chewingCrossCorr: each row is the cross corr b/w chewing jaw traces and the corresponding LFP for one 'food' trial. Length: [-maxLag:maxLag]
%       - allCrunchLFP: each row is the smoothed LFP centered around the big crunch for one 'food' trial. Length: [-s.crunchTimeWindow : s.crunchTimeWindow]*s.ephysFs
% (2) Plot series of crunch and chewing plots for each trial.
% (3) Plots are saved in the 'trialFigs' folder in the session folder.

% settings
s.hasLFP = true; % whether this session's recording contains a good LFP channel
s.hasMic = true; % whether this session contains good microphone recordings
s.hasJawTrace = true; % whether this session contains video tracking for jaw trace
s.supressFigure = false; % whether to only process the data but supress plotting figures
s.crunchSearchTimeWindow = 4; % sec
s.crunchTimeWindow = 0.02; % sec
s.chewingSearchTimeWindow = 10; % sec
s.chewingTimeWindow = 3; % sec
s.chewingSearchStepLength = 0.1; % sec
s.fs = 150000; % sampling rate for mic signal
s.ephysFs = 30000; % sampling rate for ephys


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
    load(fullfile(sessionFolder, 'spikeAnalyzed.mat'));
end

% load cameraAnalyzed.mat
if ~exist(fullfile(sessionFolder, 'cameraAnalyzed.mat'), 'file')
    error('cameraAnalyzed.mat NOT exist! Run analyzeSession.m first!');
else
    disp('loading cameraAnalyzed.mat...')
    load(fullfile(sessionFolder, 'cameraAnalyzed.mat'));
    if ~any(strcmp('jawDistance', videoTracking.Properties.VariableNames))
        warning('JawDistance does NOT exist in videoTracking! Will NOT analyze chewing LFP!');
        s.hasJawTrace = false;
    end
end


% get LFP channel data
ephysInfo = readtable(fullfile(rootFolder, 'Spreadsheets', 'ephysInfo.xlsx'));
LFPChannel = ephysInfo.LFP(strcmp(ephysInfo.session, session));
if isnan(LFPChannel) || LFPChannel == 0
    warning('No LFP channel in this session! Will NOT analyze crunch LFP and chewing LFP!');
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
    rmsv = sqrt(movmean(LFPVoltage.^2, 100));
end

% format LFP data into trial structure
if s.hasLFP && s.hasMic && s.hasJawTrace
    trialTotalCount = spike.totalFoodNum;
    trial = table();
    
    maxLag = 10000;
    chewingCrossCorr = nan(trialTotalCount, maxLag*2+1);
    allCrunchLFP = nan(trialTotalCount, s.crunchTimeWindow*2*s.ephysFs);
    allCrunchChunkLFP = nan(trialTotalCount, s.crunchSearchTimeWindow*s.ephysFs);
    allChewingChunkLFP = nan(trialTotalCount, s.chewingSearchTimeWindow*s.ephysFs);
    for i = 1:trialTotalCount
        
        fprintf('trial %d/%d \n', i, trialTotalCount);
        
        %%%%%%%%%%%%%%%%%%%%% Processing crunch %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % locate big crunch chunk
        trial.foodTriggerTime(i) = spike.foodTriggerTimes(i);
        crunchStartTime = spike.foodTriggerTimes(i) - 0.1;
        crunchEndTime = crunchStartTime + s.crunchSearchTimeWindow;
        trial.trialStartTime(i) = crunchStartTime;
        trial.trialEndTime(i) = crunchEndTime;
        
        % crunch chunk microphone
        crunchChunkMicStartInd = find(spike.micSignalTimes >= crunchStartTime, 1, 'first');
        crunchChunkMicEndInd = find(spike.micSignalTimes <= crunchEndTime, 1, 'last');
        crunchChunkMic = spike.micSignal(crunchChunkMicStartInd : crunchChunkMicEndInd);
        crunchChunkMic = highpass(crunchChunkMic, 100, 150000);
        crunchChunkMicRMSV = sqrt(movmean(crunchChunkMic.^2, 100));
        
        % crunch chunk ephys(LFP)
        crunchChunkEphysStartInd = find(sessionEphysInfo.convertedEphysTimestamps >= crunchStartTime, 1, 'first');
        crunchChunkEphysEndInd = find(sessionEphysInfo.convertedEphysTimestamps <= crunchEndTime, 1, 'last');
        trial.crunchChunkEphysStartInd(i) =  crunchChunkEphysStartInd;
        trial.crunchChunkEphysEndInd(i) = crunchChunkEphysEndInd;
        crunchChunkLFP = rmsv(crunchChunkEphysStartInd:crunchChunkEphysEndInd);
        smoothedCrunchChunkLFP = smooth(crunchChunkLFP, 0.002);
        
        allCrunchChunkLFP(i, 1:length(crunchChunkLFP)) = crunchChunkLFP;
        
        % locate and analyze the big crunch
        crunchMicThreshold = min(0.2, max(crunchChunkMic)*0.8);
        bigCrunchTime = spike.micSignalTimes(find(crunchChunkMic(1+s.fs*0.1:end) >= crunchMicThreshold, 1, 'first') + crunchChunkMicStartInd);
        if any(bigCrunchTime)
            bigCrunch = true;
            trial.CrunchTime(i) = bigCrunchTime;
            [trial, allCrunchLFP] = analyzeBigCrunch(bigCrunchTime, i, smoothedCrunchChunkLFP, trial, sessionEphysInfo, s.crunchTimeWindow);
        else
            bigCrunch = false;
            disp('Did not detect the big crunch in mic recording, cannot analyze big crunch for this trial...')
        end
     
        
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
        
        % get LFP RMSV of the chewing period
        chewingChunkEphysStartInd = find(sessionEphysInfo.convertedEphysTimestamps >= chewingChunkStartTime, 1, 'first');
        chewingChunkEphysEndInd = find(sessionEphysInfo.convertedEphysTimestamps <= chewingChunkEndTime, 1, 'last');
        chewingChunkLFP = rmsv(chewingChunkEphysStartInd:chewingChunkEphysEndInd);
        smoothedChewingChunkLFP = smooth(chewingChunkLFP, 0.009);
        
        % calculate the cross corr b/w jaw distance and LFP
        cameraTimestamps = videoTracking.frameTimestamps(chewingChunkVideoStartInd:chewingChunkVideoEndInd);
        ephysTimestamps = sessionEphysInfo.convertedEphysTimestamps(chewingChunkEphysStartInd:chewingChunkEphysEndInd);
        chewingChunkJaw_interp = interp1(cameraTimestamps, chewingChunkJaw, ephysTimestamps);
        chewingChunkLFP_xcorr = smoothedChewingChunkLFP(~isnan(chewingChunkJaw_interp));
        chewingChunkJaw_xcorr = chewingChunkJaw_interp(~isnan(chewingChunkJaw_interp));
        
        [c, lags] = xcorr(chewingChunkJaw_xcorr, chewingChunkLFP_xcorr, maxLag);
        chewingCrossCorr(i, :) = c;
        [maxCorr, maxCorrLoc] = max(c);
        
        % save key variables into trial table
        trial.chewingChunkStartTime(i) = chewingChunkStartTime;
        trial.chewingChunkVideoStartInd(i) = chewingChunkVideoStartInd;
        trial.chewingChunkEphysStartInd(i) = chewingChunkEphysStartInd;
        trial.chewingChunkEndTime(i) = chewingChunkEndTime;
        trial.chewingChunkVideoEndInd(i) = chewingChunkVideoEndInd;
        trial.chewingChunkEphysEndInd(i) = chewingChunkEphysEndInd;
        trial.chewingChunkMaxCorr(i) = maxCorr; % normalized
        trial.chewingChunkMaxCorrLag(i) = (maxCorrLoc-maxLag)/s.ephysFs; % unit: sec
        
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
            
            % LFP rmsv (less smoothed)
            plotInd = plotInd + 1;
            subplot(rows, cols, plotMatrix(plotInd));
            
            plot(sessionEphysInfo.convertedEphysTimestamps(crunchChunkEphysStartInd:crunchChunkEphysEndInd),...
                smoothedCrunchChunkLFP, '-', 'Color', colorMatrix{plotInd});
            hold on; box off; axis tight;
            plot([trial.foodTriggerTime(i), trial.foodTriggerTime(i)], [min(smoothedCrunchChunkLFP), max(smoothedCrunchChunkLFP)], 'k-');
            if bigCrunch
                plot(sessionEphysInfo.convertedEphysTimestamps(trial.crunchLFP_MaxValLoc(i)), trial.crunchLFP_MaxVal(i), 'r.', 'MarkerSize', 13);
                plot(sessionEphysInfo.convertedEphysTimestamps(trial.crunchLFP_MinValLoc(i)), trial.crunchLFP_MinVal(i), 'y.', 'MarkerSize', 13);
            end
            ylabel('rms value');
            legend('LFP RMSV, span = 0.002', 'Location', 'northwest'); legend boxoff;
            h = gca; h.XAxis.Visible = 'off';
            
            % LFP rmsv (more smoothed)
            plotInd = plotInd + 1;
            subplot(rows, cols, plotMatrix(plotInd));
            smoothedCrunchChunkLFP2 = smooth(crunchChunkLFP, 0.009);
            
            plot(sessionEphysInfo.convertedEphysTimestamps(crunchChunkEphysStartInd:crunchChunkEphysEndInd),...
                smoothedCrunchChunkLFP2, '-', 'Color', colorMatrix{plotInd});
            hold on; box off; axis tight;
            plot([trial.foodTriggerTime(i), trial.foodTriggerTime(i)], [min(smoothedCrunchChunkLFP2), max(smoothedCrunchChunkLFP2)], 'k-');
            ylabel('rms value');
            legend('LFP RMSV, span = 0.009', 'Location', 'northwest'); legend boxoff;
            xlabel('time (s)');
            
            %%%%%%%%%%% series of plots for chewing period %%%%%%%%%%%
            plotMatrix = [2, 4, 6, 8, 12, 10];
            colorMatrix = {' ', '#0072BD', 	'#EDB120', 	'#77AC30', 	'#7E2F8E', 	'#A2142F'};
                   
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
            
            
            % LFP rmsv (less smoothed)
            plotInd = plotInd + 1;
            subplot(rows, cols, plotMatrix(plotInd));
            
            smoothedChewingChunkLFP2 = smooth(chewingChunkLFP, 0.002);
            
            plot(sessionEphysInfo.convertedEphysTimestamps(chewingChunkEphysStartInd:chewingChunkEphysEndInd),...
                smoothedChewingChunkLFP2, '-', 'Color', colorMatrix{plotInd});
            box off; axis tight;
            ylabel('rms value');
            xlabel('time (s)');
            legend('LFP RMSV, span = 0.002', 'Location', 'northwest'); legend boxoff;
           
            
            % LFP rmsv (more smoothed)
            plotInd = plotInd + 1;
            subplot(rows, cols, plotMatrix(plotInd));
            
            plot(sessionEphysInfo.convertedEphysTimestamps(chewingChunkEphysStartInd:chewingChunkEphysEndInd),...
                smoothedChewingChunkLFP, '-', 'Color', colorMatrix{plotInd});
            hold on; box off; axis tight;
            ylabel('rms value');
            legend('LFP RMSV, span = 0.009', 'Location', 'northwest'); legend boxoff;            
            h = gca; h.XAxis.Visible = 'off';
            
            saveas(gcf, fullfile(sessionFolder, 'trialFigs', ['trial ', num2str(i), '.png']));
        end
        
    end   
    
    save(fullfile(sessionFolder, 'trialAnalyzedLFP.mat'), 'trial', 'chewingCrossCorr', 'allCrunchLFP');
    
    % plot the cross corr curve b/w chewing jaw distance and LFP over all the trials in
    % this session
    figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
    x = [-maxLag:maxLag]/s.ephysFs;
    y = mean(chewingCrossCorr);
    err = std(chewingCrossCorr)/sqrt(size(chewingCrossCorr, 1));
    
    shadedErrorBar(x, y, err, 'lineProps', {'-', 'Color', '#0072BD'});
    hold on; axis tight;
    xlabel('time (s)');
    title([session, ' Trial ', num2str(i), ' cross corr b/w chewing jaw trace & LFP'], 'Interpreter','none');
    
    % plot the line for max cross corr and find its corresponding time point
    maxCorrLoc = find(y == max(y));
    maxCorrLoc = (maxCorrLoc - maxLag)/s.ephysFs;
    yLimits = get(gca,'YLim');
    plot([maxCorrLoc, maxCorrLoc], [min(yLimits), max(yLimits)], 'r-');
    legend('chewing & LFP cross corr', ['max corr loc = ', num2str(maxCorrLoc), ' s']);
    legend boxoff;
    
    % save the cross corr plot
    saveas(gcf, fullfile(sessionFolder, 'trialFigs', 'chewing jaw trace & LFP cross corr.png'));
    
    
    % plot the avg LFP RMVS across all the crunches
    figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
    x = (1:length(crunchLFP))/s.ephysFs - s.crunchTimeWindow;
    y = mean(allCrunchLFP);
    err = std(allCrunchLFP)/sqrt(size(allCrunchLFP, 1));
    
    shadedErrorBar(x, y, err, 'lineProps', {'-', 'Color', '#0072BD'});
    hold on; axis tight;
    xlabel('time (s)');
    title([session, ' Trial ', num2str(i), ' avg crunch LFP'], 'Interpreter','none');
    
    % save the avg crunch LFP plot
    saveas(gcf, fullfile(sessionFolder, 'trialFigs', 'session avg crunch LFP.png'));
    
    disp('Trial Analysis Finished!');
else   
    disp('Crunches and chewings are NOT analyzed, because no good LFP or no good JawTrace or no good Mic');
end

    function  [trial, allCrunchLFP] = analyzeBigCrunch(bigCrunchTime, i, smoothedCrunchChunkLFP, trial, sessionEphysInfo, crunchTimeWindow)
        
        crunchChunkEphysStartInd = trial.crunch
        
        % get LFP around the big crunch
        crunchEphysStartTime = bigCrunchTime - crunchTimeWindow;
        crunchEphysEndTime = bigCrunchTime + crunchTimeWindow;
        crunchEphysStartInd = find(sessionEphysInfo.convertedEphysTimestamps >= crunchEphysStartTime, 1, 'first');
        crunchEphysEndInd = find(sessionEphysInfo.convertedEphysTimestamps <= crunchEphysEndTime, 1, 'last');
        
        crunchLFP = smoothedCrunchChunkLFP(crunchEphysStartInd-crunchChunkEphysStartInd : crunchEphysEndInd-crunchChunkEphysStartInd);
        allCrunchLFP(i, :) = crunchLFP;
        
        % calculate LFP rms value for the big crunch
        [maxVal, maxloc] = max(crunchLFP);
        [minVal, minloc] = min(crunchLFP);
        
        maxPeakLoc = maxloc + crunchEphysStartInd;
        minPeakLoc = minloc + crunchEphysStartInd;
        
        trial.crunchLFP_rmsVal(i) = maxVal - minVal;
        trial.crunchLFP_MaxVal(i) = maxVal;
        trial.crunchLFP_MinVal(i) = minVal;
        trial.crunchLFP_MaxValLoc(i) = maxPeakLoc;
        trial.crunchLFP_MaxValTime(i) = sessionEphysInfo.convertedEphysTimestamps(maxPeakLoc);
        trial.crunchLFP_MinValLoc(i) = minPeakLoc;
        trial.crunchLFP_MinValTime(i) = sessionEphysInfo.convertedEphysTimestamps(minPeakLoc);
        
        
    end




end