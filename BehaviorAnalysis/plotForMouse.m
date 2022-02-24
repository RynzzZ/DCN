%% extract jaw movement and overlay it with DCN LFP

session = '20211111_002';

% initializations
rootFolder = 'Z:\Qianyun\DCN\';
gitFolder = 'D:\DCN_Project\Github\DCN\';
sessionFolder = fullfile(rootFolder, 'Data', session);

% load analyzed mat files
load(fullfile(sessionFolder, 'spikeAnalyzed.mat'));
load(fullfile(sessionFolder, 'cameraAnalyzed.mat'));
load(fullfile(sessionFolder, 'sessionEphysInfo.mat'));

% get food trigger times
randomSampleNum = 5;
randomFoodTriggerTimes = sort(randsample(spike.foodTriggerTimes, randomSampleNum));

for i = 1:randomSampleNum
    
    figure('Color', 'white', 'position', get(0,'ScreenSize')); clf; hold on;
    
    startTime = randomFoodTriggerTimes(i) - 0.3;
    stopTime = randomFoodTriggerTimes(i) + 10;
    
    % plot mic signal
    startInd = find(spike.micSignalTimes > startTime, 1, 'first');
    stopInd = find(spike.micSignalTimes < stopTime, 1, 'last');
    ind = startInd:stopInd;
    
    plot(spike.micSignalTimes(ind) - randomFoodTriggerTimes(i), spike.micSignal(ind)/max(spike.micSignal(ind))); 
    hold on
    
    % plot jaw movement
    yoffset = 2;
    
    startInd = find(videoTracking.frameTimestamps > startTime, 1, 'first');
    stopInd = find(videoTracking.frameTimestamps < stopTime, 1, 'last');
    ind = startInd:stopInd;
        
    plot(videoTracking.frameTimestamps(ind) - randomFoodTriggerTimes(i),...
        videoTracking.jawDistance(ind)/max(videoTracking.jawDistance(ind)) + yoffset);
    
    
    % plot tongue movement
    yoffset = 3;
    
    startInd = find(videoTracking.frameTimestamps > startTime, 1, 'first');
    stopInd = find(videoTracking.frameTimestamps < stopTime, 1, 'last');
    ind = startInd:stopInd;
    
    plot(videoTracking.frameTimestamps(ind) - randomFoodTriggerTimes(i), videoTracking.tongueExist(ind) + yoffset);
    plot(videoTracking.frameTimestamps(ind) - randomFoodTriggerTimes(i), videoTracking.tongue_confidence(ind) + yoffset);
    
    % plot DCN LFP
    yoffset = 1;
    
    LFPChannel = 30;
    trace = showChannelsTimeChunk(session, startTime, stopTime, 'selectedChannels', LFPChannel, 'dataOnly', true);
    plot(linspace(startTime, stopTime, length(trace)) - randomFoodTriggerTimes(i), ((trace - min(trace))/(max(trace) - min(trace))) + yoffset);
    
    xlim([startTime , stopTime] - randomFoodTriggerTimes(i));
    ylim([-1, 4]);
    
    xlabel('time (sec, 0 = food dispenser trigger)');
    legend('mic signal', 'jaw movement', 'tongue', 'tongue confidence', 'DCN LFP');
    ax = gca;
    ax.YAxis.Visible = 'off'; % remove y-axis
end

%% trial plot: crunch analyses, plot Units and/or LFP and/or behavior

session = '20211216_000';
cwcChannel = 23;
opcChannel = 28;
LFPChannel = [];

% settings
s.crunchSearchTimeWindow = 4; % sec
s.crunchTimeWindow = 0.02; % sec
s.chewingSearchTimeWindow = 10; % sec
s.chewingTimeWindow = 3; % sec
s.chewingSearchStepLength = 0.1; % sec
s.fs = 150000; % sampling rate for mic signal
s.ephysFs = 30000; % sampling rate for ephys
s.spkRateFs = 1000; % sampling rate for unit spkRate
s.supressFigure = false;


% initializations
rootFolder = 'Z:\Qianyun\DCN\';
gitFolder = 'D:\DCN_Project\Github\DCN\';
sessionFolder = fullfile(rootFolder, 'Data', session);
saveFolder = fullfile(rootFolder, 'Data', session, 'trialFigs', 'crunch');

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

getVoltage = @(data) double(data)*sessionEphysInfo.bitVolts; % extract voltage from memmapfile, converting to votlage, highpassing, and only return specific channel
contFiles = dir(fullfile(sessionEphysInfo.ephysFolder, '*.continuous'));
data = memmapfile(fullfile(sessionEphysInfo.ephysFolder, [contFiles(1).name(1:end-12), 's.dat']), ...
    'Format', {'int16', [sessionEphysInfo.channelNum sessionEphysInfo.smps], 'Data'}, 'Writable', false);

cwcVoltage = getVoltage(data.Data.Data(channelNum_OpenEphys(cwcChannel), :));
opcVoltage = getVoltage(data.Data.Data(channelNum_OpenEphys(opcChannel), :));
if any(LFPChannel)
    LFPVoltage = getVoltage(data.Data.Data(channelNum_OpenEphys(LFPChannel), :));
    rmsv = sqrt(movmean(LFPVoltage.^2, 100));
end

% process trials
trialTotalCount = spike.totalFoodNum;
for i = 1:trialTotalCount
        
        fprintf('trial %d/%d \n', i, trialTotalCount);

        %%%%%%%%%%%%%%%%%%%%% Processing crunch %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get big crunch chunk 
        crunchStartTime = spike.foodTriggerTimes(i) - 0.1;
        crunchEndTime = crunchStartTime + s.crunchSearchTimeWindow;
        
        crunchChunkMicStartInd = find(spike.micSignalTimes >= crunchStartTime, 1, 'first');
        crunchChunkMicEndInd = find(spike.micSignalTimes <= crunchEndTime, 1, 'last');
        crunchChunkMic = spike.micSignal(crunchChunkMicStartInd : crunchChunkMicEndInd);
        crunchChunkMic = highpass(crunchChunkMic, 100, 150000);
        crunchChunkMicRMSV = sqrt(movmean(crunchChunkMic.^2, 100));
        
        % get ephysData around the big crunch        
        crunchChunkEphysStartInd = find(sessionEphysInfo.convertedEphysTimestamps >= crunchStartTime, 1, 'first');
        crunchChunkEphysEndInd = find(sessionEphysInfo.convertedEphysTimestamps <= crunchEndTime, 1, 'last');
        crunchChunkEphysTime = sessionEphysInfo.convertedEphysTimestamps(crunchChunkEphysStartInd:crunchChunkEphysEndInd);
        
        if any(cwcChannel)
            cwcCrunchChunkData = cwcVoltage(crunchChunkEphysStartInd:crunchChunkEphysEndInd); % cwc raw data
        end
        
        if any(opcChannel)
            opcCrunchChunkData = opcVoltage(crunchChunkEphysStartInd:crunchChunkEphysEndInd); % opc raw data
        end
        
        if any(LFPChannel)
            LFPCrunchChunkData = rmsv(crunchChunkEphysStartInd:crunchChunkEphysEndInd); % LFP rmsv data
        end
        
        % get the jaw trace
        crunchChunkVideoStartInd = find(videoTracking.frameTimestamps >= crunchStartTime, 1, 'first');
        crunchChunkVideoEndInd = find(videoTracking.frameTimestamps <= crunchEndTime, 1, 'last');
        jawTrace = videoTracking.jawDistance(crunchChunkVideoStartInd:crunchChunkVideoEndInd);
    
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % PLOT!
        if ~s.supressFigure
            figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
            rows = 6; cols = 1;
            
            % series of plots for big crunch
            colorMatrix = {' ', '#0072BD', 		'#7E2F8E', 	'#A2142F', '#EDB120', 	'#77AC30', };
            
            
            % plot mic signal spectrogram
            plotInd = 1;
            subplot(rows, cols, plotInd);
            
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
            subplot(rows, cols, plotInd);
            
            plot(spike.audioSignalTimes(crunchChunkMicStartInd:crunchChunkMicEndInd), crunchChunkMicRMSV, ...
                '-', 'Color', colorMatrix{plotInd});
            hold on; box off; axis tight;
            plot([spike.foodTriggerTimes(i), spike.foodTriggerTimes(i)], [min(crunchChunkMicRMSV), max(crunchChunkMicRMSV)], 'k-');
            ylabel('RMS Value');
            legend('Mic RMSV', 'Food Trigger', 'Location', 'northwest'); legend boxoff;
            h = gca; h.XAxis.Visible = 'off';
            
            
            
            
            % opc Raw Data
            if any(opcChannel)
                plotInd = plotInd + 1;
                subplot(rows, cols, plotInd);
                
                plot(crunchChunkEphysTime, opcCrunchChunkData, '-', 'Color', colorMatrix{plotInd});
                hold on; box off; axis tight;
                plot([spike.foodTriggerTimes(i), spike.foodTriggerTimes(i)], [min(opcCrunchChunkData), max(opcCrunchChunkData)], 'k-');
                ylabel('microvolt');
                legend('opc raw data', 'Location', 'northwest'); legend boxoff;
                h = gca; h.XAxis.Visible = 'off';
            end
            
            % cwc raw data
            if any(cwcChannel)
                plotInd = plotInd + 1;
                subplot(rows, cols, plotInd);
                
                plot(crunchChunkEphysTime, cwcCrunchChunkData, '-', 'Color', colorMatrix{plotInd});
                hold on; box off; axis tight;
                plot([spike.foodTriggerTimes(i), spike.foodTriggerTimes(i)], [min(cwcCrunchChunkData), max(cwcCrunchChunkData)], 'k-');
                ylabel('microvolt');
                legend('cwc raw data', 'Location', 'northwest'); legend boxoff;
                h = gca; h.XAxis.Visible = 'off';
            end
            
            
            % plot tongue (confidence) trace 
            plotInd = plotInd + 1;
            subplot(rows, cols, plotInd);
            plot(videoTracking.frameTimestamps(crunchChunkVideoStartInd:crunchChunkVideoEndInd), ...
                videoTracking.tongue_confidence(crunchChunkVideoStartInd:crunchChunkVideoEndInd), ...
                '-', 'Color', colorMatrix{plotInd});
            hold on; box off; axis tight;
            plot([spike.foodTriggerTimes(i), spike.foodTriggerTimes(i)], [0, 1], 'k-');
            ylabel('confidence');
            ylim([0, 1]);
            legend('Tongue Existence','Location', 'northwest'); legend boxoff;
            h = gca; h.XAxis.Visible = 'off';
            
            % plot jaw trace
            plotInd = plotInd + 1;
            subplot(rows, cols, plotInd);
            plot(videoTracking.frameTimestamps(crunchChunkVideoStartInd:crunchChunkVideoEndInd), jawTrace, ...
                '-', 'Color', colorMatrix{plotInd});
            hold on; box off; axis tight;
            plot([spike.foodTriggerTimes(i), spike.foodTriggerTimes(i)], [min(jawTrace), max(jawTrace)], 'k-');
            ylabel('pixel');
            legend('Jaw Distance', 'Location', 'northwest'); legend boxoff;
            xlabel('time (s)');  
            
            % save figs
            saveas(gcf, fullfile(saveFolder, ['trial_', num2str(i), '_crunch.png']))
        end
        
        
end



















