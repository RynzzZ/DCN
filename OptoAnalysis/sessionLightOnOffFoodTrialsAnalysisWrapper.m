%% session 20220610_000, opto effect, ch17(output cell)
session = '20220610_000';
outputChannel = 17;
LFPChannel = 32;
cwcChannel = [];
timeWindow = [-1, 15]; % sec
chewingStartTimeShift = 4; % sec
chewingSearchLength = 15; % sec
chewingChunkLength = 4; % sec
chewingSearchStepLength = 0.1; % sec
cwcThresh = [];
outputThresh = -250;
toneMapCues = 'KLHUOP';
fs = 30000;
smoothLevel = 0.05;

% initializations
rootFolder = 'Z:\Qianyun\DCN\';
gitFolder = 'D:\DCN_Project\Github\DCN\';
sessionFolder = fullfile(rootFolder, 'Data', session);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%processing data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load spikeAnalyzed.mat
if ~exist(fullfile(sessionFolder, 'spikeAnalyzed.mat'), 'file')
    error('spikeAnalyzed.mat NOT exist! Run analyzeSession.m first!');
else
    disp('loading spikeAnalyzed.mat...')
    load(fullfile(sessionFolder, 'spikeAnalyzed.mat'), 'spike');
    
    % get sampling frequency for audio sig channel and mic channel in spike
    s.audioFs = round(length(spike.audioSignalTimes)/(spike.audioSignalTimes(end) - spike.audioSignalTimes(1)));
    s.micFs = round(length(spike.micSignalTimes)/(spike.micSignalTimes(end) - spike.micSignalTimes(1)));
    disp('finish loading spikeAnalyzed.mat');
end


if ~exist(fullfile(sessionFolder, 'cameraAnalyzed.mat'), 'file')
    warning('cameraAnalyzed.mat NOT exist! Run analyzeSession.m first!');
else
    disp('loading cameraAnalyzed.mat...')
    load(fullfile(sessionFolder, 'cameraAnalyzed.mat'));
    
    jawY = videoTracking.jaw_lower_y(videoTracking.jaw_lower_confidence>0.8);
    frameTimestamps = videoTracking.frameTimestamps(videoTracking.jaw_lower_confidence>0.8);
    disp('finish loading spikeAnalyzed.mat');
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
   
    
    % determine target channel
    if contains(sessionEphysInfo.ephysFolder, '101')
        outputChannel = channelNum_OpenEphys(outputChannel);
        cwcChannel = channelNum_OpenEphys(cwcChannel);
        LFPChannel = channelNum_OpenEphys(LFPChannel);
    end
    
    LFPData = getVoltage(data.Data.Data(LFPChannel, :));
    LFPDataFiltered = medfilt1(LFPData,20);
    LFPrmsv = sqrt(movmean(LFPDataFiltered.^2, 100)); % turn into rmsv
    
%     startTime = 544;
%     stopTime = 549;
%     startInd = find(ephysTime >= startTime, 1, 'first');
%     stopInd = find(ephysTime <= stopTime, 1, 'last');
%     
%     figure;
%     rows = 4;
%     subplot(rows, 1, 1);
%     plot(LFPData(1, startInd:stopInd));
%     subplot(rows, 1, 2);
%     plot(LFPDataFiltered(1, startInd:stopInd));
%     subplot(rows, 1, 3);
%     plot(LFPrmsv(1, startInd:stopInd));
%     subplot(rows, 1, 4);
%     plot(smooth(LFPrmsv(1, startInd:stopInd), 0.02));
%        
    
    disp('finish loading ephys data');
    
end

%% quality check
figure;
subplot(2, 1, 1)
startTime = 144.5;
endTime = 150.5;


startInd = find(ephysTime >= startTime, 1, 'first');
endInd = find(ephysTime <= endTime, 1, 'last');

cameraStartInd = find(frameTimestamps >= startTime, 1, 'first');
cameraStopInd = find(frameTimestamps <= endTime, 1, 'last');

% if ~isempty(cwcChannel)    
%     x = ephysTime(startInd:endInd);
%     y = getVoltage(data.Data.Data(channelNum_OpenEphys(cwcChannel), startInd:endInd));
%     
%     plot(x, y); hold on;
%     plot([x(1), x(end)], [cwcThresh, cwcThresh], '-r');
%     axis tight
% end

if ~isempty(outputChannel)
    
    subplot(2, 1, 1)
    x = ephysTime(startInd:endInd);
    y = getVoltage(data.Data.Data(outputChannel, startInd:endInd));
    
    plot(x, y); hold on;
    plot([x(1), x(end)], [outputThresh, outputThresh], '-r'); hold on;
    yyaxis right
    x2 = frameTimestamps(cameraStartInd:cameraStopInd);
    y2 = jawY(cameraStartInd:cameraStopInd);
    plot(x2, y2, 'LineWidth', 1.5);
    axis tight
    title('example opto food trial')
    
    subplot(2, 1, 2)
    startTime = 342;
    endTime = 349;
    
    
    startInd = find(ephysTime >= startTime, 1, 'first');
    endInd = find(ephysTime <= endTime, 1, 'last');
    
    cameraStartInd = find(frameTimestamps >= startTime, 1, 'first');
    cameraStopInd = find(frameTimestamps <= endTime, 1, 'last');
    
    
    x = ephysTime(startInd:endInd);
    y = getVoltage(data.Data.Data(outputChannel, startInd:endInd));
    
    plot(x, y); hold on;
    plot([x(1), x(end)], [outputThresh, outputThresh], '-r'); hold on;
    yyaxis right
    yyaxis right
    x2 = frameTimestamps(cameraStartInd:cameraStopInd);
    y2 = jawY(cameraStartInd:cameraStopInd);
    plot(x2, y2, 'LineWidth', 1.5);
    axis tight
    
    title('example normal food trial')
end

%% preprocess the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Find start time for all the normal and opto food trials (normal =
% light off; opto = light on)
% 2. Put normal and opto food trial ephys data + LFP data into cells
% 3. Calculate output cell & cwc firing rate for each normal and opto food trial

close all
 
% locate opto train start
optoTrialNum = 50;
optoFoodTrialStartTime = nan(optoTrialNum, 1);
normalFoodTrialStartTime = nan(optoTrialNum, 1);
outputEphysData = cell(optoTrialNum, 2);
cwcEphysData = cell(optoTrialNum, 2);
cwcSpikeInds = cell(optoTrialNum, 1);
outputSpkRates = cell(optoTrialNum, 2);
outputSpkTimes = cell(optoTrialNum, 2);
LFPEphysData = cell(optoTrialNum, 2);
ephysTimeData = cell(optoTrialNum, 2);

normalFoodColumn = 1;
optoFoodColum = 2;

normalFoodTrialCount = 0;
optoFoodTrialCount = 0;

slideWindow = 3000;
keyboardString = convertCharsToStrings(spike.keyboardInput);
GInds = strfind(keyboardString, 'G');

for i = 1:length(spike.keyboardInput)-1
    
    disp(['i = ', num2str(i)]);
    if strcmp(spike.keyboardInput(i), 'F')
        if sum(spike.optoTrainTimes >= spike.keyboardTimes(i) & spike.optoTrainTimes <= spike.keyboardTimes(i) + 15) > 450
            foodTrigerTime = spike.keyboardTimes(i);
            
            % determine light trigger time
            timeDiff = spike.keyboardTimes(GInds) - spike.keyboardTimes(i);
            [~, ind] = min(abs(timeDiff));
            lightTriggerTime = spike.keyboardTimes(GInds(ind));
                       
            optoFoodTrialCount = optoFoodTrialCount + 1;
            disp(['opto food trial: ', num2str(optoFoodTrialCount)]);
            optoFoodTrialStartTime(optoFoodTrialCount) = foodTrigerTime;
            
            % determine trial start and stop time
            ephysStartTime = foodTrigerTime + timeWindow(1);
            ephysStopTime = ephysStartTime + timeWindow(2);
    
            
            % get ephys time inds
            ephysStartInd = find(ephysTime >= ephysStartTime, 1, 'first');
            ephysStopInd = find(ephysTime <= ephysStopTime, 1, 'last');
            
            % cwc cell
            if ~isempty(cwcChannel)
                % voltage data
                cwcEphysData{optoFoodTrialCount, optoFoodColum} = getVoltage(data.Data.Data(channelNum_OpenEphys(cwcChannel), ephysStartInd:ephysStopInd));
                
%                 % get cwc spike times
%                 [~, cwcSpikeInds{optoFoodTrialCount, optoFoodColum}, ~] = crossdet(cwcEphysData{count, 1}, cwcThresh, 'thresholdDown');
            end
            
            % LFP
            LFPtemp = LFPrmsv(ephysStartInd:ephysStopInd);
            LFPsmoothed = smooth(LFPtemp, smoothLevel);
            LFPEphysData{optoFoodTrialCount, optoFoodColum} = LFPsmoothed;
            ephysTimeData{optoFoodTrialCount, optoFoodColum} = ephysTime(ephysStartInd:ephysStopInd);
            
            % output cell
            % voltage data
            outputEphysData{optoFoodTrialCount, optoFoodColum} = getVoltage(data.Data.Data(outputChannel, ephysStartInd:ephysStopInd));
            
            % get output cell firing rate
            [~, outputSpikeInds, ~] = crossdet(outputEphysData{optoFoodTrialCount, optoFoodColum}, outputThresh, 'thresholdDown');
            spikes = zeros(size(outputEphysData{optoFoodTrialCount, optoFoodColum}));
            spikes(outputSpikeInds) = 1;
            
            temp = 1:slideWindow/2:length(spikes);
            FR = nan(size(spikes));
            
            for j = 1:length(temp)
                ind = temp(j);
                startInd = ind - slideWindow/2;
                stopInd = ind + slideWindow/2;
                
                if startInd <= 0
                    startInd = 1;
                end
                
                if stopInd > length(spikes)
                    stopInd = length(spikes);
                end
                
                FR(j) = sum(spikes(startInd:stopInd)) / ((stopInd - startInd + 1)/fs);
            end
            
            outputSpkRates{optoFoodTrialCount, optoFoodColum} = FR(~isnan(FR));
            Time = linspace(timeWindow(1), ephysStopTime - ephysStartTime + timeWindow(1), length(spikes));
            outputSpkTimes{optoFoodTrialCount, optoFoodColum} = Time(temp);
            
            figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
            % LFP rmsv data
            x = ephysTimeData{optoFoodTrialCount, optoFoodColum};
            y1 = LFPEphysData{optoFoodTrialCount, optoFoodColum};
            y1 = (y1 - min(y1))/range(y1)*2; % normalize to [0, 1] range
            plot(x, y1);
            hold on;
            % output cell raw data
            y2 = outputEphysData{optoFoodTrialCount, optoFoodColum};
            y2 = ((y2 - min(y2)) / range(y2))*2 + 2;
            plot(x, y2);
            % spikes
            y3 = ones(size(outputSpikeInds))*2.7;
            plot(x(outputSpikeInds), y3, '.r', 'MarkerSize', 8);
            % opto train (if any)
            x4 = find(spike.optoTrainTimes >= ephysStartTime & spike.optoTrainTimes <= ephysStopTime);
            y4 = ones(size(x4))*4;
            ylim([0, 5]);
            plot(spike.optoTrainTimes(x4), y4, '.m', 'MarkerSize', 8);
            
            legend('LFPrmsv', 'output cell raw data');
            title(['opto food trial ', num2str(optoFoodTrialCount)]);
            
                     
        else
            normalFoodTrialCount = normalFoodTrialCount + 1;
            disp(['normal food trial: ', num2str(normalFoodTrialCount)]);
            
            % determine trial start and stop time
            foodTrigerTime = spike.keyboardTimes(i);
            normalFoodTrialStartTime(normalFoodTrialCount) = foodTrigerTime;
            ephysStartTime = foodTrigerTime + timeWindow(1);
            ephysStopTime = ephysStartTime + timeWindow(2);
            
            % get ephys time inds
            ephysStartInd = find(ephysTime >= ephysStartTime, 1, 'first');
            ephysStopInd = find(ephysTime <= ephysStopTime, 1, 'last');
            
            % cwc cell
            if ~isempty(cwcChannel)
                % voltage data
                cwcEphysData{normalFoodTrialCount, normalFoodColumn} = getVoltage(data.Data.Data(cwcChannel, ephysStartInd:ephysStopInd));
                 
%                 % get cwc spike times
%                 [~, cwcSpikeInds{count, 1}, ~] = crossdet(cwcEphysData{count, 1}, cwcThresh, 'thresholdDown');
            end
            
            % LFP
            LFPtemp = LFPrmsv(ephysStartInd:ephysStopInd);
            LFPsmoothed = smooth(LFPtemp, 0.01);
            LFPEphysData{normalFoodTrialCount, normalFoodColumn} = LFPsmoothed;
            ephysTimeData{normalFoodTrialCount, normalFoodColumn} = ephysTime(ephysStartInd:ephysStopInd);
            
            % output cell
            % voltage data
            outputEphysData{normalFoodTrialCount, normalFoodColumn} = getVoltage(data.Data.Data(outputChannel, ephysStartInd:ephysStopInd));
            
            % get output cell firing rate
            [~, outputSpikeInds, ~] = crossdet(outputEphysData{normalFoodTrialCount, normalFoodColumn}, outputThresh, 'thresholdDown');
            spikes = zeros(size(outputEphysData{normalFoodTrialCount, normalFoodColumn}));
            spikes(outputSpikeInds) = 1;
            
            temp = 1:slideWindow/2:length(spikes);
            FR = nan(size(spikes));
            
            for j = 1:length(temp)
                ind = temp(j);
                startInd = ind - slideWindow/2;
                stopInd = ind + slideWindow/2;
                
                if startInd <= 0
                    startInd = 1;
                end
                
                if stopInd > length(spikes)
                    stopInd = length(spikes);
                end
                
                FR(j) = sum(spikes(startInd:stopInd)) / ((stopInd - startInd + 1)/fs);
            end
            
            outputSpkRates{normalFoodTrialCount, normalFoodColumn} = FR(~isnan(FR));
            Time = linspace(timeWindow(1), ephysStopTime - ephysStartTime + timeWindow(1), length(spikes));
            outputSpkTimes{normalFoodTrialCount, normalFoodColumn} = Time(temp);
            
            figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
            % LFP rmsv data
            x = ephysTimeData{normalFoodTrialCount, normalFoodColumn};
            y1 = LFPEphysData{normalFoodTrialCount, normalFoodColumn};
            y1 = (y1 - min(y1))/range(y1)*2; % normalize to [0, 1] range
            plot(x, y1);
            hold on;
            % output cell raw data
            y2 = outputEphysData{normalFoodTrialCount, normalFoodColumn};
            y2 = ((y2 - min(y2)) / range(y2))*2 + 2;
            plot(x, y2);
            % spikes
            y3 = ones(size(outputSpikeInds))*2.7;
            plot(x(outputSpikeInds), y3, '.r', 'MarkerSize', 8);
            % opto train (if any)
            x4 = find(spike.optoTrainTimes >= ephysStartTime & spike.optoTrainTimes <= ephysStopTime);
            y4 = ones(size(x4))*4;
            ylim([0, 5]);
            plot(spike.optoTrainTimes(x4), y4, '.m', 'MarkerSize', 8);
            
            legend('LFPrmsv', 'output cell raw data');
            title(['normal food trial ', num2str(normalFoodTrialCount)]);
            
            
        end
    end
end

%% process the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Algorithm to pick the best rhythmic chewing periods in each food trial
% 2. Calculate the phase (according to jaw movement) 
% 3. Save data into the optoTrialsChewingConstruct and
% normalTrialsChewingConstruct


% load cameraAnalyzed.mat
if ~exist(fullfile(sessionFolder, 'cameraAnalyzed.mat'), 'file')
    error('cameraAnalyzed.mat NOT exist! Run analyzeSession.m first!');
else
    disp('loading cameraAnalyzed.mat...')
    load(fullfile(sessionFolder, 'cameraAnalyzed.mat'));
    if ~any(strcmp('jawDistance', videoTracking.Properties.VariableNames))
        error('JawDistance does NOT exist in videoTracking! Run analyzeSession.m first!');
    end
    disp('finish loading cameraAnalyzed.mat');
end

% calculate good chewing periods & its jaw phase for opto food trials
% [optoTrialsChewingStruct] = getGoodChewingPeriodsJaw(videoTracking.jawDistance, videoTracking.frameTimestamps, rmmissing(optoFoodTrialStartTime));
optoLFPEphysData = nan(optoFoodTrialCount, length(LFPEphysData{1, 1}));
optoEphysTimeData = nan(optoFoodTrialCount, length(LFPEphysData{1, 1}));
for i = 1:optoFoodTrialCount
    optoLFPEphysData(i, :) = LFPEphysData{i, optoFoodColum}'; 
    optoEphysTimeData(i, :) = ephysTimeData{i, optoFoodColum};
end
[optoTrialsChewingStruct] = getGoodChewingPeriodsLFP(optoLFPEphysData, optoEphysTimeData, rmmissing(optoFoodTrialStartTime));

% calculate good chewing periods & its jaw phase for normal food trials
% [normalTrialsChewingStruct] = getGoodChewingPeriodsJaw(videoTracking.jawDistance, videoTracking.frameTimestamps, rmmissing(normalFoodTrialStartTime));
normalLFPEphysData = nan(normalFoodTrialCount, length(LFPEphysData{1, 1}));
normalEphysTimeData = nan(normalFoodTrialCount, length(LFPEphysData{1, 1}));
for i = 1:optoFoodTrialCount
    normalLFPEphysData(i, :) = LFPEphysData{i, normalFoodColumn}'; 
    normalEphysTimeData(i, :) = ephysTimeData{i, normalFoodColumn};
end
[normalTrialsChewingStruct] = getGoodChewingPeriodsLFP(LFPEphysData, ephysTimeData, rmmissing(normalFoodTrialStartTime));

%% if automatic algorithm fails, use this part of code to do it manually!
% manually define the start and stop time for good chewing chunk
% selectedOptoFoodTrials = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14];
% optoChewingChunkStartTimes = [43.8, 129.7, 236.1, 314.1, 436.2, 563.1, 705.8, 879.6, 1414.8, 1506.6, 1762.7, 2855.9, 2878.4, 3033.1];
% optoChewingChunkStopTimes = [45.1, 130.9, 237.8, 317.8, 438.7, 566, 708.3, 881.1, 1415.8, 1508.1, 1764.3, 2859.2, 2882, 3035.2];

selectedOptoFoodTrials = [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 13, 14];
optoChewingChunkStartTimes = [43.8, 129.7, 236.1, 314.1, 436.2, 563.1, 705.8, 879.6, 1506.6, 2855.9, 2878.4, 3033.1];
optoChewingChunkStopTimes = [45.1, 130.9, 237.8, 317.8, 438.7, 566, 708.3, 881.1, 1508.1, 2859.2, 2882, 3035.2];


optoTrialsChewingStruct.selectedOptoTrials = selectedOptoFoodTrials;
optoTrialsChewingStruct.optoChewingChunkStartTimes = optoChewingChunkStartTimes;
optoTrialsChewingStruct.optoChewingChunkStopTimes = optoChewingChunkStopTimes;


selectedNormalFoodTrials = [1, 2, 3, 4, 6, 7, 8, 9, 11, 13];
normalChewingChunkStartTimes = [15.1, 166, 262, 369.8, 662.2, 766.2, 939.7, 1440.9, 1794.7, 2917.9, ];
normalChewingChunkStopTimes = [18.0 169, 264.8, 373.1, 665.75, 770.3, 942.5, 1443, 1797.1, 2919.7, ];

normalTrialsChewingStruct.selectedNormalTrials = selectedNormalFoodTrials;
normalTrialsChewingStruct.normalChewingChunkStartTimes = normalChewingChunkStartTimes;
normalTrialsChewingStruct.normalChewingChunkStopTimes = normalChewingChunkStopTimes;

%% quality check for selected chewing chunks (plotted with LFP or Jaw)

flag = 'Jaw';
ephysChannel = outputChannel;
ephysThresh = outputThresh;
LFPSmoothFactor = 0.02;

for i = 1:length(optoChewingChunkStartTimes)
    
    figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
    % plot opto trial
    subplot(2, 1, 1)
    
    if i <= length(optoChewingChunkStartTimes)
        startTime = optoChewingChunkStartTimes(i);
        endTime = optoChewingChunkStopTimes(i);
        
        startInd = find(ephysTime >= startTime, 1, 'first');
        endInd = find(ephysTime <= endTime, 1, 'last');
        
        if strcmp(flag, 'Jaw')
            cameraStartInd = find(frameTimestamps >= startTime, 1, 'first');
            cameraStopInd = find(frameTimestamps <= endTime, 1, 'last');
        end
        
        % plot ephys raw data
        x = ephysTime(startInd:endInd);
        y = getVoltage(data.Data.Data(ephysChannel, startInd:endInd));
        plot(x, y); hold on;
        plot([x(1), x(end)], [ephysThresh, ephysThresh], '-r');
        hold on; axis tight;
        title(['opto trial ', num2str(i)]);
        
        % plot LFP/Jaw trace
        yyaxis right
        if strcmp(flag, 'Jaw')
            x2 = frameTimestamps(cameraStartInd:cameraStopInd);
            y2 = jawY(cameraStartInd:cameraStopInd);
            plot(x2, y2, 'LineWidth', 1.5);
            axis tight;
            
        elseif strcmp(flag, 'LFP')
            y2 = LFPrmsv(startInd:endInd);
            y2 = smooth(y2, LFPSmoothFactor);
            plot(x, y2);
            axis tight;
        end
    end
    
    % plot normal food trials
    subplot(2, 1, 2)
    
    if i <= length(normalChewingChunkStartTimes)
     
        startTime = normalChewingChunkStartTimes(i);
        endTime = normalChewingChunkStopTimes(i);
        
        startInd = find(ephysTime >= startTime, 1, 'first');
        endInd = find(ephysTime <= endTime, 1, 'last');
        
        if strcmp(flag, 'Jaw')
            cameraStartInd = find(frameTimestamps >= startTime, 1, 'first');
            cameraStopInd = find(frameTimestamps <= endTime, 1, 'last');
        end
        
        % plot ephys raw data
        x = ephysTime(startInd:endInd);
        y = getVoltage(data.Data.Data(ephysChannel, startInd:endInd));
        plot(x, y); hold on;
        plot([x(1), x(end)], [ephysThresh, ephysThresh], '-r');
        hold on; axis tight;
        title(['normal trial ', num2str(i)]);
        
        % plot LFP/Jaw trace
        yyaxis right
        if strcmp(flag, 'Jaw')
            x2 = frameTimestamps(cameraStartInd:cameraStopInd);
            y2 = jawY(cameraStartInd:cameraStopInd);
            plot(x2, y2, 'LineWidth', 1.5);
            axis tight;
            
        elseif strcmp(flag, 'LFP')
            y2 = LFPrmsv(startInd:endInd);
            y2 = smooth(y2, LFPSmoothFactor);
            plot(x, y2);
            axis tight;
        end
    end
    
    
end

%% data processing for phase histogram for opto food trials

close all

flag = 'Jaw'; % use 'LFP' or 'Jaw', this is to determine which signal to use for calculating the phase of chewing cycles.

selectedOptoFoodTrials = optoTrialsChewingStruct.selectedOptoTrials;

selectedChewingChunkLFPSmoothed = cell(length(selectedOptoFoodTrials), 1);
selectedChewingChunkLFPTime = cell(length(selectedOptoFoodTrials), 1);
selectedChewsStartTime = nan(500, 1);
selectedChewsStopTime = nan(500, 1);
selectedChewsLFPSmoothed = cell(500, 1);
selectedChewsJaw = cell(500, 1);
selectedChewsJawTime = cell(500, 1);
selectedChewsCount = 1;

optoChewingChunkStartTimes = optoTrialsChewingStruct.optoChewingChunkStartTimes;
optoChewingChunkStopTimes = optoTrialsChewingStruct.optoChewingChunkStopTimes;


% select good chews from good chewing chunk
for i = 1:length(selectedOptoFoodTrials)
    
    switch(flag)
        case 'LFP'
            % get chewing chunk start and stop time
            chewingChunkStartTime = optoChewingChunkStartTimes(i);
            chewingChunkStopTime = optoChewingChunkStopTimes(i);
            
            % get chewing chunk LFP
            inds = ephysTime >= chewingChunkStartTime & ephysTime <= chewingChunkStopTime;
            chewingChunkTime = ephysTime(inds);
            chewingChunkLFPSmoothed = smooth(LFPrmsv(inds), smoothLevel);
            
            selectedChewingChunkLFPSmoothed{i, 1} = chewingChunkLFPSmoothed;
            selectedChewingChunkLFPTime{i, 1} = chewingChunkTime;
            
            % calculate LFP phase for each chewing chunk
            temp = hilbert(chewingChunkLFPSmoothed);
            chewingChunkLFPPhase = angle(temp);
            chewingChunkLFPPhase = smooth(chewingChunkLFPPhase, 0.015);
            
            [maxVals, locs] = findpeaks(chewingChunkLFPPhase);
            [~, minlocs] = findpeaks(-chewingChunkLFPPhase);
            
            % filter out bad chews, only keep good chews
            temp = chewingChunkLFPPhase(3000:length(chewingChunkLFPPhase)-3000);
            minlocs(chewingChunkLFPPhase(minlocs) >= mean(temp)) = [];
            
            avgChewLength = mean(diff(minlocs));
            minlocsBackup = minlocs;
            inds = diff(minlocs) >= 12000 | diff(minlocs) <= 3000;
            minlocs(inds) = [];
            
            for j = 1:length(minlocs)
                if j < length(minlocs)
                    selectedChewsStartTime(selectedChewsCount, 1) = chewingChunkTime(minlocs(j));
                    selectedChewsStopTime(selectedChewsCount, 1) = chewingChunkTime(minlocsBackup(find(minlocsBackup == minlocs(j))+1));
                    selectedChewsLFPSmoothed{selectedChewsCount, 1} = chewingChunkLFPSmoothed(minlocs(j):minlocsBackup(find(minlocsBackup == minlocs(j))+1));
                    selectedChewsCount = selectedChewsCount + 1;
                end
            end
            
            % sanity check - phase calculation and chews selection
            figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
            subplot(2, 1, 1);
            plot(chewingChunkTime, chewingChunkLFPSmoothed); hold on;
            plot(chewingChunkTime(locs), chewingChunkLFPSmoothed(locs), 'r.', 'MarkerSize', 10)
            plot(chewingChunkTime(minlocs), chewingChunkLFPSmoothed(minlocs), 'b.', 'MarkerSize', 10)
            box off; axis tight;
            title(['trial ', num2str(selectedOptoFoodTrials(i))]);
            
            subplot(2, 1, 2);
            plot(chewingChunkTime, chewingChunkLFPPhase, '-r');
            box off; axis tight; hold on;
            plot(chewingChunkTime(locs), chewingChunkLFPPhase(locs), 'r.', 'MarkerSize', 10)
            plot(chewingChunkTime(minlocs), chewingChunkLFPPhase(minlocs), 'b.', 'MarkerSize', 10);
            
            
            
        case 'Jaw'
            % get chewing chunk start and stop time
            chewingChunkStartTime = optoChewingChunkStartTimes(i);
            chewingChunkStopTime = optoChewingChunkStopTimes(i);
            
            % get chewing chunk Jaw trace
            inds = frameTimestamps >= chewingChunkStartTime & frameTimestamps <= chewingChunkStopTime;
            chewingChunkJawTime = frameTimestamps(inds);
            chewingChunkJaw = jawY(inds);
            
            % calculate JAW phase for each chewing chunk
            temp = hilbert(chewingChunkJaw);
            chewingChunkJawPhase = angle(temp);
            chewingChunkJawPhase = smooth(chewingChunkJaw, 0.015);
            
            [maxVals, locs] = findpeaks(chewingChunkJawPhase);
            [~, minlocs] = findpeaks(-chewingChunkJawPhase);
            
            
            % filter out bad chews, only keep good chews
            temp = chewingChunkJawPhase(15:length(chewingChunkJawPhase)-15);
            minlocs(chewingChunkJawPhase(minlocs) >= mean(temp)) = [];
            
            avgChewLength = mean(diff(minlocs));
            minlocsBackup = minlocs;
            inds = diff(minlocs) >= 40 | diff(minlocs) <= 10;
            minlocs(inds) = [];
            
            for j = 1:length(minlocs)
                if j < length(minlocs)
                    selectedChewsStartTime(selectedChewsCount, 1) = chewingChunkJawTime(minlocs(j));
                    selectedChewsStopTime(selectedChewsCount, 1) = chewingChunkJawTime(minlocsBackup(find(minlocsBackup == minlocs(j))+1));
                    selectedChewsJaw{selectedChewsCount, 1} = chewingChunkJaw(minlocs(j):minlocsBackup(find(minlocsBackup == minlocs(j))+1));
                    selectedChewsCount = selectedChewsCount + 1;
                end
            end
            
            % sanity check - phase calculation and chews selection
            figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
            subplot(2, 1, 1);
            plot(chewingChunkJawTime, chewingChunkJaw); hold on;
            plot(chewingChunkJawTime(locs), chewingChunkJaw(locs), 'r.', 'MarkerSize', 10)
            plot(chewingChunkJawTime(minlocs), chewingChunkJaw(minlocs), 'b.', 'MarkerSize', 10)
            box off; axis tight;
            title(['trial ', num2str(selectedOptoFoodTrials(i))]);
            
            subplot(2, 1, 2);
            plot(chewingChunkJawTime, chewingChunkJawPhase, '-r');
            box off; axis tight; hold on;
            plot(chewingChunkJawTime(locs), chewingChunkJawPhase(locs), 'r.', 'MarkerSize', 10)
            plot(chewingChunkJawTime(minlocs), chewingChunkJawPhase(minlocs), 'b.', 'MarkerSize', 10);
    end
end

selectedChewsStartTime = rmmissing(selectedChewsStartTime);
selectedChewsStopTime = rmmissing(selectedChewsStopTime);

% give warnings if selected chews are way too long
if any(selectedChewsStopTime - selectedChewsStartTime > 0.3)
    warning('Contain chews > 0.3s length, need double check!');
end


% detect spikes & calculate firing rate for all good chews!
selectedChewsOutputFR = cell(length(selectedChewsStartTime), 1);
selectedChewsLFPInterp = nan(length(selectedChewsStartTime), 0.3*fs);
selectedChewsJawInterp = nan(length(selectedChewsStartTime), 0.3*fs);
selectedChewsOutputFRInterp = nan(length(selectedChewsStartTime), 0.3*fs);
selectedChewsNormalizedSpikeInds = nan(length(selectedChewsStartTime)*20, 1);
tempInd = 1;
slideWindow = 900;
figure;
for i = 1:length(selectedChewsStartTime)
    
    selectedChewEphysStartInd = find(ephysTime >= selectedChewsStartTime(i), 1, 'first');
    selectedChewEphysStopInd = find(ephysTime <= selectedChewsStopTime(i), 1, 'last');
    
    outputEphysChewData = getVoltage(data.Data.Data(outputChannel, selectedChewEphysStartInd:selectedChewEphysStopInd));
    
    % get output cell firing rate
    [~, outputSpikeInds, ~] = crossdet(outputEphysChewData, outputThresh, 'thresholdDown');
    spikes = zeros(size(outputEphysChewData));
    spikes(outputSpikeInds) = 1;
    
    temp = 1:slideWindow/4:length(spikes);
    FR = nan(size(spikes));
    
    for j = 1:length(temp)
        ind = temp(j);
        startInd = ind - slideWindow/4;
        stopInd = ind + slideWindow/4;
        
        if startInd <= 0
            startInd = 1;
        end
        
        if stopInd > length(spikes)
            stopInd = length(spikes);
        end
        
        FR(j) = sum(spikes(startInd:stopInd)) / ((stopInd - startInd + 1)/fs);
    end
    FR = FR(~isnan(FR));
    selectedChewsOutputFR{i, 1} = FR;
    x = outputSpikeInds/length(outputEphysChewData);
    y = ones(size(x))*i;
    plot(x, y, '.', 'MarkerSize', 8); hold on;
    
    selectedChewsNormalizedSpikeInds(tempInd:tempInd+length(x)-1, 1) = x;
    tempInd = tempInd + length(x);
    
    % interp LFP and output cell FR onto the same time strech
    switch flag
        case 'LFP'
            x = selectedChewsLFPSmoothed{i, 1};
            selectedChewsLFPInterp(i, :) = interp1(1:length(x), x, ...
                linspace(1, length(x), size(selectedChewsLFPInterp, 2)));
            
        case 'Jaw'
            x = selectedChewsJaw{i, 1};
            selectedChewsJawInterp(i, :) = interp1(1:length(x), x, ...
                linspace(1, length(x), size(selectedChewsJawInterp, 2)));
    end
    
    selectedChewsOutputFRInterp(i, :) = interp1(1:length(FR), FR, ...
        linspace(1, length(FR), size(selectedChewsOutputFRInterp, 2)));
    
end

selectedChewsNormalizedSpikeInds = selectedChewsNormalizedSpikeInds(~isnan(selectedChewsNormalizedSpikeInds));


% save data to opto trial struct
switch flag
    case 'LFP'
        optoTrialsChewingStruct.selectedChewsLFPSmoothed = selectedChewsLFPSmoothed;
        optoTrialsChewingStruct.selectedChewsLFPInterp = selectedChewsLFPInterp;
        
    case 'Jaw'
        optoTrialsChewingStruct.selectedChewsJawInterp = selectedChewsJawInterp;      
end
optoTrialsChewingStruct.selectedChewsStartTime = selectedChewsStartTime;
optoTrialsChewingStruct.selectedChewsStopTime = selectedChewsStopTime;
optoTrialsChewingStruct.selectedChewsOutputFR = selectedChewsOutputFR;
optoTrialsChewingStruct.selectedChewsOutputFRInterp = selectedChewsOutputFRInterp;
optoTrialsChewingStruct.selectedChewsNormalizedSpikeInds = selectedChewsNormalizedSpikeInds;
optoTrialsChewingStruct.outputThresh = outputThresh;

% calculate the vector strength for normal food trials
optoTrialChewsPhase = optoTrialsChewingStruct.selectedChewsNormalizedSpikeInds*2*pi;
x = sum(cos(optoTrialChewsPhase));
y = sum(sin(optoTrialChewsPhase));
optoVS = sqrt(x^2 + y^2) / length(optoTrialChewsPhase);

optoTrialsChewingStruct.VS = optoVS;


%% data processing for phase histogram for normal food trials

flag = 'Jaw'; % use 'LFP' or 'Jaw', this is to determine which signal to use for calculating the phase of chewing cycles.

selectedNormalFoodTrials = normalTrialsChewingStruct.selectedNormalTrials;

selectedChewingChunkLFPSmoothed = cell(length(selectedOptoFoodTrials), 1);
selectedChewingChunkLFPTime = cell(length(selectedOptoFoodTrials), 1);
selectedChewsStartTime = nan(500, 1);
selectedChewsStopTime = nan(500, 1);
selectedChewsLFPSmoothed = cell(500, 1);
selectedChewsJaw = cell(500, 1);
selectedChewsJawTime = cell(500, 1);
selectedChewsCount = 1;

normalChewingChunkStartTimes = normalTrialsChewingStruct.normalChewingChunkStartTimes;
normalChewingChunkStopTimes = normalTrialsChewingStruct.normalChewingChunkStopTimes;

% select good chews from good chewing chunk
for i = 1:length(selectedNormalFoodTrials)

    chewingChunkStartTime = normalChewingChunkStartTimes(i);
    chewingChunkStopTime = normalChewingChunkStopTimes(i);
    
    switch flag
        case 'LFP'
            inds = ephysTime >= chewingChunkStartTime & ephysTime <= chewingChunkStopTime;
            chewingChunkTime = ephysTime(inds);
            chewingChunkLFPSmoothed = smooth(LFPrmsv(inds), smoothLevel);
            
            selectedChewingChunkLFPSmoothed{i, 1} = chewingChunkLFPSmoothed;
            selectedChewingChunkLFPTime{i, 1} = chewingChunkTime;
            
            
            % calculate LFP phase for each chewing chunk
            temp = hilbert(chewingChunkLFPSmoothed);
            chewingChunkLFPPhase = angle(temp);
            chewingChunkLFPPhase = smooth(chewingChunkLFPPhase, 0.015);
            
            [maxVals, locs] = findpeaks(chewingChunkLFPPhase);
            [~, minlocs] = findpeaks(-chewingChunkLFPPhase);           
            
            % filter out bad chews, only keep good chews
            temp = chewingChunkLFPPhase(3000:length(chewingChunkLFPPhase)-3000);
            minlocs(chewingChunkLFPPhase(minlocs) >= mean(temp)) = [];
            
            avgChewLength = mean(diff(minlocs));
            minlocsBackup = minlocs;
            inds = find(diff(minlocs) >= 12000 | diff(minlocs) <= 3000);
            minlocs(inds) = [];
            
            for j = 1:length(minlocs)
                if j < length(minlocs)
                    selectedChewsStartTime(selectedChewsCount, 1) = chewingChunkTime(minlocs(j));
                    selectedChewsStopTime(selectedChewsCount, 1) = chewingChunkTime(minlocsBackup(find(minlocsBackup == minlocs(j))+1));
                    selectedChewsLFPSmoothed{selectedChewsCount, 1} = chewingChunkLFPSmoothed(minlocs(j):minlocsBackup(find(minlocsBackup == minlocs(j))+1));
                    selectedChewsCount = selectedChewsCount + 1;
                end
            end
            
            % sanity check - phase calculation and chews selection
            figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
            subplot(2, 1, 1);
            plot(chewingChunkTime, chewingChunkLFPSmoothed); hold on;
            plot(chewingChunkTime(locs), chewingChunkLFPSmoothed(locs), 'r.', 'MarkerSize', 10)
            plot(chewingChunkTime(minlocs), chewingChunkLFPSmoothed(minlocs), 'b.', 'MarkerSize', 10)
            box off; axis tight;
            title(['trial ', num2str(selectedNormalFoodTrials(i))]);
            
            subplot(2, 1, 2);
            plot(chewingChunkTime, chewingChunkLFPPhase, '-r');
            box off; axis tight; hold on;
            plot(chewingChunkTime(locs), chewingChunkLFPPhase(locs), 'r.', 'MarkerSize', 10)
            plot(chewingChunkTime(minlocs), chewingChunkLFPPhase(minlocs), 'b.', 'MarkerSize', 10);
            
            
        case 'Jaw'
            inds = frameTimestamps >= chewingChunkStartTime & frameTimestamps <= chewingChunkStopTime;
            chewingChunkJawTime = frameTimestamps(inds);
            chewingChunkJaw = jawY(inds);
            
            % calculate LFP phase for each chewing chunk
            temp = hilbert(chewingChunkJaw);
            chewingChunkJawPhase = angle(temp);
            chewingChunkJawPhase = smooth(chewingChunkJawPhase, 0.015);
            
            [maxVals, locs] = findpeaks(chewingChunkJawPhase);
            [~, minlocs] = findpeaks(-chewingChunkJawPhase);
            
                               
            % filter out bad chews, only keep good chews
            temp = chewingChunkJawPhase(15:length(chewingChunkJawPhase)-15);
            minlocs(chewingChunkJawPhase(minlocs) >= mean(temp)) = [];
            minlocs = minlocs - round(mean(diff(minlocs))*0.4);
            
            avgChewLength = mean(diff(minlocs));
            minlocsBackup = minlocs;
            inds = diff(minlocs) >= 40 | diff(minlocs) <= 10;
            minlocs(inds) = [];
            minlocs(minlocs < 0) = [];
            
            for j = 1:length(minlocs)
                if j < length(minlocs)
                    selectedChewsStartTime(selectedChewsCount, 1) = chewingChunkJawTime(minlocs(j));
                    selectedChewsStopTime(selectedChewsCount, 1) = chewingChunkJawTime(minlocsBackup(find(minlocsBackup == minlocs(j))+1));
                    selectedChewsJaw{selectedChewsCount, 1} = chewingChunkJaw(minlocs(j):minlocsBackup(find(minlocsBackup == minlocs(j))+1));
                    selectedChewsCount = selectedChewsCount + 1;
                end
            end
            
            % sanity check - phase calculation and chews selection
            figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
            subplot(2, 1, 1);
            plot(chewingChunkJawTime, chewingChunkJaw); hold on;
            plot(chewingChunkJawTime(locs), chewingChunkJaw(locs), 'r.', 'MarkerSize', 10)
            plot(chewingChunkJawTime(minlocs), chewingChunkJaw(minlocs), 'b.', 'MarkerSize', 10)
            box off; axis tight;
            title(['trial ', num2str(selectedNormalFoodTrials(i))]);
            
            subplot(2, 1, 2);
            plot(chewingChunkJawTime, chewingChunkJawPhase, '-r');
            box off; axis tight; hold on;
            plot(chewingChunkJawTime(locs), chewingChunkJawPhase(locs), 'r.', 'MarkerSize', 10)
            plot(chewingChunkJawTime(minlocs), chewingChunkJawPhase(minlocs), 'b.', 'MarkerSize', 10);
    end
end

selectedChewsStartTime = rmmissing(selectedChewsStartTime);
selectedChewsStopTime = rmmissing(selectedChewsStopTime);

if any(selectedChewsStopTime - selectedChewsStartTime > 0.3)
    warning('Contain chews > 0.3s length, need double check!');
end

selectedChewsOutputFR = cell(length(selectedChewsStartTime), 1);
selectedChewsLFPInterp = nan(length(selectedChewsStartTime), 0.3*fs);
selectedChewsJawInterp = nan(length(selectedChewsStartTime), 0.3*fs);
selectedChewsOutputFRInterp = nan(length(selectedChewsStartTime), 0.3*fs);
selectedChewsNormalizedSpikeInds = nan(length(selectedChewsStartTime)*20, 1);
tempInd = 1;
slideWindow = 900;
figure;


for i = 1:length(selectedChewsStartTime)
    
    selectedChewEphysStartInd = find(ephysTime >= selectedChewsStartTime(i), 1, 'first');
    selectedChewEphysStopInd = find(ephysTime <= selectedChewsStopTime(i), 1, 'last');
    
    outputEphysChewData = getVoltage(data.Data.Data(outputChannel, selectedChewEphysStartInd:selectedChewEphysStopInd));
    
    % get output cell firing rate
    [~, outputSpikeInds, ~] = crossdet(outputEphysChewData, outputThresh, 'thresholdDown');
    spikes = zeros(size(outputEphysChewData));
    spikes(outputSpikeInds) = 1;
    
    temp = 1:slideWindow/4:length(spikes);
    FR = nan(size(spikes));
    
    for j = 1:length(temp)
        ind = temp(j);
        startInd = ind - slideWindow/4;
        stopInd = ind + slideWindow/4;
        
        if startInd <= 0
            startInd = 1;
        end
        
        if stopInd > length(spikes)
            stopInd = length(spikes);
        end
        
        FR(j) = sum(spikes(startInd:stopInd)) / ((stopInd - startInd + 1)/fs);
    end
    FR = FR(~isnan(FR));
    selectedChewsOutputFR{i, 1} = FR;
    x = outputSpikeInds/length(outputEphysChewData);
    y = ones(size(x))*i;
    plot(x, y, '.', 'MarkerSize', 8); hold on;
    
    selectedChewsNormalizedSpikeInds(tempInd:tempInd+length(x)-1, 1) = x;
    tempInd = tempInd + length(x);
    
    % interp LFP and output cell FR onto the same time strech
    switch flag
        case 'LFP'
            x = selectedChewsLFPSmoothed{i, 1};
            selectedChewsLFPInterp(i, :) = interp1(1:length(x), x, ...
                linspace(1, length(x), size(selectedChewsLFPInterp, 2)));            
        case 'Jaw'
            x = selectedChewsJaw{i, 1};
            selectedChewsJawInterp(i, :) = interp1(1:length(x), x, ...
                linspace(1, length(x), size(selectedChewsJawInterp, 2)));          
    end
    
    selectedChewsOutputFRInterp(i, :) = interp1(1:length(FR), FR, ...
        linspace(1, length(FR), size(selectedChewsOutputFRInterp, 2)));
    
end

selectedChewsNormalizedSpikeInds = selectedChewsNormalizedSpikeInds(~isnan(selectedChewsNormalizedSpikeInds));


% save data to opto trial struct

switch flag
    case 'LFP'
        normalTrialsChewingStruct.selectedChewsLFPSmoothed = selectedChewsLFPSmoothed;
        normalTrialsChewingStruct.selectedChewsLFPInterp = selectedChewsLFPInterp;        
    case 'Jaw'
        normalTrialsChewingStruct.selectedChewsJawInterp = selectedChewsJawInterp;      
end

normalTrialsChewingStruct.selectedChewsStartTime = selectedChewsStartTime;
normalTrialsChewingStruct.selectedChewsStopTime = selectedChewsStopTime;
normalTrialsChewingStruct.selectedChewsOutputFR = selectedChewsOutputFR;
normalTrialsChewingStruct.selectedChewsOutputFRInterp = selectedChewsOutputFRInterp;
normalTrialsChewingStruct.selectedChewsNormalizedSpikeInds = selectedChewsNormalizedSpikeInds;
normalTrialsChewingStruct.outputThresh = outputThresh;

% calculate the vector strength for normal food trials
normalTrialChewsPhase = normalTrialsChewingStruct.selectedChewsNormalizedSpikeInds*2*pi;
x = sum(cos(normalTrialChewsPhase));
y = sum(sin(normalTrialChewsPhase));
normalVS = sqrt(x^2 + y^2) / length(normalTrialChewsPhase);

normalTrialsChewingStruct.VS = normalVS;


%% plot phase histogram!

figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
subplot(2, 2, 1);
x = linspace(0, 1, 0.3*fs);
shift = 0.6;
switch flag
    case 'LFP'
        chewsLFPPhase = optoTrialsChewingStruct.selectedChewsLFPInterp;
        ind = round(size(chewsLFPPhase, 2)*shift);
        if ind ~= 0
            shiftedChewsLFPPhase = [chewsLFPPhase(:, ind:end), chewsLFPPhase(:, 1:ind-1)];
        else
            shiftedChewsLFPPhase = chewsLFPPhase;
        end
        shadedErrorBar(x, mean(shiftedChewsLFPPhase), std(shiftedChewsLFPPhase), ...
            'lineProps', {'-', 'LineWidth', 2, 'Color', '#0072BD'});
        ylabel('rmsv');
        legend('LFP rmsv');
        
    case 'Jaw'
        chewsJawPhase = optoTrialsChewingStruct.selectedChewsJawInterp;
        ind = round(size(chewsJawPhase, 2)*shift);
        if ind ~= 0
            shiftedChewsJawPhase = [chewsJawPhase(:, ind:end), chewsJawPhase(:, 1:ind-1)];
        else
            shiftedChewsJawPhase = chewsJawPhase;
        end
        shadedErrorBar(x, mean(shiftedChewsJawPhase), std(shiftedChewsJawPhase), ...
            'lineProps', {'-', 'LineWidth', 2, 'Color', '#0072BD'});
        ylabel('pixel');
        legend('Lower Jaw Trace');
end
        
xlabel('phase');
title('Light On Trials');
box off;
axis tight;

subplot(2, 2, 2);
binNum = 20;
normalizedSpikeInds = optoTrialsChewingStruct.selectedChewsNormalizedSpikeInds;
if shift ~= 0
    normalizedSpikeInds = normalizedSpikeInds + shift;
    normalizedSpikeInds(normalizedSpikeInds > 1) = normalizedSpikeInds(normalizedSpikeInds > 1) - 1;
end
histogram(normalizedSpikeInds, binNum, ...
    'Normalization', 'probability', 'EdgeAlpha', 0.2);
box off;
axis tight;
% ytix = get(gca, 'YTick');
% set(gca, 'YTick',ytix, 'YTickLabel',ytix*100);
xlabel('phase');
ylabel('spike prob.');
ylim([0, 0.1]);
title('Light On Trials');

subplot(2, 2, 3);
x = linspace(0, 1, 0.3*fs);
switch flag
    case 'LFP'
        chewsLFPPhase = normalTrialsChewingStruct.selectedChewsLFPInterp;
        ind = round(size(chewsLFPPhase, 2)*shift);
        if ind ~= 0
            shiftedChewsLFPPhase = [chewsLFPPhase(:, ind:end), chewsLFPPhase(:, 1:ind-1)];
        else
            shiftedChewsLFPPhase = chewsLFPPhase;
        end
        shadedErrorBar(x, mean(shiftedChewsLFPPhase), std(shiftedChewsLFPPhase), ...
            'lineProps', {'-', 'LineWidth', 2, 'Color', [0.93, 0.596, 0.596]});
        ylabel('rmsv');
        legend('LFP rmsv');
        
    case 'Jaw'
        chewsJawPhase = normalTrialsChewingStruct.selectedChewsJawInterp;
        ind = round(size(chewsJawPhase, 2)*shift);
        if ind ~= 0
            shiftedChewsJawPhase = [chewsJawPhase(:, ind:end), chewsJawPhase(:, 1:ind-1)];
        else
            shiftedChewsJawPhase = chewsJawPhase;
        end
        shadedErrorBar(x, mean(shiftedChewsJawPhase), std(shiftedChewsJawPhase), ...
            'lineProps', {'-', 'LineWidth', 2, 'Color', [0.93, 0.596, 0.596]});
        ylabel('pixel');
        legend('Lower Jaw Trace');
end
xlabel('phase');
title('Light Off Trials');
box off
axis tight

subplot(2, 2, 4);
binNum = 20;
normalizedSpikeInds = normalTrialsChewingStruct.selectedChewsNormalizedSpikeInds;
if shift ~= 0
    normalizedSpikeInds = normalizedSpikeInds + shift;
    normalizedSpikeInds(normalizedSpikeInds > 1) = normalizedSpikeInds(normalizedSpikeInds > 1) - 1;
end

histogram(normalizedSpikeInds, binNum, 'Normalization', 'probability', ...
    'FaceColor',  [0.93, 0.596, 0.596], 'EdgeAlpha', 0.2);
box off;
axis tight;
% ytix = get(gca, 'YTick');
% set(gca, 'YTick',ytix, 'YTickLabel',ytix*100);
xlabel('phase');
ylabel('spike prob.');
ylim([0, 0.1]);
title('Light Off Trials');

%% save the data!

save(fullfile(sessionFolder, 'TrialsChewingStruct.mat'), 'normalTrialsChewingStruct', 'optoTrialsChewingStruct');


