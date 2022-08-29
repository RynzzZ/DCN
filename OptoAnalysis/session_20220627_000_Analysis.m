%% session 20220627_000, opto effect, ch12(potential output cell), ch22(LFP)
session = '20220627_000';
outputChannel = 12;
LFPChannel = 22;
cwcChannel = [];
timeWindow = [-1, 20];
cwcThresh = -200;
outputThresh = -500;
toneMapCues = 'KLHUOP';
fs = 30000;

% initializations
rootFolder = 'Z:\Qianyun\DCN\';
gitFolder = 'D:\DCN_Project\Github\DCN\';
sessionFolder = fullfile(rootFolder, 'Data', session);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%processing data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load spikeAnalyzed.mat
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
    disp('finish loading ephys data');
end

if contains(sessionEphysInfo.ephysFolder, '101')
    if any(outputChannel)
        outputChannel = channelNum_OpenEphys(outputChannel);
    end
    
    if any(cwcChannel)
        cwcChannel = channelNum_OpenEphys(cwcChannel);
    end
    
    if any(LFPChannel)
        LFPChannel = channelNum_OpenEphys(LFPChannel);
    end
    
end

%% quality check
figure;
subplot(2, 1, 1)
startTime = 120;
endTime = 135;


startInd = find(ephysTime >= startTime, 1, 'first');
endInd = find(ephysTime <= endTime, 1, 'last');

if ~isempty(cwcChannel)    
    x = ephysTime(startInd:endInd);
    y = getVoltage(data.Data.Data(cwcChannel, startInd:endInd));
    
    plot(x, y); hold on;
    plot([x(1), x(end)], [cwcThresh, cwcThresh], '-r');
    axis tight
end

if ~isempty(outputChannel)
    
    subplot(2, 1, 1)
    x = ephysTime(startInd:endInd);
    y = getVoltage(data.Data.Data(outputChannel, startInd:endInd));
    
    plot(x, y); hold on;
    plot([x(1), x(end)], [outputThresh, outputThresh], '-r');
    axis tight
    title('example opto food trial')
    
    subplot(2, 1, 2)
    startTime = 40;
    endTime = 55;
    
    
    startInd = find(ephysTime >= startTime, 1, 'first');
    endInd = find(ephysTime <= endTime, 1, 'last');
    
    x = ephysTime(startInd:endInd);
    y = getVoltage(data.Data.Data(outputChannel, startInd:endInd));
    
    plot(x, y); hold on;
    plot([x(1), x(end)], [outputThresh, outputThresh], '-r');
    axis tight
    title('example normal food trial')
end

%% process the data

dataStartTime = 1;
dataStopTime = 896;
dataStartInd = find(ephysTime >= dataStartTime, 1, 'first');
dataStopInd = find(ephysTime <= dataStopTime, 1, 'last');
ephysData = getVoltage(data.Data.Data(outputChannel, dataStartInd:dataStopInd));
LFPData = getVoltage(data.Data.Data(LFPChannel, dataStartInd:dataStopInd));
LFPData = LFPData - mean(LFPData); % dc remove
LFPrmsv = sqrt(movmean(LFPData.^2, 100)); % turn into rmsv

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
            
            if ephysStopTime >= dataStopTime
                disp('break');
                break;
            end
            
            
            % get ephys time inds
            ephysStartInd = find(ephysTime >= ephysStartTime, 1, 'first');
            ephysStopInd = find(ephysTime <= ephysStopTime, 1, 'last');
            
            % cwc cell
            if ~isempty(cwcChannel)
                % voltage data
                cwcEphysData{count, 1} = getVoltage(data.Data.Data(channelNum_OpenEphys(cwcChannel), ephysStartInd:ephysStopInd));
                cwcEphysData{count, 2} = ephysTime(ephysStartInd:ephysStopInd) - ephysTime(ephysStartInd) - 5;
                
                % get cwc spike times
                [~, cwcSpikeInds{count, 1}, ~] = crossdet(cwcEphysData{count, 1}, cwcThresh, 'thresholdDown');
            end
            
            % LFP
            LFPtemp = LFPrmsv(ephysStartInd:ephysStopInd);
            LFPsmoothed = smooth(LFPtemp, 0.01);
            LFPEphysData{optoFoodTrialCount, optoFoodColum} = LFPsmoothed;
            ephysTimeData{optoFoodTrialCount, optoFoodColum} = ephysTime(ephysStartInd:ephysStopInd);
            
            % output cell
            % voltage data
            ephysData = getVoltage(data.Data.Data(outputChannel, ephysStartInd:ephysStopInd));
            ephysData = removeOptoArtifacts(ephysData, [-400, -100], []);
            outputEphysData{optoFoodTrialCount, optoFoodColum} = ephysData;
            
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
            
            % output cell spk rate
            y5 = outputSpkRates{optoFoodTrialCount, optoFoodColum};
            xx = linspace(x(1), x(end), length(y5));
            y5 = (y5 - min(y5))/range(y5)*2;
            r = ksr(1:length(y5), y5, 1.33, length(y5));
            plot(xx, r.f)
            
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
            
            if ephysStopTime >= dataStopTime
                disp('break');
                break;
            end
            
            % get ephys time inds
            ephysStartInd = find(ephysTime >= ephysStartTime, 1, 'first');
            ephysStopInd = find(ephysTime <= ephysStopTime, 1, 'last');
            
            % cwc cell
            if ~isempty(cwcChannel)
                % voltage data
                cwcEphysData{count, 1} = getVoltage(data.Data.Data(channelNum_OpenEphys(cwcChannel), ephysStartInd:ephysStopInd));
                cwcEphysData{count, 2} = ephysTime(ephysStartInd:ephysStopInd) - ephysTime(ephysStartInd) - 5;
                
                % get cwc spike times
                [~, cwcSpikeInds{count, 1}, ~] = crossdet(cwcEphysData{count, 1}, cwcThresh, 'thresholdDown');
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
            
            % output cell spk rate
            y5 = outputSpkRates{normalFoodTrialCount, normalFoodColumn};
            xx = linspace(x(1), x(end), length(y5));
            y5 = (y5 - min(y5))/range(y5)*2;
            r = ksr(1:length(y5), y5, 1.33, length(y5));
            plot(xx, r.f)
            
            
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






















