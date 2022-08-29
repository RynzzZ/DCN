%% session 20220609_000, opto effect, ch10(output cell), ch4(cwc)
session = '20220609_000';
outputChannel = 8;
cwcChannel = [];
timeWindow = [-1, 20];
cwcThresh = -200;
outputThresh = -200;
toneMapCues = 'KLHUOP';
fs = 30000;

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
    ephysTime = [1:size(data.Data.Data, 2)]/sessionEphysInfo.fs;
    disp('finish loading ephys data');
end

%% quality check
figure;
subplot(2, 1, 1)
startTime = 2089;
endTime = 2116;


startInd = find(ephysTime >= startTime, 1, 'first');
endInd = find(ephysTime <= endTime, 1, 'last');

if ~isempty(cwcChannel)    
    x = ephysTime(startInd:endInd);
    y = getVoltage(data.Data.Data(channelNum_OpenEphys(cwcChannel), startInd:endInd));
    
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
    startTime = 1681;
    endTime = 1708;
    
    
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%processing data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% locate opto train start
optoTrialNum = 50;
optoFoodTrialStartTime = nan(optoTrialNum, 1);
normalFoodTrialStartTime = nan(optoTrialNum, 1);
outputEphysData = cell(optoTrialNum, 2);
cwcEphysData = cell(optoTrialNum, 2);
cwcSpikeInds = cell(optoTrialNum, 1);
outputSpkRates = cell(optoTrialNum, 2);
outputSpkTimes = cell(optoTrialNum, 2);

normalFoodColumn = 1;
optoFoodColum = 2;

foodTrialCount = 0;
optoFoodTrialCount = 0;

slideWindow = 7500;

for i = 1:length(spike.keyboardInput)-2
    if strcmp(spike.keyboardInput(i), 'F')
        if strcmp(spike.keyboardInput(i+1), 'G')
            foodTriggerTime = spike.keyboardTimes(i);
            lightTriggerTime = spike.keyboardTimes(i+1);
            
            if lightTriggerTime - foodTriggerTime > 3 || sum(spike.optoTrainTimes >= lightTriggerTime & spike.optoTrainTimes <= lightTriggerTime + 10)<100
                continue
            else
                optoFoodTrialCount = optoFoodTrialCount + 1;
                disp(['opto food trial: ', num2str(optoFoodTrialCount)]);
                optoFoodTrialStartTime(optoFoodTrialCount) = foodTriggerTime;
                
                % determine trial start and stop time
                ephysStartTime = foodTriggerTime + timeWindow(1);
                if strcmp(spike.keyboardInput(i+2), 'G')
                    ephysStopTime = spike.keyboardTimes(i+2);
                else
                    ephysStopTime = lightTriggerTime + timeWindow(2);
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
                
                
                % quality check
                % if mod(i, 5) == 0
                %                         figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
                %                         x = linspace(-5, ephysStopTime - ephysStartTime, ephysStopInd - ephysStartInd + 1);
                %                         plot(x, cwcEphysData{count, 1}); axis tight; box off; hold on;
                %                         y = ones(size(cwcSpikeInds{count, 1}))*cwcThresh;
                %                         scatter(x(cwcSpikeInds{count, 1}), y, '.', 'r');
                
                %                         figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
                %                         subplot(2, 1, 1);
                %                         x = linspace(-5, ephysStopTime - ephysStartTime, ephysStopInd - ephysStartInd + 1);
                %                         plot(x, outputEphysData{count, 1}); axis tight; box off; hold on;
                %                         subplot(2, 1, 2);
                %                         plot(outputSpkTimes{optoFoodTrialCount, 1}, outputSpkRates{optoFoodTrialCount, 1}, '-r', 'LineWidth', 2);
                %                         axis tight; box off;
                %                         ylim([0, 60]);
                % end
                
            end
            
        elseif spike.keyboardTimes(i+1) - spike.keyboardTimes(i) > 20
            foodTrialCount = foodTrialCount + 1;
            disp(['normal food trial: ', num2str(foodTrialCount)]);
            
           
            % determine trial start and stop time
            foodTriggerTime = spike.keyboardTimes(i);
            normalFoodTrialStartTime(foodTrialCount) = foodTriggerTime;
            ephysStartTime = foodTriggerTime + timeWindow(1);
            ephysStopTime = foodTriggerTime + timeWindow(2);
            
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
           
            % output cell
            % voltage data
            outputEphysData{foodTrialCount, normalFoodColumn} = getVoltage(data.Data.Data(outputChannel, ephysStartInd:ephysStopInd));
            
            % get output cell firing rate
            [~, outputSpikeInds, ~] = crossdet(outputEphysData{foodTrialCount, normalFoodColumn}, outputThresh, 'thresholdDown');
            spikes = zeros(size(outputEphysData{foodTrialCount, normalFoodColumn}));
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
            
            outputSpkRates{foodTrialCount, normalFoodColumn} = FR(~isnan(FR));
            Time = linspace(timeWindow(1), ephysStopTime - ephysStartTime + timeWindow(1), length(spikes));
            outputSpkTimes{foodTrialCount, normalFoodColumn} = Time(temp);
            
        
        end
    end
end


%% make plots!


figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
rows = 6;
cols = 1;
normalFoodTrial = 1;
optoFoodTrial = 1;
plotInd = 1;

for i = 1:36
    subplot(rows, cols, plotInd);
    if mod(i, 2) == 1
        plot(outputSpkTimes{normalFoodTrial, normalFoodColumn}, outputSpkRates{normalFoodTrial, normalFoodColumn}, '-k', 'LineWidth', 2);
        axis tight;
        normalFoodTrial = normalFoodTrial + 1;   
    else
        plot(outputSpkTimes{optoFoodTrial, optoFoodColum}, outputSpkRates{optoFoodTrial, optoFoodColum}, '-b', 'LineWidth', 2);
        axis tight;
        optoFoodTrial = optoFoodTrial + 1;
    end
    
    plotInd = plotInd + 1;
    if plotInd > rows && i <36
        plotInd = 1;
        figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
    end

end




figure('Color', 'white', 'position', get(0,'ScreenSize')); clf; 
violinplot(dataMatrix, names); 
box off;
ylabel('normalized corss corr val')








