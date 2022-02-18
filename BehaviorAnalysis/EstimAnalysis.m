%% initializations and data loading
session = '20220127_001';
ephysChannel = 29;
ephysThresh = -220; % unit: uV;

% settings
estimWindow = [-0.02, 0.20]; % unit: sec;
A.currentLevel = [0.2, 0.2, 0.3, 0.5, 0.3, 0.3]; % unit:mA; -> 20220127_001
S.currentLevel = [0.2, repmat(0.3, 1, 15)]; % unit:mA; -> 20220127_001
% A.currentLevel = repmat(0.8, 1, 10); % unit:mA; -> 20220128_003
% S.currentLevel = repmat(0.8, 1, 19); % unit:mA; -> 20220128_003
A.currentLevelToPlot = 0.5;
S.currentLevelToPlot = 0.3;

A.pulseRepeat = 10;
A.trainRepeat = 10;
S.pulseRepeat = 10;
S.trainRepeat = 10;

% initializations
rootFolder = 'Z:\Qianyun\DCN\';
gitFolder = 'D:\DCN_Project\Github\DCN\';
sessionFolder = fullfile(rootFolder, 'Data', session);

% load data
[spike, ~, sessionEphysInfo, getVoltage, data, channelNum_OpenEphys] = loadSessionData(session);
disp('Data loading finished!');

%% process the data - estim PSTH analysis

A.times = spike.keyboardTimes(strfind(convertCharsToStrings(spike.keyboardInput), 'A'));
S.times = spike.keyboardTimes(strfind(convertCharsToStrings(spike.keyboardInput), 'S'));

if length(A.times) ~= length(A.currentLevel) || length(S.times) ~= length(S.currentLevel)
    error('times does NOT match currentLevels!!! (A or S)');
end

% Start with A
currenLevelOccurInds = find(A.currentLevel == A.currentLevelToPlot);
currenLevelOccurCount = length(currenLevelOccurInds);
spikeRate = cell(currenLevelOccurCount*A.trainRepeat, 1);
for i = 1:length(currenLevelOccurInds)
    keyboardTime = A.times(currenLevelOccurInds(i));
    estimTimeInd = find(spike.EstimTimes >= keyboardTime, 1, 'first');
    estimStartTime = spike.EstimTimes(estimTimeInd);
    if (estimStartTime - keyboardTime > 0.1) 
        error('did NOT find the correct estim time!')
    end
    
    estimTrainStartTime(1) = estimStartTime + estimWindow(1);
    estimTrainEndTime(1) = spike.EstimTimes(estimTimeInd + A.pulseRepeat-1) + estimWindow(2);
    ephysTrainStartInd(1) = find(sessionEphysInfo.convertedEphysTimestamps >= estimTrainStartTime(1), 1, 'first');
    ephysTrainEndInd(1) = find(sessionEphysInfo.convertedEphysTimestamps <= estimTrainEndTime(1), 1, 'last');
    
    for j = 1:A.trainRepeat
        ephysTrainData = getVoltage(data.Data.Data(channelNum_OpenEphys(ephysChannel), ephysTrainStartInd(j):ephysTrainEndInd(j)));
        [~, spikeInds, ~] = crossdet(ephysTrainData, ephysThresh, 'thresholdDown');
        spikeTimes = (spikeInds/sessionEphysInfo.fs) + estimWindow(1);
        [spikeRate{(i-1)*A.trainRepeat+j}, spikeRateTimes] = getFiringRate(spikeTimes, 'tLims', [estimWindow(1), estimTrainEndTime(j) - estimTrainStartTime(j) + estimWindow(1)]);
        
        if j == 9
            figure('Color', 'white','WindowState','maximized'); clf
            x = linspace(estimWindow(1), estimTrainEndTime(j) - estimTrainStartTime(j) + estimWindow(1), length(ephysTrainData));
            % x = sessionEphysInfo.convertedEphysTimestamps(ephysTrainStartInd(j) : ephysTrainEndInd(j));
            plot(x, ephysTrainData); hold on;
            % startLine = sessionEphysInfo.convertedEphysTimestamps(ephysTrainStartInd(j)) - estimWindow(1);
            % endLine = sessionEphysInfo.convertedEphysTimestamps(ephysTrainEndInd(j)) - estimWindow(2);
            plot([0, 0], [min(ephysTrainData), max(ephysTrainData)], '-c', 'LineWidth', 2);
            plot([0.1, 0.1], [min(ephysTrainData), max(ephysTrainData)], '-r', 'LineWidth', 2);
            % plot(sessionEphysInfo.convertedEphysTimestamps(spikeInds+ephysTrainStartInd(j)), repmat(-220, 1, length(spikeInds)), '.', 'MarkerSize', 12);
            plot(spikeRateTimes, spikeRate{(i-1)*A.trainRepeat + j});
            ylabel('uV');
            xlabel('time');
            axis tight;
        end
        
        if j < 10
            estimTrainStartTime(j+1) = spike.EstimTimes(estimTimeInd + j*A.pulseRepeat) + estimWindow(1);
            estimTrainEndTime(j+1) = spike.EstimTimes(estimTimeInd + j*A.pulseRepeat + A.pulseRepeat-1) + estimWindow(2);
            ephysTrainStartInd(j+1) = find(sessionEphysInfo.convertedEphysTimestamps >= estimTrainStartTime(j+1), 1, 'first');
            ephysTrainEndInd(j+1) = find(sessionEphysInfo.convertedEphysTimestamps <= estimTrainEndTime(j+1), 1, 'last');
        end
        
    end
    
    
end

% plot PSTH!
figure('Color', 'white','WindowState','maximized'); clf
tempLength = length(spikeRate{1});
for i = 1:length(spikeRate)
    trialLength = length(spikeRate{i});
    tempLength = min(trialLength, tempLength);
end
minCommonLength = tempLength;

spikeRateArray = nan(length(spikeRate), minCommonLength);
for i = 1:length(spikeRate)
    spikeRateArray(i, :) = spikeRate{i}(1:minCommonLength);
end

x = linspace(-0.02, 0.29, size(spikeRateArray, 2));
shadedErrorBar(x, nanmean(spikeRateArray), nanstd(spikeRateArray), ...
     'lineprops', {'linewidth', 3, 'color', '#EDB120'}, 'patchSaturation', .1);
xlabel('time (s)');
ylabel('spike rate (spikes/s)');

hold on;
ylim([-30, 160]);
xlim([-0.02, 0.29]);
plot([-0.02, 0.29], [-20, -20], '-k');
for i = 1:A.pulseRepeat
    plot([0.01*(i-1), 0.01*(i-1)], [-20, -17], '-k', 'HandleVisibility', 'off');
end

plot([0, 0], [-30, 160], '--m');
plot([0.1, 0.1], [-30, 160], '--r');

shadedErrorBar(x, nanmean(OutputCellSpikeRateArray), nanstd(OutputCellSpikeRateArray), ...
     'lineprops', {'linewidth', 3, 'color', '#0072BD'}, 'patchSaturation', .1);


legend('cwc avg spike rate (0.5mA estim)', 'estim pulses (100Hz, 100ms duration)', 'pulse start', 'pulse end', 'output cell avg spike rate (0.8mA estim)');

