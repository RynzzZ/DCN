%% session 20220329_000, opto effect, ch10(output cell), ch4(cwc)
session = '20220329_000';
outputChannel = 10;
cwcChannel = 4;
timeWindow = [-5, 5];
cwcThresh = -100;
outputThresh = -100;
toneMapCues = 'JKLHUOPY';
fs = 30000;

% initializations
rootFolder = 'Z:\Qianyun\DCN\';
gitFolder = 'D:\DCN_Project\Github\DCN\';
sessionFolder = fullfile(rootFolder, 'Data', session);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%processing data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load spikeAnalyzed.mat
if ~exist(fullfile(sessionFolder, 'spikeAnalyzed.mat'), 'file')
    error('spikeAnalyzed.mat NOT exist! Run analyzeSession.m first!');
else
    disp('loading spikeAnalyzed.mat...')
    load(fullfile(sessionFolder, 'spikeAnalyzed.mat'), 'spike');
    
    % get sampling frequency for audio sig channel and mic channel in spike
    s.audioFs = round(length(spike.audioSignalTimes)/(spike.audioSignalTimes(end) - spike.audioSignalTimes(1)));
    s.micFs = round(length(spike.micSignalTimes)/(spike.micSignalTimes(end) - spike.micSignalTimes(1)));
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
end

% quality check
figure;
subplot(2, 1, 1)
startTime = 915;
endTime = 986;

startInd = find(ephysTime >= startTime, 1, 'first');
endInd = find(ephysTime <= endTime, 1, 'last');

x = ephysTime(startInd:endInd);
y = getVoltage(data.Data.Data(channelNum_OpenEphys(cwcChannel), startInd:endInd));

plot(x, y); hold on;
plot([x(1), x(end)], [cwcThresh, cwcThresh], '-r');
axis tight


subplot(2, 1, 2)
x = ephysTime(startInd:endInd);
y = getVoltage(data.Data.Data(channelNum_OpenEphys(outputChannel), startInd:endInd));

plot(x, y); hold on;
plot([x(1), x(end)], [outputThresh, outputThresh], '-r');
axis tight

%% process the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%processing data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% locate opto train start
optoTrainNum = 14;
selectedOptoStartTime = nan(optoTrainNum, 1);
outputEphysData = cell(optoTrainNum, 2);
cwcEphysData = cell(optoTrainNum, 2);
cwcSpikeInds = cell(optoTrainNum, 1);
outputSpkRates = cell(optoTrainNum, 1);
outputSpkTimes = cell(optoTrainNum, 1);
count = 0;
for i = 1:length(spike.keyboardInput)-1
    if strcmp(spike.keyboardInput(i), 'G')
        switch (spike.keyboardInput(i+1))
            case 'G'
                startTime = spike.keyboardTimes(i);
                stopTime = spike.keyboardTimes(i+1);
                
                if (sum(spike.optoTrainTimes >= startTime & spike.optoTrainTimes <= stopTime)<100)...
                        || (stopTime - startTime < 20) ...
                        || any(strfind(toneMapCues, spike.keyboardInput(i-1)))
                    continue
                else
                    count = count + 1;
                    disp(count);
                    selectedOptoStartTime(count) = startTime;
                    
                    % get ephys data
                    ephysStartTime = startTime + timeWindow(1);
                    ephysStopTime = stopTime + timeWindow(2);
                    
                    ephysStartInd = find(ephysTime >= ephysStartTime, 1, 'first');
                    ephysStopInd = find(ephysTime <= ephysStopTime, 1, 'last');
                    
                    outputEphysData{count, 1} = getVoltage(data.Data.Data(channelNum_OpenEphys(outputChannel), ephysStartInd:ephysStopInd));
                    cwcEphysData{count, 1} = getVoltage(data.Data.Data(channelNum_OpenEphys(cwcChannel), ephysStartInd:ephysStopInd));
                    cwcEphysData{count, 2} = ephysTime(ephysStartInd:ephysStopInd) - ephysTime(ephysStartInd) - 5;
                    
                    % get cwc spike times
                    [~, cwcSpikeInds{count, 1}, ~] = crossdet(cwcEphysData{count, 1}, cwcThresh, 'thresholdDown');
                    
                    % get output cell firing rate
                    [~, outputSpikeInds, ~] = crossdet(outputEphysData{count, 1}, outputThresh, 'thresholdDown');
                    spikes = zeros(size(outputEphysData{count, 1}));
                    spikes(outputSpikeInds) = 1;
                    
                    slideWindow = 15000;
                    temp = 1:7500:length(spikes);
                    FR = nan(size(spikes));
                    
                    parfor j = 1:length(temp)
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
                    
                    outputSpkRates{count, 1} = FR(~isnan(FR));
                    Time = linspace(-5, ephysStopTime - ephysStartTime, length(spikes));
                    outputSpkTimes{count, 1} = Time(temp);
                    
         
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
%                         plot(outputSpkTimes{count, 1}, outputSpkRates{count, 1}, '-r', 'LineWidth', 2);
%                         axis tight; box off;
%                         ylim([0, 60]);
                    % end
                    
                end
        end
    end
    
end

% interp the output cell firing rate
outputFR = zeros(optoTrainNum, 112);
inds1 = [];
inds2 = [];
inds3 = [];
for i = 1:size(outputSpkRates, 1)

% 5s before the opto light turn on
ind1 = find(outputSpkTimes{i, 1} < 0);
ind2 = find(outputSpkTimes{i, 1} > 0 & outputSpkTimes{i, 1} < 20);
ind3 = find(outputSpkTimes{i, 1} > outputSpkTimes{i, 1}(end) - 5);

inds1(i) = length(ind1);
inds2(i) = length(ind2);
inds3(i) = length(ind3);

outputFR(i, 1:20) = interp1(outputSpkTimes{i, 1}(ind1), outputSpkRates{i, 1}(ind1), linspace(-5, 0, 20), 'linear','extrap');
outputFR(i, 21:92) = interp1(outputSpkTimes{i, 1}(ind2), outputSpkRates{i, 1}(ind2), linspace(0, 20, 72), 'linear','extrap');
outputFR(i, 93:112) = interp1(outputSpkTimes{i, 1}(ind3), outputSpkRates{i, 1}(ind3), ...
    linspace(outputSpkTimes{i, 1}(end)-5, outputSpkTimes{i, 1}(end), 20), 'linear','extrap');
end

%% make the plot!!!
outputCell_vertical_offset = 40;
figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
x = linspace(-5, 25, length(outputFR));
shadedErrorBar(x, mean(outputFR) + outputCell_vertical_offset, std(outputFR), 'lineProps', {'-', 'Color', '#7E2F8E', 'LineWidth', 2});
ax = gca;
ax.YAxis.Visible = 'off';

hold on;
for i = 1:length(cwcSpikeInds)
    x = cwcEphysData{i, 2}(cwcSpikeInds{i, 1});
    x(x > 20 & x < x(end) - 5) = [];
    x(x>x(end) - 5) = x(end) - x(x>x(end)-5) + 20;
    y = ones(size(x))*i*1.5;
    plot(x, y, '.', 'Color', '#A2142F', 'MarkerSize', 10);
end
ylim([-5, 80]); hold on;
xlim([-5.5, 25]);

% light ON and OFF line
plot([0, 0], [-5, 80], '--', 'Color', '#0072BD');
text(-0.8, 83, 'Light ON', 'Color', '#0072BD', 'FontWeight', 'Bold')
plot([20, 20], [-5, 80], '--', 'Color', '#0072BD');
text(19.2, 83, 'Light OFF', 'Color', '#0072BD', 'FontWeight', 'Bold');

% scale bar for FR
plot([-5, -5], [0, 30] + outputCell_vertical_offset, '-k');

plot([-5.2, -5], [0, 0] + outputCell_vertical_offset, '-k');
text(-6, 0 + outputCell_vertical_offset, '0');
plot([-5.2, -5], [10, 10] + outputCell_vertical_offset, '-k');
text(-6, 10 + outputCell_vertical_offset, '10');
plot([-5.2, -5], [20, 20] + outputCell_vertical_offset, '-k');
text(-6, 20 + outputCell_vertical_offset, '20');
plot([-5.2, -5], [30, 30] + outputCell_vertical_offset, '-k');
text(-6, 30 + outputCell_vertical_offset, '30');
% plot([-5.2, -5], [80, 80], '-k');
% text(-6, 80, '40');
% plot([-5.2, -5], [90, 90], '-k');
% text(-6, 90, '50');

h = text(-6.7, 10 + outputCell_vertical_offset, 'spk rate (sp/s)');
set(h,'Rotation',90);

legend('Ch10 output cell spk rate', 'Ch4 cwc spikes raster');
xlabel('time (sec)')


