% settings
session = '20220816_003';
ephysChannel = 30;
Thresh = -100;
toneMapCues = 'KLHUOPY';
fs = 30000;

expType = 'awake';

%% initializations & loading the data
rootFolder = 'Z:\Qianyun\DCN\';
gitFolder = 'D:\DCN_Project\Github\DCN\';
sessionFolder = fullfile(rootFolder, 'Data', session);

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
    
    % determine target channel
    if contains(sessionEphysInfo.ephysFolder, '101')
        ephysChannel = channelNum_OpenEphys(ephysChannel);
    end
    
    disp('finish loading ephys data');   
end

% load responses.mat and optoResponses.mat
if ~exist(fullfile(sessionFolder, 'responses.mat'), 'file')
    error('responses.mat NOT exist!');
else
    load(fullfile(sessionFolder, 'responses.mat'), 'responses');
end

if ~exist(fullfile(sessionFolder, 'optoResponses.mat'), 'file')
    error('optoResponses.mat NOT exist!');
else
    load(fullfile(sessionFolder, 'optoResponses.mat'), 'optoResponses');
end


%% get the opto on off times (only for spontaneous firing) 

optoOnOffTimes = getOptoSpFROnOffTimes(spike, 'expType', expType);
inds = find(optoOnOffTimes(:, 2) - optoOnOffTimes(:, 1) < 20);
optoOnOffTimes(inds, :) = [];


%% calculate the firing rate and plot!

slideWindow = 45000;
outputFR = cell(length(optoOnOffTimes), 1);
FRTimes = cell(length(optoOnOffTimes), 1);
outputFRInterped = nan(301, length(optoOnOffTimes));
figure('Color', 'white','WindowState','maximized'); clf;
colors = cool(length(optoOnOffTimes));

for i = 1:length(optoOnOffTimes)
    
    ephysStartTime = optoOnOffTimes(i,1) - 5;
    ephysStopTime = optoOnOffTimes(i,2) + 5;
    ephysStartInd = find(ephysTime >= ephysStartTime, 1, 'first');
    ephysStopInd = find(ephysTime <= ephysStopTime, 1, 'last');
    
    ephysChunkData = getVoltage(data.Data.Data(ephysChannel, ephysStartInd:ephysStopInd));
    ephysChunkTime = ephysTime(ephysStartInd:ephysStopInd);
    
    [~, spikeInds, ~] = crossdet(ephysChunkData, Thresh, 'thresholdDown');

    
    spikes = zeros(size(ephysChunkData));
    spikes(spikeInds) = 1;
    
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
    
    FR = FR(~isnan(FR));
    outputFR{i, 1} = FR;
    Time = linspace(-5, ephysStopTime - ephysStartTime - 5, length(spikes));
    Time = Time(temp);
    FRTimes{i, 1} = Time;
    
    % interperate the FR data
    inds1 = find(Time >= -5 & Time < 0);
    inds2 = find(Time >= 0 & Time <= 20);
    inds3 = find(Time >= Time(end) - 5);
    
    outputFRInterped(i, 1:50) = interp1(1:length(inds1), FR(inds1), linspace(1, length(inds1), 50));
    outputFRInterped(i, 51:251) = interp1(1:length(inds2), FR(inds2), linspace(1, length(inds2), 201));
    outputFRInterped(i, 252:301) = interp1(1:length(inds3), FR(inds3), linspace(1, length(inds3), 50));
    
    plot(linspace(-5, 25, 301), outputFRInterped(i, :), 'Color', colors(i, :)); hold on;
end

% if having different laser power levels
optoLevel = [1, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4];
laserPower = {'3%', '10%', '2.0mW', '3.0mW',};
colors = cool(max(optoLevel));

% % if having only one laswer power levels
% optoLevel = ones(1, length(optoOnOffTimes));
% laserPower = {'5%'};
figure('Color', 'white','WindowState','maximized'); clf;


for i = 3:max(optoLevel)
    
    inds = find(optoLevel == i); 
    levelFR = outputFRInterped(inds, :);
    
    if length(inds) == 1
        plot(linspace(-5, 25, 301), levelFR, '-', 'LineWidth', 2, 'Color', colors(i, :));
        hold on;
    else
        shadedErrorBar(linspace(-5, 25, 301), mean(levelFR), std(levelFR),...
            'lineProps', {'-', 'LineWidth', 2, 'Color', colors(i, :)}); hold on
    end
    
%     shadedErrorBar(linspace(-5, 25, 301), mean(levelFR), std(levelFR),...
%         'lineProps', {'-', 'LineWidth', 2, 'Color', [0.7961, 0.5765, 0.9294]}); hold on  
%     
    if i == 1
        maxFR = max(mean(levelFR));
    else
        temp = max(mean(levelFR));
        if temp > maxFR
            maxFR = temp;
        end
    end
    
end

for i = 3:max(optoLevel)
    
%     text(21.5, maxFR+10-(i-1)*0.8, ['laser power = ', laserPower{i}],...
%         'Color', [0.7961, 0.5765, 0.9294], 'FontSize', 10)
    
    text(21.5, maxFR+12-(i-2-1)*0.9, ['laser power = ', laserPower{i}],...
        'Color', colors(i, :), 'FontSize', 10);
    box off;
        
end


plot([0, 20], [maxFR+12, maxFR+12], '-', 'LineWidth', 2, 'Color',  [0.3, 0.74, 0.933]);
text(9, maxFR+13, 'OPTO ON', 'Color', [0.3, 0.74, 0.933], 'FontWeight', 'bold');
ylim([0, maxFR+17]);
ylabel('Spike rate (spks/s)');
xlabel('Time(s)');
title([session, ' Ch', num2str(ephysChannel), ' Cartwheel Cell, Opto Effect on Spon. FR'], 'Interpreter', 'none');

