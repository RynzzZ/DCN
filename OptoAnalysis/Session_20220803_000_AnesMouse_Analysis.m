% settings
session = '20220803_000';
outputChannel = 29;
outputThresh = -75;
toneMapCues = 'KLHUOPY';
fs = 30000;

%% initializations
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
        outputChannel = channelNum_OpenEphys(outputChannel);
        % cwcChannel = channelNum_OpenEphys(cwcChannel);
        % LFPChannel = channelNum_OpenEphys(LFPChannel);
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

%% plot the BBN RLF data (opto off vs. opto on at different light intensity level)

optoLevel = [0, 0, 1, 0, 1, 0, 1, 2, 2, 2, 3, 3, 3, 4, 0, 4, 0, 4, 3, 3, 2, 2, 1, 0, 1]; 
% 0 = opto off; 1 = 11% 1.5mW; 2 = 20%, 3.0mW; 3 = 30%, 5.0mW; 4 = 42%, 7mW. 

% opto off
inds = find(optoLevel == 0);
inds = inds(3:4);

BBNFR = mean(responses.B(inds, :));
BBNSTD = std(responses.B(inds, :));
[sortedLoudnessBBN, sortInds] = sort(loudness.RLFBBN);
sortedOptoOffBBNFR = BBNFR(sortInds);
sortedOptoOffBBNSTD = BBNSTD(sortInds);

% apply a 3-point trigular filter
filteredBBNFR = nan(size(sortedOptoOffBBNFR));
for i = 1:length(sortedOptoOffBBNFR)
    if i == 1
        filteredBBNFR(i) = (2*sortedOptoOffBBNFR(i) + sortedOptoOffBBNFR(i+1))/3;
    elseif i == length(sortedOptoOffBBNFR)
        filteredOptoOffBBNFR(i) = (2*sortedOptoOffBBNFR(i) + sortedOptoOffBBNFR(i-1))/3;
    else
        filteredOptoOffBBNFR(i) = (2*sortedOptoOffBBNFR(i) + sortedOptoOffBBNFR(i-1) + sortedOptoOffBBNFR(i+1))/4;
    end
end

% opto on
sortedOptoBBNFR = cell(max(optoLevel), 1);
sortedOptoBBNFRSTD = cell(max(optoLevel), 1);
filteredOptoBBNFR = cell(max(optoLevel), 1);
for i = 1:3
    inds = find(optoLevel == i);
    inds = inds(1:3);
    BBNFR = mean(optoResponses.B(inds, :));
    BBNSTD = std(optoResponses.B(inds, :));
    [sortedLoudnessBBN, sortInds] = sort(loudness.RLFBBN);
    sortedBBNFR = BBNFR(sortInds);
    sortedOptoBBNFR{i, 1} = BBNFR(sortInds);
    sortedOptoBBNFRSTD{i, 1} = BBNSTD(sortInds);
    
    % apply a 3-point trigular filter
    for j = 1:length(sortedOptoBBNFR{i, 1})
        if j == 1
            filteredBBNFR(j) = (2*sortedBBNFR(j) + sortedBBNFR(j+1))/3;
        elseif j == length(sortedOptoBBNFR{i, 1})
            filteredBBNFR(j) = (2*sortedBBNFR(j) + sortedBBNFR(j-1))/3;
        else
            filteredBBNFR(j) = (2*sortedBBNFR(j) + sortedBBNFR(j-1) + sortedBBNFR(j+1))/4;
        end
    end
    
    filteredOptoBBNFR{i, 1} = filteredBBNFR;
    
end

%% Plot - BBN RLF opto off + opto on at different light intensity level

% PLOT!!!
% Plot 2 - BBN RLF Plots, with shaded error bar. OPTO ON vs. OPTO OFF
figure('Color', 'white','WindowState','maximized'); clf;
shadedErrorBar(sortedLoudnessBBN, sortedOptoOffBBNFR, sortedOptoOffBBNSTD,...
    'lineProps', {'-', 'LineWidth', 2, 'Color', [0.65 0.65 0.65]});
hold on; box off; axis tight;

color = cool(max(optoLevel));
for i = 1:3
    shadedErrorBar(sortedLoudnessBBN, sortedOptoBBNFR{i, 1}, sortedOptoBBNFRSTD{i, 1},...
        'lineProps', {'-', 'LineWidth', 2, 'Color', color(i, :)}); hold on
end
plot([sortedLoudnessBBN(1), sortedLoudnessBBN(end)], [responses.baselineFR, responses.baselineFR],...
    '--', 'Color', [0.8 0.8 0.8], 'LineWidth', 1.5);

xlim([sortedLoudnessBBN(1), sortedLoudnessBBN(end)]);
xlabel('loudness (dB SPL)');
ylabel('rate (spk/s)');
legend('Opto Off', 'Opto On (1.5mW)', 'Opto On (3.0mW)', 'Opto On (5.0mW)',  'baseline');
title([session, ' Ch', num2str(ephysChannel), ' BBN RLF'], 'Interpreter', 'none');


%% Plot - comparing BBN RLF with the same opto level but at different time of the session
figure('Color', 'white','WindowState','maximized'); clf;
rows = 2; cols = 2;

for i = 1:3
    % plot opto off traces
    subplot(rows, cols, i)
    shadedErrorBar(sortedLoudnessBBN, sortedOptoOffBBNFR, sortedOptoOffBBNSTD, 'lineProps', {'-', 'LineWidth', 2, 'Color', [0.65, 0.65, 0.65]});
    hold on; box off; axis tight;
    
    % (early) opto on trace
    inds = find(optoLevel == i);
    inds = inds(1:3);
    BBNFR = mean(optoResponses.B(inds, :));
    BBNSTD = std(optoResponses.B(inds, :));
    [sortedLoudnessBBN, sortInds] = sort(loudness.RLFBBN);
    sortedBBNFR = BBNFR(sortInds);
    sortedBBNSTD = BBNSTD(sortInds);
    
    % plot early opto on traces
    shadedErrorBar(sortedLoudnessBBN, sortedBBNFR, sortedBBNSTD,...
        'lineProps', {'-', 'LineWidth', 2, 'Color', [0.3020, 0.7451, 0.9333]}); hold on
    
    % (late) opto on trace
    inds = find(optoLevel == i);
    inds = inds(4:5);
    BBNFR = mean(optoResponses.B(inds, :));
    BBNSTD = std(optoResponses.B(inds, :));
    [sortedLoudnessBBN, sortInds] = sort(loudness.RLFBBN);
    sortedBBNFR = BBNFR(sortInds);
    sortedBBNSTD = BBNSTD(sortInds);
    
    % plot late opto on traces
    shadedErrorBar(sortedLoudnessBBN, sortedBBNFR, sortedBBNSTD,...
        'lineProps', {'-', 'LineWidth', 2, 'Color', [0.9490, 0.3333, 0.3333]}); hold on
    
    legend('opto off', 'opto on (early)', 'opto on (late)');
end


%% Plot - different level opto effect on spontaneous firing rate 

optoOnOffTimes = getOptoSpFROnOffTimes(spike);

slideWindow = 30000;
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
    
    ephysChunkData = getVoltage(data.Data.Data(outputChannel, ephysStartInd:ephysStopInd));
    ephysChunkTime = ephysTime(ephysStartInd:ephysStopInd);
    
    [~, outputSpikeInds, ~] = crossdet(ephysChunkData, outputThresh, 'thresholdDown');

    
    spikes = zeros(size(ephysChunkData));
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

optoLevel = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4];
laserPower = [1.5, 3.0, 5.0, 7.0];
figure('Color', 'white','WindowState','maximized'); clf;
colors = cool(max(optoLevel));

for i = 1:max(optoLevel)
    
    inds = find(optoLevel == i);
    
    levelFR = outputFRInterped(inds, :);
    
    shadedErrorBar(linspace(-5, 25, 301), mean(levelFR), std(levelFR),...
        'lineProps', {'-', 'LineWidth', 2, 'Color', colors(i, :)}); hold on
    
    text(21.5, 62-(i-1)*2, ['laser power = ', num2str(laserPower(i)), 'mW'],...
        'Color', colors(i, :), 'FontSize', 10)
    
    
end

plot([0, 20], [62, 62], '-', 'LineWidth', 2, 'Color',  [0.3, 0.74, 0.933]);
text(9, 64, 'OPTO ON', 'Color', [0.3, 0.74, 0.933], 'FontWeight', 'bold');
ylim([0, 70]);
ylabel('Spike rate (spks/s)');
xlabel('Time(s)');
title([session, ' Ch', num2str(outputChannel), ' Output Cell, Opto Effect on Spon. FR'], 'Interpreter', 'none');
