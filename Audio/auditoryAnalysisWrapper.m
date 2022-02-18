%% settings

sessionList = {'20211111_002', '20211111_003'};
ephysChannelList = [ 24, 24];

crudeStartFreq = 4500; 
crudeStopFreq = 72000; 
crudeSteps = 41;
crudeFreqs = logspace(log10(crudeStartFreq), log10(crudeStopFreq), crudeSteps); 

fineStartFreq = 7500;
fineStopFreq = 30000;
fineSteps = 101;
fineFreqs = logspace(log10(fineStartFreq), log10(fineStopFreq), fineSteps); 

freq.RLFBF = 21000;
load('C:\Users\qyzzz\OneDrive\AudioFiles\RateLevelFunction\BF21K_10to70DBSPL_1DBSPLInterval_RateLevelFunction_LoudnessLevel.mat');
loudness.RLFBF = loudnessLevel;
clear loudnessLevel;

load('C:\Users\qyzzz\OneDrive\AudioFiles\RateLevelFunction\BBNNEW_10to70DBSPL_1DBSPLInterval_RateLevelFunction_BBNLoudnessLevel.mat');
loudness.RLFBBN = BBNLoudnessLevel;


freq.J = crudeFreqs;
freq.K = crudeFreqs;
freq.L = crudeFreqs;
freq.H = crudeFreqs;

loudness.J = 70; % unit: DBSPL;
loudness.K = 60;
loudness.L = 50;
loudness.H = 40;

freq.U = fineFreqs;
freq.O = fineFreqs;
freq.P = fineFreqs;
freq.Y = fineFreqs;

loudness.U = 40;
loudness.O = 30;
loudness.P = 20;
loudness.Y = 15;

%% process data 

for i = 1:length(sessionList)
    session = sessionList{i};
    ephysChannel = ephysChannelList(i);
    
    if ~exist('responses', 'var')
        responses = processSessionAuditoryData(session, freq, loudness, ephysChannel);
    else
        responses = processSessionAuditoryData(session, freq, loudness, ephysChannel, responses);
    end
    
end

%% plot response map

% plot response map!!
figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
colorMatrix = {'#0072BD', 	'#EDB120', 	'#77AC30', 	'#7E2F8E', 	'#A2142F', '#D95319', '#4DBEEE' };
rows = length(colorMatrix);
cols = 1;

keys = fieldnames(responses);
keys(5) = [];

for i = 1:length(keys)
    
    subplot(rows, cols, i);
    
    x = freq.(keys{i})/1000;
    y = responses.(keys{i});

    semilogx(x, y, '-', 'LineWidth', 1.5, 'Color', colorMatrix{i}); hold on;
    % plot(x, repmat(meanBaselineFR, size(x)), '--', 'LineWidth', 0.5, 'Color', colorMatrix{i});
    box off; 
    legend([num2str(loudness.(keys{i})) ' dB SPL']);
    legend boxoff;
    
    if i ~= length(keys)
        xlim([freq.J(1)/1000, freq.J(end)/1000]);
        h = gca;
        h.YAxis.Visible = 'off';
        h.XAxis.Visible = 'off';
        
    else
        xticks([4.5, 9, 18, 36, 72]);
        xlim([freq.J(1)/1000, freq.J(end)/1000]);
        ylim([0, 120])
        xlabel('frequency (kHz)');
        ylabel('spikes/s');
    end
    
        
    if i == 1
        title(['Response Map ', session, ' Ch ', num2str(ephysChannel)],...
            'Interpreter', 'none');
    end
    
end

%% plot heatmaps ver. of the response map
% dim 1 -> freq
% dim 2 -> loudness
% dim 3 -> response

% crude search
keys = {'J', 'K', 'L', 'H'};

responseData = nan(crudeSteps*length(keys), 3);
for i = 1:length(keys)
    for j = 1:crudeSteps
        responseData(crudeSteps*(i-1)+j, :) = [freq.(keys{i})(j), loudness.(keys{i}),...
            mean(responses.(keys{i})(:, j))];
    end
end

Freqs = logspace(log10(crudeStartFreq), log10(crudeStopFreq), crudeSteps*10);
Levels = 70:-1:40;
[F,L]=meshgrid(Freqs, Levels);

Vq = griddata(responseData(:,1), responseData(:,2),responseData(:,3),F,L,'cubic');

figure('Color', 'white','WindowState','maximized'); clf
imagesc(Freqs(:)/1000,Levels(:),Vq);
colorbar
set(gca,'YDir','normal');
set(gca, 'XScale', 'log');

xticks([4, 8, 16, 32, 64])
xticklabels({'4', '8', '16', '32', '64'})
h=gca; h.XAxis.TickLength = [0 0];
ylabel('Sound Level (dB SPL)');
xlabel('Frequency (kHz)');
title([session, ' Output Cell Ch', num2str(ephysChannel)], 'Interpreter', 'none');


% fine search
keys = {'U', 'O', 'P', 'Y'};

responseData = nan(fineSteps*length(keys), 3);
for i = 1:length(keys)
    for j = 1:crudeSteps
        responseData(fineSteps*(i-1)+j, :) = [freq.(keys{i})(j), loudness.(keys{i}),...
            mean(responses.(keys{i})(:, j))];
    end
end

Freqs = logspace(log10(fineStartFreq), log10(fineStopFreq), fineSteps*10);
Levels = 40:-1:15;
[F,L]=meshgrid(Freqs, Levels);

Vq = griddata(responseData(:,1), responseData(:,2),responseData(:,3),F,L,'cubic');

figure('Color', 'white','WindowState','maximized'); clf
imagesc(Freqs(:)/1000,Levels(:),Vq);
colorbar
set(gca,'YDir','normal');
set(gca, 'XScale', 'log');

xticks([15 20 25 30])
xticklabels({'15', '20', '25', '30'})
h=gca; h.XAxis.TickLength = [0 0];
ylabel('Sound Level (dB SPL)');
xlabel('Frequency (kHz)');
title([session, ' Output Cell Ch', num2str(ephysChannel)], 'Interpreter', 'none');

%% plot LFPs during chewing mimic 
clear all; 
clc;

session = '20211212_002';
LFPChannel = 30;

% settings
totalDuration = 5*1000; % ms
oneMimicDuration = 30; % ms
fastBBNDuration = 2; % ms
mimicRepeat = 20;

DBAttenuation = 0:-3:-27;

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
    % function to extract voltage from binary file
    getVoltage = @(data) ...
        double(data)*sessionEphysInfo.bitVolts; % extract voltage from memmapfile, converting to votlage, highpassing, and only return specific channel
    
    % load raw ephys voltage
    contFiles = dir(fullfile(sessionEphysInfo.ephysFolder, '*.continuous'));
    data = memmapfile(fullfile(sessionEphysInfo.ephysFolder, [contFiles(1).name(1:end-12), 's.dat']), ...
        'Format', {'int16', [sessionEphysInfo.channelNum sessionEphysInfo.smps], 'Data'}, 'Writable', false);
end

%---------------------------processing the data---------------------------%
stimuliTimes = [];
ephysStartInds = [];
ephysStopInds = [];
for i = 1:length(spike.keyboardInput)
    if strcmp(spike.keyboardInput(i), 'C')
        stimuliTimes = [stimuliTimes, spike.keyboardTimes(i)];
        ind = find(sessionEphysInfo.convertedEphysTimestamps >= spike.keyboardTimes(i), 1, 'first');
        ephysStartInds = [ephysStartInds, ind];
        ind2 = find(sessionEphysInfo.convertedEphysTimestamps <= spike.keyboardTimes(i)+5, 1, 'last');
        ephysStopInds = [ephysStopInds, ind2];
    end
end


if length(stimuliTimes) ~= length(DBAttenuation)
    error('Stimuli times does NOT match attenuation levels!');
end

figure('Color', 'white','WindowState','maximized'); clf
cols = 5;
rows = 2;
plotInd = 1;
for i = 1:length(stimuliTimes)
    chunkEphysData = getVoltage(data.Data.Data(channelNum_OpenEphys(LFPChannel),  ephysStartInds(i):ephysStopInds(i)));
    
    subplot(rows, cols, plotInd);
    x = sessionEphysInfo.convertedEphysTimestamps(ephysStartInds(i):ephysStopInds(i));
    plot(x, chunkEphysData);
    box off; axis tight;
    ylim([-300, 300]);
    xlabel('time (sec)');
    ylabel('microVolt');
    title({[' LFP Ch', num2str(LFPChannel)]}, {['chewing mimic (DB Attenuation = ', num2str(DBAttenuation(i)), ' )']},...
        'Interpreter', 'none'); 
    plotInd = plotInd + 1;
end

    











