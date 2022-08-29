%% settings

sessionList = { '20220823_001' };
ephysChannelList = [31, 31];

crudeStartFreq = 4500; 
crudeStopFreq = 72000; 
crudeSteps = 101;
crudeFreqs = logspace(log10(crudeStartFreq), log10(crudeStopFreq), crudeSteps); 

fineStartFreq = 4500;
fineStopFreq = 72000;
fineSteps = 101;
fineFreqs = logspace(log10(fineStartFreq), log10(fineStopFreq), fineSteps); 

freq.RLFBF = 19500;
load('Z:\Qianyun\DCN\AudioFiles\RateLevelFunction\BF19.5KHz_5to75DBSPL_1DBSPLInterval_RateLevelFunction.mat');
loudness.RLFBF = randomLevel;
clear loudnessLevel;
% loudness.RLFBF = 10:1:70;

% load('C:\Users\Qianyun Zhang\OneDrive\AudioFiles\RateLevelFunction\BBNNEW_10to70DBSPL_1DBSPLInterval_RateLevelFunction_BBNLoudnessLevel.mat');
load('Z:\Qianyun\DCN\AudioFiles\RateLevelFunction\BBNNEW_5to75DBSPL_1DBSPLInterval_RateLevelFunction.mat');
loudness.RLFBBN = randomLevel;
% loudness.RLFBBN = 10:1:70;

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

loudness.U = 30;
loudness.O = 20;
loudness.P = 10;
loudness.Y = 5;


if ~exist('responses', 'var')
    responses = [];
end

if ~exist('optoResponses', 'var')
    optoResponses = [];
end

if ~exist('intervalResponses', 'var')
    intervalResponses = [];
end

if ~exist('optoIntervalResponses', 'var')
    optoIntervalResponses = [];
end



%% process data 

ephysThresh = -190;
for i = 1:length(sessionList)
    session = sessionList{i};
    ephysChannel = ephysChannelList(i);
    
    [responses, optoResponses, intervalResponses, optoIntervalResponses] = processSessionAuditoryData(session, freq, loudness, ...
        ephysChannel, responses, optoResponses, intervalResponses, optoIntervalResponses,...
        'threshold', ephysThresh, 'plotType', [], 'thresholdType', 'thresholdDown');  
end


%% save the responses & optoResponses data

fileFolder = fullfile('Z:\Qianyun\DCN\Data', session);

% save responses.mat
if exist('responses', 'var')
    save(fullfile(fileFolder, 'responses.mat'), 'responses');
end

% save optoResponses.mat
if exist('optoResponses', 'var')
    save(fullfile(fileFolder, 'optoResponses.mat'), 'optoResponses');
end

% save intervalResponses.mat
if exist('intervalResponses', 'var')
    save(fullfile(fileFolder, 'intervalResponses.mat'), 'intervalResponses');
end

% save optoIntervalResponses.mat
if exist('optoIntervalResponses', 'var')
    save(fullfile(fileFolder, 'optoIntervalResponses.mat'), 'optoIntervalResponses');
end

%% plot heatmaps ver. of the response map (normal exp, no opto at all; if opto exp, use the next block of code!)
% dim 1 -> freq
% dim 2 -> loudness
% dim 3 -> response

% crude search
keys = {'K', 'L', 'H', 'U', 'O', 'P'};
%
responseData = nan(crudeSteps*length(keys), 3);
for i = 1:length(keys)
    for j = 1:crudeSteps
        responseData(crudeSteps*(i-1)+j, :) = [freq.(keys{i})(j), loudness.(keys{i}),...
            mean(responses.(keys{i})(:, j))];
    end
end

Freqs = logspace(log10(crudeStartFreq), log10(crudeStopFreq), crudeSteps*10);
Levels = 60:-1:10;
[F,L]=meshgrid(Freqs, Levels);

Vq = griddata(responseData(:,1), responseData(:,2),responseData(:,3),F,L,'cubic');

figure('Color', 'white','WindowState','maximized'); clf;
% rows = 2; cols = 2;
% plotInd = 1;
% subplot(rows, cols, plotInd);
imagesc(Freqs(:)/1000,Levels(:),Vq);
colorbar
set(gca,'YDir','normal');
set(gca, 'XScale', 'log'); hold on;

xticks([4, 8, 16, 32, 64])
xticklabels({'4', '8', '16', '32', '64'})
h=gca; h.XAxis.TickLength = [0 0];
ylabel('Sound Level (dB SPL)');
xlabel('Frequency (kHz)');
title([session, ' Output Cell Ch', num2str(ephysChannel), ' Opto Off'], 'Interpreter', 'none');

% % fine search
% keys = {'U', 'O', 'P', 'Y'};
% 
% responseData = nan(fineSteps*length(keys), 3);
% for i = 1:length(keys)
%     for j = 1:crudeSteps
%         responseData(fineSteps*(i-1)+j, :) = [freq.(keys{i})(j), loudness.(keys{i}),...
%             mean(responses.(keys{i})(:, j))];
%     end
% end
% 
% Freqs = logspace(log10(fineStartFreq), log10(fineStopFreq), fineSteps*10);
% Levels = 35:-1:10;
% [F,L]=meshgrid(Freqs, Levels);
% 
% Vq = griddata(responseData(:,1), responseData(:,2),responseData(:,3),F,L,'cubic');
% 
% figure('Color', 'white','WindowState','maximized'); clf
% imagesc(Freqs(:)/1000,Levels(:),Vq);
% colorbar
% set(gca,'YDir','normal');
% set(gca, 'XScale', 'log');
% 
% xticks([4, 6, 8, 10])
% xticklabels({'4', '6', '8', '10'})
% h=gca; h.XAxis.TickLength = [0 0];
% ylabel('Sound Level (dB SPL)');
% xlabel('Frequency (kHz)');
% % title([session, ' Output Cell Ch', num2str(ephysChannel)], 'Interpreter', 'none');
% 

% RLF for BBN
% responses = respnosesBackup;
BBNFR = mean(rmmissing(responses.B));
BBNSTD = std(rmmissing(responses.B));
[sortedLoudnessBBN, inds] = sort(loudness.RLFBBN);
sortedBBNFR = BBNFR(inds);
sortedBBNSTD = BBNSTD(inds);

% apply a 3-point trigular filter
filteredBBNFR = nan(size(sortedBBNFR));
for i = 1:length(sortedBBNFR)
    if i == 1
        filteredBBNFR(i) = (2*sortedBBNFR(i) + sortedBBNFR(i+1))/3;
    elseif i == length(sortedBBNFR)
        filteredBBNFR(i) = (2*sortedBBNFR(i) + sortedBBNFR(i-1))/3;
    else
        filteredBBNFR(i) = (2*sortedBBNFR(i) + sortedBBNFR(i-1) + sortedBBNFR(i+1))/4;
    end
end

% RLF plot for BF
% responses = responsesBFRLF2;
if ~isempty(responses.E)
    inds = [];
    BFFR = responses.E;
    BFFR = mean(responses.E);
    BFSTD = std(responses.E);
    [sortedLoudnessBF, inds] = sort(loudness.RLFBF);
    sortedBFFR = BFFR(inds);
    sortedBFSTD = BFSTD(inds);
    
    % apply a 3-point trigular filter
    filteredBFFR = nan(size(sortedBFFR));
    for i = 1:length(sortedBFFR)
        if i == 1
            filteredBFFR(i) = (2*sortedBFFR(i) + sortedBFFR(i+1))/3;
        elseif i == length(sortedBFFR)
            filteredBFFR(i) = (2*sortedBFFR(i) + sortedBFFR(i-1))/3;
        else
            filteredBFFR(i) = (2*sortedBFFR(i) + sortedBFFR(i-1) + sortedBFFR(i+1))/4;
        end
    end
end



% Plot 2 - RLF Plots, with shaded error bar
baselineFR = responses.baselineFR;

figure('Color', 'white','WindowState','maximized'); clf;
shadedErrorBar(sortedLoudnessBBN, sortedBBNFR, sortedBBNSTD, 'lineProps', {'-', 'LineWidth', 2});
% plot(sortedLoudnessBBN, responses.B, '-', 'LineWidth', 2)

hold on; box off; axis tight;
plot([sortedLoudnessBBN(1), sortedLoudnessBBN(end)], [baselineFR, baselineFR],...
    '--', 'Color', [0.8 0.8 0.8], 'LineWidth', 1.5); 
xlim([sortedLoudnessBBN(1), sortedLoudnessBBN(end)]);
xlabel('loudness (dB SPL)');
ylabel('rate (spk/s)');
legend('BBN', 'baseline');
title([session, ' Ch', num2str(ephysChannel), ' BBN RLF'], 'Interpreter', 'none');

if ~isempty(responses.E)
    shadedErrorBar(sortedLoudnessBF, sortedBFFR, sortedBFSTD, 'lineProps', {'-k', 'LineWidth', 2});
    xlim([sortedLoudnessBF(1), sortedLoudnessBF(end)]);
    legend('BBN', 'baseline', ['BF = ', num2str(freq.RLFBF/1000) ,'KHz']);
end

% Plot 3 - RLF Plots with 3-point triangular filter

figure('Color', 'white','WindowState','maximized'); clf;
plot(sortedLoudnessBBN, filteredBBNFR, '-', 'LineWidth', 2);

hold on; box off;
plot([sortedLoudnessBBN(1), sortedLoudnessBBN(end)], [baselineFR, baselineFR],...
    '--', 'Color', [0.8 0.8 0.8], 'LineWidth', 1.5); hold on;
xlim([sortedLoudnessBBN(1), sortedLoudnessBBN(end)]);
xlabel('loudness (dB SPL)');
ylabel('rate (spk/s)');
ylim([0, 160]);
xlim([sortedLoudnessBBN(1), sortedLoudnessBBN(end)])
legend('BBN', 'baseline');
% title([session, ' Ch', num2str(ephysChannel), ' RLFs'], 'Interpreter', 'none');

if ~isempty(responses.E)
    plot(sortedLoudnessBF, filteredBFFR, '-k', 'LineWidth', 2); 
    legend('BBN', 'baseline', ['BF = ', num2str(freq.RLFBF/1000) ,'KHz']);
end


%% plot OPTO ON response map 
% dim 1 -> freq
% dim 2 -> loudness
% dim 3 -> response

% crude search
keys = { 'K', 'L', 'H', 'U', 'O', 'P'};
%
responseData = nan(crudeSteps*length(keys), 3);
for i = 1:length(keys)
    for j = 1:crudeSteps
        responseData(crudeSteps*(i-1)+j, :) = [freq.(keys{i})(j), loudness.(keys{i}),...
            nanmean(optoResponses.(keys{i})(:, j))];
    end
end

Freqs = logspace(log10(crudeStartFreq), log10(crudeStopFreq), crudeSteps*10);
Levels = 60:-1:10;
[F,L]=meshgrid(Freqs, Levels);

Vq = griddata(responseData(:,1), responseData(:,2),responseData(:,3),F,L,'cubic');

figure('Color', 'white','WindowState','maximized'); clf;
% rows = 2; cols = 2;
% plotInd = 1;
% subplot(rows, cols, plotInd);
imagesc(Freqs(:)/1000,Levels(:),Vq);
colorbar
set(gca,'YDir','normal');
set(gca, 'XScale', 'log');

xticks([4, 8, 16, 32, 64])
xticklabels({'4', '8', '16', '32', '64'})
h=gca; h.XAxis.TickLength = [0 0];
ylabel('Sound Level (dB SPL)');
xlabel('Frequency (kHz)');
title([session, ' Output Cell Ch', num2str(ephysChannel), ' Opto On'], 'Interpreter', 'none');


%% process BBN RLF data for both opto on and off condition
% BBN - OPTO OFF
baselineFR = responses.baselineFR;

BBNFR = mean(rmmissing(responses.B));
BBNSTD = std(rmmissing(responses.B));
[sortedLoudnessBBN, inds] = sort(loudness.RLFBBN);
sortedBBNFR = BBNFR(inds);
sortedBBNSTD = BBNSTD(inds);

% apply a 3-point trigular filter
filteredBBNFR = nan(size(sortedBBNFR));
for i = 1:length(sortedBBNFR)
    if i == 1
        filteredBBNFR(i) = (2*sortedBBNFR(i) + sortedBBNFR(i+1))/3;
    elseif i == length(sortedBBNFR)
        filteredBBNFR(i) = (2*sortedBBNFR(i) + sortedBBNFR(i-1))/3;
    else
        filteredBBNFR(i) = (2*sortedBBNFR(i) + sortedBBNFR(i-1) + sortedBBNFR(i+1))/4;
    end
end

% BBN - OPTO ON
OPTOBBNFR = mean(rmmissing(optoResponses.B));
OPTOBBNSTD = std(rmmissing(optoResponses.B));
[sortedLoudnessBBN, inds] = sort(loudness.RLFBBN);
sortedOPTOBBNFR = OPTOBBNFR(inds);
sortedOPTOBBNSTD = OPTOBBNSTD(inds);

% apply a 3-point trigular filter
filteredOPTOBBNFR = nan(size(sortedOPTOBBNFR));
for i = 1:length(sortedOPTOBBNFR)
    if i == 1
        filteredOPTOBBNFR(i) = (2*sortedOPTOBBNFR(i) + sortedOPTOBBNFR(i+1))/3;
    elseif i == length(sortedOPTOBBNFR)
        filteredOPTOBBNFR(i) = (2*sortedOPTOBBNFR(i) + sortedOPTOBBNFR(i-1))/3;
    else
        filteredOPTOBBNFR(i) = (2*sortedOPTOBBNFR(i) + sortedOPTOBBNFR(i-1) + sortedOPTOBBNFR(i+1))/4;
    end
end


% intervals - OPTO OFF
BBNFRi = mean(rmmissing(intervalResponses.B));
BBNSTDi = std(rmmissing(intervalResponses.B));
[sortedLoudnessBBN, inds] = sort(loudness.RLFBBN);
sortedBBNFRi = BBNFRi(inds);
sortedBBNSTDi = BBNSTDi(inds);

% apply a 3-point trigular filter
filteredBBNFRi = nan(size(sortedBBNFRi));
for i = 1:length(sortedBBNFRi)
    if i == 1
        filteredBBNFRi(i) = (2*sortedBBNFRi(i) + sortedBBNFRi(i+1))/3;
    elseif i == length(sortedBBNFRi)
        filteredBBNFRi(i) = (2*sortedBBNFRi(i) + sortedBBNFRi(i-1))/3;
    else
        filteredBBNFRi(i) = (2*sortedBBNFRi(i) + sortedBBNFRi(i-1) + sortedBBNFRi(i+1))/4;
    end
end


% intervals - OPTO ON
OPTOBBNFRi = mean(rmmissing(optoIntervalResponses.B));
OPTOBBNSTDi = std(rmmissing(optoIntervalResponses.B));
[sortedLoudnessBBN, inds] = sort(loudness.RLFBBN);
sortedOPTOBBNFRi = OPTOBBNFRi(inds);
sortedOPTOBBNSTDi = OPTOBBNSTDi(inds);

% apply a 3-point trigular filter
filteredOPTOBBNFRi = nan(size(sortedOPTOBBNFRi));
for i = 1:length(sortedOPTOBBNFRi)
    if i == 1
        filteredOPTOBBNFRi(i) = (2*sortedOPTOBBNFRi(i) + sortedOPTOBBNFRi(i+1))/3;
    elseif i == length(sortedOPTOBBNFRi)
        filteredOPTOBBNFRi(i) = (2*sortedOPTOBBNFRi(i) + sortedOPTOBBNFRi(i-1))/3;
    else
        filteredOPTOBBNFRi(i) = (2*sortedOPTOBBNFRi(i) + sortedOPTOBBNFRi(i-1) + sortedOPTOBBNFRi(i+1))/4;
    end
end

%% process BF RLF data for both opto on and off condition
% BF - OPTO OFF
if ~isempty(responses.E)
    clear inds
    BFFR = mean(rmmissing(responses.E));
    BFSTD = std(rmmissing(responses.E));
    [sortedLoudnessBF, inds] = sort(loudness.RLFBF);
    sortedBFFR = BFFR(inds);
    sortedBFSTD = BFSTD(inds);
    
    % apply a 3-point trigular filter
    filteredBFFR = nan(size(sortedBFFR));
    for i = 1:length(sortedBFFR)
        if i == 1
            filteredBFFR(i) = (2*sortedBFFR(i) + sortedBFFR(i+1))/3;
        elseif i == length(sortedBFFR)
            filteredBFFR(i) = (2*sortedBFFR(i) + sortedBFFR(i-1))/3;
        else
            filteredBFFR(i) = (2*sortedBFFR(i) + sortedBFFR(i-1) + sortedBFFR(i+1))/4;
        end
    end
end




% BF - OPTO ON
if ~isempty(optoResponses.E)
    clear inds
    OPTOBFFR = mean(rmmissing(optoResponses.E));
    OPTOBFSTD = std(rmmissing(optoResponses.E));
    [sortedLoudnessBF, inds] = sort(loudness.RLFBF);
    sortedOPTOBFFR = OPTOBFFR(inds);
    sortedOPTOBFSTD = OPTOBFSTD(inds);
    
    % apply a 3-point trigular filter
    filteredOPTOBFFR = nan(size(sortedOPTOBFFR));
    for i = 1:length(sortedOPTOBFFR)
        if i == 1
            filteredOPTOBFFR(i) = (2*sortedOPTOBFFR(i) + sortedOPTOBFFR(i+1))/3;
        elseif i == length(sortedOPTOBFFR)
            filteredOPTOBFFR(i) = (2*sortedOPTOBFFR(i) + sortedOPTOBFFR(i-1))/3;
        else
            filteredOPTOBFFR(i) = (2*sortedOPTOBFFR(i) + sortedOPTOBFFR(i-1) + sortedOPTOBFFR(i+1))/4;
        end
    end
end

startLoudness = 5;
stopLoudness =  75;
tempInds = sortedLoudnessBBN >= startLoudness & sortedLoudnessBBN <= stopLoudness;
sortedLoudnessBBN = sortedLoudnessBBN(tempInds);
sortedBBNFR = sortedBBNFR(tempInds);
sortedBBNSTD = sortedBBNSTD(tempInds);
sortedOPTOBBNFR = sortedOPTOBBNFR(tempInds);
sortedOPTOBBNSTD = sortedOPTOBBNSTD(tempInds);
filteredBBNFR = filteredBBNFR(tempInds);
filteredBFFR = filteredBFFR(tempInds);
clear tempInds
tempInds = sortedLoudnessBF >= startLoudness & sortedLoudnessBF <= stopLoudness;
sortedLoudnessBF = sortedLoudnessBF(tempInds);
sortedBFFR = sortedBFFR(tempInds);
sortedBFSTD = sortedBFSTD(tempInds);
sortedOPTOBFFR = sortedOPTOBFFR(tempInds);
sortedOPTOBFSTD = sortedOPTOBFSTD(tempInds);
filteredOPTOBBNFR = filteredOPTOBBNFR(tempInds);
filteredOPTOBFFR = filteredOPTOBFFR(tempInds);

%%

% PLOT!!!
% Plot 2 - BBN RLF Plots, with shaded error bar. OPTO ON vs. OPTO OFF
figure('Color', 'white','WindowState','maximized'); clf;
shadedErrorBar(sortedLoudnessBBN, sortedBBNFR, sortedBBNSTD, 'lineProps', {'-k', 'LineWidth', 2});
hold on; box off; axis tight;

shadedErrorBar(sortedLoudnessBBN, sortedOPTOBBNFR, sortedOPTOBBNSTD, 'lineProps', {'-', 'LineWidth', 2, 'Color', [0.79, 0.57, 0.92]});

plot([sortedLoudnessBBN(1), sortedLoudnessBBN(end)], [baselineFR, baselineFR],...
    '--', 'Color', [0.8 0.8 0.8], 'LineWidth', 1.5); 

xlim([sortedLoudnessBBN(1), sortedLoudnessBBN(end)]);
xlabel('loudness (dB SPL)');
ylabel('rate (spk/s)');
legend('BBN (OPTO OFF)', 'BBN (OPTO ON)', 'baseline');
title([session, ' Ch', num2str(ephysChannel), ' BBN RLF'], 'Interpreter', 'none');

% PLOT!!!
% Plot 2.5 - interval BBN RLF Plots, with shaded error bar. OPTO ON vs. OPTO OFF
figure('Color', 'white','WindowState','maximized'); clf;
shadedErrorBar(sortedLoudnessBBN, sortedBBNFRi, sortedBBNSTDi, 'lineProps', {'-k', 'LineWidth', 2});
hold on; box off; axis tight;

shadedErrorBar(sortedLoudnessBBN, sortedOPTOBBNFRi, sortedOPTOBBNSTDi, 'lineProps', {'-', 'LineWidth', 2, 'Color', [0.79, 0.57, 0.92]});

plot([sortedLoudnessBBN(1), sortedLoudnessBBN(end)], [baselineFR, baselineFR],...
    '--', 'Color', [0.8 0.8 0.8], 'LineWidth', 1.5); 

xlim([sortedLoudnessBBN(1), sortedLoudnessBBN(end)]);
xlabel('loudness (dB SPL)');
ylabel('rate (spk/s)');
legend('BBN Interval Rebound (OPTO OFF)', 'BBN Interval Rebound (OPTO ON)', 'baseline');
title([session, ' Ch', num2str(ephysChannel), ' Interval BBN RLF'], 'Interpreter', 'none');


%%
% Plot 3 - BF RLF Plots, with shaded error bar. OPTO ON vs. OPTO OFF
figure('Color', 'white','WindowState','maximized'); clf;
shadedErrorBar(sortedLoudnessBF, sortedBFFR, sortedBFSTD, 'lineProps', {'-k', 'LineWidth', 2});
hold on; box off; axis tight;

shadedErrorBar(sortedLoudnessBF, sortedOPTOBFFR, sortedOPTOBFSTD, 'lineProps', {'-', 'LineWidth', 2, 'Color', [0.79, 0.57, 0.92]});

plot([sortedLoudnessBF(1), sortedLoudnessBF(end)], [baselineFR, baselineFR],...
    '--', 'Color', [0.8 0.8 0.8], 'LineWidth', 1.5); 

xlim([sortedLoudnessBF(1), sortedLoudnessBF(end)]);
xlabel('loudness (dB SPL)');
ylabel('rate (spk/s)');
legend(['BF = ', num2str(freq.RLFBF/1000) ,'KHz (OPTO OFF)'],...
    ['BF = ', num2str(freq.RLFBF/1000) ,'KHz (OPTO ON)'], 'baseline');
title([session, ' Ch', num2str(ephysChannel), 'BF RLF'], 'Interpreter', 'none');


%%
% Plot 4 - BBN RLF Plot, triangular filtered, opto off vs. opto on
figure('Color', 'white','WindowState','maximized'); clf;
plot(sortedLoudnessBBN, filteredBBNFR, '-k', 'LineWidth', 2);
hold on; box off;
plot(sortedLoudnessBBN, filteredOPTOBBNFR, '-', 'LineWidth', 2, 'Color', [0.79, 0.57, 0.92]);

plot([sortedLoudnessBBN(1), sortedLoudnessBBN(end)], [baselineFR, baselineFR],...
    '--', 'Color', [0.8 0.8 0.8], 'LineWidth', 1.5); hold on;

xlim([sortedLoudnessBBN(1), sortedLoudnessBBN(end)]);
xlabel('loudness (dB SPL)');
ylabel('rate (spk/s)');
ylim([0, 100]);
xlim([sortedLoudnessBBN(1), sortedLoudnessBBN(end)])


legend('BBN (OPTO OFF)', 'BBN(OPTO ON)', 'baseline');
title([session, ' Ch', num2str(ephysChannel), ' BBN RLF'], 'Interpreter', 'none');


if ~isempty(responses.E)
    plot(sortedLoudnessBF, filteredBFFR, '-k', 'LineWidth', 2); 
    legend('BBN', 'baseline', ['BF = ', num2str(freq.RLFBF/1000) ,'KHz']);
end

%%
% Plot 5 - BF RLF Plot, triangular filtered, opto off vs. opto on

figure('Color', 'white','WindowState','maximized'); clf;
plot(sortedLoudnessBF, filteredBFFR, '-k', 'LineWidth', 2);
hold on; box off;
plot(sortedLoudnessBF, filteredOPTOBFFR, '-', 'LineWidth', 2, 'Color', [0.79, 0.57, 0.92]);

plot([sortedLoudnessBF(1), sortedLoudnessBF(end)], [baselineFR, baselineFR],...
    '--', 'Color', [0.8 0.8 0.8], 'LineWidth', 1.5); hold on;

xlim([sortedLoudnessBF(1), sortedLoudnessBF(end)]);
xlabel('loudness (dB SPL)');
ylabel('rate (spk/s)');
ylim([0, 100]);
xlim([sortedLoudnessBF(1), sortedLoudnessBF(end)])


legend(['BF = ', num2str(freq.RLFBF/1000) ,'KHz (OPTO OFF)'],...
    ['BF = ', num2str(freq.RLFBF/1000) ,'KHz (OPTO ON)'], 'baseline');
title([session, ' Ch', num2str(ephysChannel), ' BF RLF'], 'Interpreter', 'none');

