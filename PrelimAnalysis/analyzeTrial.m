function analyzeTrial(session, varargin)

% settings
s.hasMic = false; % whether this session contains good microphone recordings
s.crunchSearchTimeWindow = 5; % sec
s.crunchTimeWindow = 0.02; % sec
s.chewingSearchTimeWindow = 12; % sec, detemine the time range (how long) of the chewing search 
s.chewingSearchStartTimeShift = 4; % sec, how many seconds after the trial start should we strat the chewing search
s.chewingTimeWindow = 4; % sec
s.chewingSearchStepLength = 0.1; % sec
s.fs = 150000; % sampling rate for mic signal
s.ephysFs = 30000; % sampling rate for ephys

if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs

% initializations
rootFolder = 'Z:\Qianyun\DCN\';
gitFolder = 'D:\DCN_Project\Github\DCN\';
sessionFolder = fullfile(rootFolder, 'Data', session);

% load spikeAnalyzed.mat
if ~exist(fullfile(sessionFolder, 'spikeAnalyzed.mat'), 'file')
    error('spikeAnalyzed.mat NOT exist! Run analyzeSession.m first!');
else
    disp('loading spikeAnalyzed.mat...')
    load(fullfile(sessionFolder, 'spikeAnalyzed.mat'));
end

% load cameraAnalyzed.mat
if ~exist(fullfile(sessionFolder, 'cameraAnalyzed.mat'), 'file')
    error('cameraAnalyzed.mat NOT exist! Run analyzeSession.m first!');
else
    disp('loading cameraAnalyzed.mat...')
    load(fullfile(sessionFolder, 'cameraAnalyzed.mat'));
    if ~any(strcmp('jawDistance', videoTracking.Properties.VariableNames))
        error('JawDistance does NOT exist in videoTracking! Run analyzeSession.m first!');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%Start Processing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trialInfo = struct();

trialInfo.totalFoodTrialNumber = spike.totalFoodNum;
trialInfo.foodTriggerTimes = spike.foodTriggerTimes;

crunchWindow = [];
crunchMicRMSV = cell(1, 1);
chewingWindow = [];
chewingJawDistance = cell(1, 1);
chewingJawTime = cell(1, 1);
chewingJawInds = cell(1, 1);
selectedTrialNum = 1;

for i = 1:trialInfo.totalFoodTrialNumber
    fprintf('%d/%d', i, trialInfo.totalFoodTrialNumber);
    %%%%%%%%%%%%%%%%%%%%% Processing crunch %%%%%%%%%%%%%%%%%%%%%%%
    % locate big crunch chunk
    crunchWindow(i, 1) = spike.foodTriggerTimes(i) - 0.1;
    crunchWindow(i, 2) = crunchWindow(i, 1) + s.crunchSearchTimeWindow;
    
    % crunch chunk microphone
    if s.hasMic
        crunchChunkMicStartInd = find(spike.micSignalTimes >= crunchWindow(i, 1), 1, 'first');
        crunchChunkMicEndInd = find(spike.micSignalTimes <= crunchWindow(i, 2), 1, 'last');
        crunchChunkMic = spike.micSignal(crunchChunkMicStartInd : crunchChunkMicEndInd);
        crunchChunkMic = highpass(crunchChunkMic, 100, s.fs);
        crunchMicRMSV{i, 1} = sqrt(movmean(crunchChunkMic.^2, 100));
    end
    
    %%%%%%%%%%%%%%%%%%%%% Processing chewing %%%%%%%%%%%%%%%%%%%%%%%%
    trialStartTime = spike.foodTriggerTimes(i) - 0.1;
    chewingEndTime = trialStartTime + s.chewingSearchStartTimeShift + s.chewingSearchTimeWindow;
    
    % locate chewing period
    chewingChunkStartTime = trialStartTime + s.chewingSearchStartTimeShift;
    chewingChunkEndTime = chewingChunkStartTime + s.chewingTimeWindow;
    fitResults = [];
    
    while(chewingChunkEndTime <= chewingEndTime)
        chewingChunkVideoStartInd = find(videoTracking.frameTimestamps >= chewingChunkStartTime, 1, 'first');
        chewingChunkVideoEndInd = find(videoTracking.frameTimestamps <= chewingChunkEndTime, 1, 'last');
        
        chewingChunkJaw = videoTracking.jawDistance(chewingChunkVideoStartInd:chewingChunkVideoEndInd);
        fitResult = sineFit(1:length(chewingChunkJaw), chewingChunkJaw', 0);
        fitResults = [fitResults; fitResult];
        
        chewingChunkStartTime = chewingChunkStartTime + s.chewingSearchStepLength;
        chewingChunkEndTime = chewingChunkEndTime + s.chewingSearchStepLength;
    end
    
    fitMSE = fitResults(:, 5);
    tempInd = find(fitMSE == min(fitMSE(fitMSE>0)));
    minFitMSE = min(fitMSE(fitMSE>0));
    disp(['minFitMSE = ', num2str(minFitMSE)]);
    
    if fitResults(tempInd, 2) > 1.5 && fitResults(tempInd, 3) > 0.03
        disp(['trial ' num2str(i) ' selected'])
        chewingWindow(selectedTrialNum, 1) = trialStartTime + s.chewingSearchStartTimeShift + s.chewingSearchStepLength*(tempInd - 1);
        chewingWindow(selectedTrialNum, 2) = chewingWindow(selectedTrialNum, 1) + s.chewingTimeWindow;
        chewingChunkVideoStartInd = find(videoTracking.frameTimestamps >= chewingWindow(selectedTrialNum, 1), 1, 'first');
        chewingChunkVideoEndInd = find(videoTracking.frameTimestamps <= chewingWindow(selectedTrialNum, 2), 1, 'last');
        chewingChunkJaw = videoTracking.jawDistance(chewingChunkVideoStartInd:chewingChunkVideoEndInd);
        % sanity check
        sineFit(1:length(chewingChunkJaw), chewingChunkJaw'); % sanity check
        
        chewingJawDistance{selectedTrialNum, 1} = chewingChunkJaw;
        chewingJawTime{selectedTrialNum, 1} = videoTracking.frameTimestamps(chewingChunkVideoStartInd:chewingChunkVideoEndInd);
        chewingJawInds{selectedTrialNum, 1} = [chewingChunkVideoStartInd:chewingChunkVideoEndInd]';
        selectedTrialNum = selectedTrialNum + 1;
    end
end

% save variables into trialInfo structure
trialInfo.crunchWindow = crunchWindow;
trialInfo.crunchMicRMSV = crunchMicRMSV;
trialInfo.chewingWindow = chewingWindow;
trialInfo.chewingJawDistance = chewingJawDistance;
trialInfo.chewingJawTime = chewingJawTime;
trialInfo.chewingJawInds = chewingJawInds;

% plot and calculate jaw phase in each chew!
avgChewCycleLength = cell(size(chewingWindow, 1), 1);
trialInfo.chewingJawDistanceLowPass = cell(size(chewingWindow, 1), 1);
goodChewsStartStopTime = [];
goodChewsStartStopInds = [];
startInds = []; stopInds = [];
startTimes = []; stopTimes = [];
% figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
% rows = 5; cols = 6; plotInd = 1;
for i = 1:size(chewingWindow, 1)    
    fprintf('%d/%d', i, size(chewingWindow, 1)); 

    y = trialInfo.chewingJawDistance{i, 1};
    x = linspace(trialInfo.chewingWindow(i, 1), trialInfo.chewingWindow(i, 2), length(trialInfo.chewingJawDistance{i, 1}));
    y2 = lowpass(y, 5, 1/(x(2) - x(1)));
    trialInfo.chewingJawDistanceLowPass{i, 1} = y2;
    
    jawPhase = hilbert(lowpass(y, 5, 1/(x(2) - x(1))));
    jawPhase = angle(jawPhase);
    
    [~, locs] = findpeaks(jawPhase);
    [~, minlocs] = findpeaks(-jawPhase);
  
    figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
    subplot(2, 1, 1);    
    plot(x, y); hold on;
    plot(x, y2);    
    box off; axis tight;
    xlabel('spike timestamps (sec)');
    ylabel('jaw distance (pixel)');
    legend('raw data', 'lowpass 5Hz');
    title(['trial ', num2str(i)]);
    
    subplot(2, 1, 2);
    plot(x, jawPhase); hold on;
    plot(x(locs), jawPhase(locs), 'r.', 'MarkerSize', 10)    
    plot(x(minlocs), jawPhase(minlocs), 'b.', 'MarkerSize', 10);
    box off; axis tight;
    title('Phase');
    
    subplot(2, 1, 1); hold on;
    plot(x(locs), y(locs), 'r.', 'MarkerSize', 10);
    plot(x(minlocs), y(minlocs), 'b.', 'MarkerSize', 10);
    
    avgChewCycleLength{i, 1} = diff(minlocs);
    inds = find(diff(minlocs) >= 20 & diff(minlocs) <= 26); % temporarily hard code
    
    startInds = trialInfo.chewingJawInds{i, 1}(minlocs(inds));
    stopInds = trialInfo.chewingJawInds{i, 1}(minlocs(inds+1))-1;
    startTimes = videoTracking.frameTimestamps(startInds);
    stopTimes = videoTracking.frameTimestamps(stopInds);
    
    goodChewsStartStopInds = [goodChewsStartStopInds;[startInds, stopInds]];
    goodChewsStartStopTime = [goodChewsStartStopTime;[startTimes, stopTimes]];
   
end

goodChewsJawDistance = cell(size(goodChewsStartStopInds, 1), 1);
for i = 1:size(goodChewsStartStopInds, 1)
    goodChewsJawDistance{i, 1} = videoTracking.jawDistance(goodChewsStartStopInds(i, 1):goodChewsStartStopInds(i, 2));
end

trialInfo.goodChewsStartStopInds = goodChewsStartStopInds;
trialInfo.goodChewsStartStopTime = goodChewsStartStopTime;
trialInfo.goodChewsJawDistance = goodChewsJawDistance;

figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
for i = 1:length(trialInfo.goodChewsJawDistance)
    plot(trialInfo.goodChewsJawDistance{i, 1}); hold on;
end

save(fullfile(sessionFolder, 'trialInfo.mat'), 'trialInfo');    

end