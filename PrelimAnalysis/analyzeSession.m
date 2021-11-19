function analyzeSession(session, varargin)
% settings
s.analyze = 'all';  % use 'spike' to only analyze spike files, use 'camera' to only analyze camera files. 'all' will analyze both.
s.foodMinInterval = 5;  % unit: sec
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs


% initializations
rootFolder = 'Z:\Qianyun\DCN\';
gitFolder = 'D:\DCN_Project\Github\DCN\';
sessionFolder = fullfile(rootFolder, 'Data', session);

% determine the content to analyze
switch s.analyze
    case 'all'
        analyzeSpike = true;
        analyzeCamera = true;
    case 'spike'
        analyzeSpike = true;
        analyzeCamera = false;
    case 'camera'
        analyzeSpike = false;
        analyzeCamera = true;
end


% format camera-related files into cameraAnalysed.mat
% including video tracking, converted camera timestamps (into spike time)
if analyzeCamera
    fprintf('%s : analyzing camera files', session)
    
    % check if the session folder has all the files needed
    if ~exist(fullfile(rootFolder, 'Data', session, 'cameraTimestamps.csv'), 'file') || ~exist(fullfile(rootFolder, 'Data', session, 'Spike.mat'), 'file')
        error('cameraTimestamps.csv or Spike.mat is missing! Can not proceed!');
    elseif ~exist(fullfile(rootFolder, 'Data', session, 'behaviorVedio00_tracking.csv'), 'file')
        error('Behavior video tracking is missing! Can not proceed!');
    end
    
    % load video tracking spreadsheet
    videoTracking = readtable(fullfile(rootFolder, 'Data', session, 'behaviorVedio00_tracking.csv'));
    videoTracking.Properties.VariableNames{1} = 'frameNum';
    videoTracking.frameNum = videoTracking.frameNum + 1;
    
    % get converted camera frame times
    frameTimestamps = convertCameraToSpikeTimes(session);
    
    if size(videoTracking, 1) ~= length(frameTimestamps)
        error('Frame num in videoTracking table does NOT match that in frameTimestamps! Cannot proceed!');
    end
    
    % add frameTimestamps into videoTracking Table
    videoTracking.frameTimestamps = frameTimestamps;
        
    save(fullfile(sessionFolder, 'cameraAnalyzed.mat'), 'videoTracking');
    fprintf('cameraAnalyzed.mat saved in %s', sessionFolder)
end


% format spike-related files into spikeAnalysed.mat
% including raw data of channel [audioSig, Food, Mic, cameraTrigger, ],
% timestamps for all the raw data
if analyzeSpike
    fprintf('%s : analyzing spike files\n', session)
    
    % check if the session folder has all the files needed
    if ~exist(fullfile(rootFolder, 'Data', session, 'Spike.mat'), 'file')
        error('Spike.mat is missing! Can not proceed!');
    end
    
    % load spike.mat
    disp('loading spike.mat...');
    spikeTemp = load(fullfile(sessionFolder, 'Spike.mat'));
    disp('finish loading, start formating...')
    
    spike = struct();
    
    spike.audioSignal = spikeTemp.AudioSig.values;
    spike.audioSignalTimes = spikeTemp.AudioSig.times;
    
    spike.micSignal = spikeTemp.Mic.values;
    spike.micSignalTimes = spikeTemp.Mic.times;
    
    spike.cameraTriggerTimes = spikeTemp.CamTrig.times;
    spike.audioTriggerTimes = spikeTemp.AudioSyn.times;
    
    % clean the food trigger - get rid of false trigger (trigger but no
    % food delivered)
    spike.foodTriggerTimes = spikeTemp.Food.times;
    inds = find(diff(spike.foodTriggerTimes)<s.foodMinInterval);
    spike.foodTriggerTimes(inds) = [];
    spike.totalFoodNum = length(spike.foodTriggerTimes);
    
    spike.keyboardInput = char(spikeTemp.Keyboard.codes(:, 1));
    spike.keyboardTimes = spikeTemp.Keyboard.times;
    
    save(fullfile(sessionFolder, 'spikeAnalyzed.mat'), 'spike', '-v7.3');
    fprintf('spikeAnalyzed.mat saved in %s\n', sessionFolder)    
end


end