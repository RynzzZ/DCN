function analyzeSession(session, varargin)
% settings
s.analyze = 'all';  % use 'spike' to only analyze spike files, use 'camera' to only analyze camera files. 'all' will analyze both.
s.foodMinInterval = 5;  % unit: sec
s.jawConfidenceThreshold = 0.5;
s.tongueConfidenceThreshold = 0.5;
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs


% initializations
rootFolder = 'Z:\Qianyun\DCN\';
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
    elseif ~exist(fullfile(rootFolder, 'Data', session, 'behaviorVideo00_tracking.csv'), 'file')
        error('Behavior video tracking is missing! Can not proceed!');
    end
    
    % load video tracking spreadsheet
    videoTracking = readtable(fullfile(rootFolder, 'Data', session, 'behaviorVideo00_tracking.csv'));
    videoTracking.Properties.VariableNames{1} = 'frameNum';
    videoTracking.frameNum = videoTracking.frameNum + 1;
    
    % get converted camera frame times
    frameTimestamps = convertCameraToSpikeTimes(session);
    
    if size(videoTracking, 1) ~= length(frameTimestamps)
        error('Frame num in videoTracking table does NOT match that in frameTimestamps! Cannot proceed!');
    end
    
    % add frameTimestamps into videoTracking Table
    videoTracking.frameTimestamps = frameTimestamps;
    
    % calculate jaw movement (pixel distance b/w upper jaw and lower jaw);
    % only include points with good tracking confidence
    ind = find(videoTracking.jaw_lower_confidence > s.jawConfidenceThreshold);
    videoTracking.jawDistance = nan(size(videoTracking, 1), 1);
    videoTracking.jawDistance(ind) = sqrt((videoTracking.jaw_upper_x(ind) - videoTracking.jaw_lower_x(ind)).^2 +...
        (videoTracking.jaw_upper_y(ind) - videoTracking.jaw_lower_y(ind)).^2 );
    
    % interpolate the jaw movement to get rid of nans, method = 'spline'
    inds = find(~isnan(videoTracking.jawDistance));
    videoTracking.jawDistance = interp1(inds, videoTracking.jawDistance(inds), 1:length(videoTracking.jawDistance), 'nearest')';
    
    % make the tongue movement into a binary format (exist in the frame or
    % not)
    videoTracking.tongueExist = zeros(size(videoTracking, 1), 1);
    ind = find(videoTracking.tongue_confidence > s.tongueConfidenceThreshold);
    videoTracking.tongueExist(ind) = 1;
        
    save(fullfile(sessionFolder, 'cameraAnalyzed.mat'), 'videoTracking');
    fprintf('cameraAnalyzed.mat saved in %s\n', sessionFolder)
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
    if any(spikeTemp.Food.times)
        spike.foodTriggerTimes = spikeTemp.Food.times;
        inds = find(diff(spike.foodTriggerTimes)<s.foodMinInterval);
        spike.foodTriggerTimes(inds) = [];
        spike.totalFoodNum = length(spike.foodTriggerTimes);
    end
    
    % get and save keyboard inputs
    spike.keyboardInput = char(spikeTemp.Keyboard.codes(:, 1));
    spike.keyboardTimes = spikeTemp.Keyboard.times;
    
    % get rid of false 'F' and cooresponding keyboard time according to
    % minimal food trial invervals
    FInds = nan(200, 1);
    FCount = 1;
    for i = 1:length(spike.keyboardInput)
        if strcmp(spike.keyboardInput(i), 'F')
            FInds(FCount, 1) = i;
            FCount = FCount + 1;
        end
    end
    
    if exist('inds', 'var')
        spike.keyboardTimes(FInds(inds)) = [];
        spike.keyboardInput(FInds(inds)) = [];
    end
    
    % get rid of false 'F' and cooresponding keybord time according to
    % manule keyboard input indication
    
    keyboardInds = nan(200, 1);
    FInds = nan(200, 1);
    FCount = 0;
    tempInd = 1;
    for i = 1:length(spike.keyboardInput)
        if strcmp(spike.keyboardInput(i), 'F') 
            FCount = FCount + 1;
            if strcmp(spike.keyboardInput(i+1), 'Z') || strcmp(spike.keyboardInput(i+1), '9')
                keyboardInds(tempInd, 1) = i;
                FInds(tempInd, 1) = FCount;
                tempInd = tempInd + 1;
            end
        end
    end
    
    if any(keyboardInds)
        spike.keyboardInput(rmmissing(keyboardInds)) = [];
        spike.keyboardTimes(rmmissing(keyboardInds)) = [];
        spike.foodTriggerTimes(rmmissing(FInds)) = [];
    end
    
    % only process the estim signal if this is an estim session
    if any(strcmp(fieldnames(spikeTemp), 'EStim'))
        spike.EstimTimes = spikeTemp.EStim.times;        
    end
    
    % only process the opto signal if this is an opto session
    if any(strcmp(fieldnames(spikeTemp), 'OptoTrain'))
        spike.optoTrainTimes = spikeTemp.OptoTrain.times;
    end
    
    save(fullfile(sessionFolder, 'spikeAnalyzed.mat'), 'spike', '-v7.3');
    fprintf('spikeAnalyzed.mat saved in %s\n', sessionFolder)    
end


end