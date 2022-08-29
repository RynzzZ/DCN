function makeUnitVid(session, fileName, varargin)

% creates videos of mouse feeding behavior with the raw neural trace plotted in real time
% along with the sound of the neuron // tell it the session, the unit_id
% of the neuron to plot, and a matrix with the start and stop times
% (columns) of trials to plot (rows) 


% TO DO:
% gate sound of spikes



% settings
s.trialsToShow = 4; % how many random trials to show if specific trial numbers are not indicated
s.vidType = 'showFeedingTrials'; % 
s.specificTrials = []; % pick specific trials to show. 
s.specificTimeWindows = []; % put in specific start times and end times for generating videos.
s.feedingTrialTimeWindow = 10; % the duration of each food trial.
s.timeBuffer = [0.5, 0]; % how many seconds before and after reward delivery to show. 

s.ephysChannel = 11;
s.ephysThresh = -200;

s.contrastLims = [.1 .9]; % pixels at these proportional values are mapped to 0 and 255
s.playbackSpeed = 0.5;
s.voltageWindow = .8;
s.audioGain = 15;
s.yLims = [-500 500];
s.lowPassFreq = 6000; % 6000 // set to false to turn off lowpass

s.compressVideo = false;

s.spkScatterColor = [1 1 0];
s.lineColors = [.8 .4 1];

% reassign settings contained in opts
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initializations
rootFolder = 'Z:\Qianyun\DCN\';
gitFolder = 'D:\DCN_Project\Github\DCN\';
sessionFolder = fullfile(rootFolder, 'Data', session);

% set up video readers / writer
disp('initializing...')
% vidName = fullfile(getenv('OBSDATADIR'), 'sessions', session, 'run_originalDimensions.mp4');
vidName = fullfile(sessionFolder, 'behaviorVideo00.avi');
if ~exist(vidName, 'file'); error('behavior vid NOT found in the session folder!'); end
vid = VideoReader(vidName);
% determine frame dimensions
frameDim = [vid.Height, vid.Width];

initialFs = vid.FrameRate;
vidWriter = vision.VideoFileWriter(fileName, ...
    'AudioInputPort', true, ...
    'FrameRate', round(initialFs*s.playbackSpeed));
% vidWriter.VideoCompressor = 'MJPEG Compressor';


% load spikeAnalyzed.mat
if ~exist(fullfile(sessionFolder, 'spikeAnalyzed.mat'), 'file')
    error('spikeAnalyzed.mat NOT exist! Run analyzeSession.m first!');
else
    disp('loading spikeAnalyzed.mat...')
    load(fullfile(sessionFolder, 'spikeAnalyzed.mat'), 'spike');
end

% load cameraAnalyzed.mat
if ~exist(fullfile(sessionFolder, 'cameraAnalyzed.mat'), 'file')
    error('cameraAnalyzed.mat NOT exist! Run analyzeSession.m first!');
else
    disp('loading cameraAnalyzed.mat...')
    load(fullfile(sessionFolder, 'cameraAnalyzed.mat'), 'videoTracking');
    if ~any(strcmp('jawDistance', videoTracking.Properties.VariableNames))
        error('JawDistance does NOT exist in videoTracking! Run analyzeSession.m first!');
    end
    frameTimeStamps = videoTracking.frameTimestamps;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get neural data
if ~exist(fullfile(sessionFolder, 'sessionEphysInfo.mat'), 'file')
    error('sessioinEphysInfo.mat NOT exist!');
else
    % load sessionEphysInfo
    disp('loading ephys data...')
    load(fullfile(rootFolder, 'Data', session, 'sessionEphysInfo.mat'), 'sessionEphysInfo');
    mapFile = sessionEphysInfo.mapFile;
    load(fullfile('Z:\obstacleData\ephys\channelMaps\kilosort', [mapFile, '.mat']), 'channelNum_OpenEphys');
    
    ephysInfo = sessionEphysInfo;
    ephysTime = sessionEphysInfo.convertedEphysTimestamps;
    timeStampsMapped = ephysTime;
    audioSmpsPerFrame = round((1/initialFs) * sessionEphysInfo.fs);
   
    % function to extract voltage from binary file
    getVoltage = @(data) ...
        double(data)*sessionEphysInfo.bitVolts; % extract voltage from memmapfile, converting to votlage, highpassing, and only return specific channel

    
    % load raw ephys data and turn it into voltage data
    contFiles = dir(fullfile(sessionEphysInfo.ephysFolder, '*.continuous'));
    data = memmapfile(fullfile(sessionEphysInfo.ephysFolder, [contFiles(1).name(1:end-12), 's.dat']), ...
        'Format', {'int16', [sessionEphysInfo.channelNum sessionEphysInfo.smps], 'Data'}, 'Writable', false);
    ephysVoltage = getVoltage(data.Data.Data(channelNum_OpenEphys(s.ephysChannel), :));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up figure
disp('setting up fig for vid...');
fig = figure('color', [0 0 0], 'position', [500, 500, frameDim(2), frameDim(1)], 'menubar', 'none');
traceLength = s.voltageWindow*ephysInfo.fs;

axes('position', [0 .2 1 .8], 'CLim', [0 255]); colormap gray
im = image(uint8(zeros(frameDim)), 'CDataMapping', 'scaled'); hold on;
set(gca, 'Visible', 'off')

plotAxis = axes('position', [0 0 1 .2], 'Color', 'black');
tracePlot = plot(plotAxis, 1:traceLength, nan(1,traceLength), 'color', 'white'); set(gca, 'color', 'black'); hold on
set(gca, 'Visible', 'off', 'YLimMode', 'manual', 'YLim', s.yLims)

foodTriggerLine = line(plotAxis, [0 0], s.yLims, 'linewidth', 2, 'color', s.lineColors);
foodTriggerText = text(plotAxis, 0, s.yLims(2), 'food trigger', 'Color', 'white');
currentTimeLine = line(plotAxis, [0, 0], s.yLims, 'linewidth', 2, 'Color', 'r');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get timeEpochs for trials to display in the vid
% current setup only supports 'showObsEvents' and 'showRewardEvents'
switch s.vidType
    case 'showFeedingTrials'
        timeEpochs = cat(2, spike.foodTriggerTimes - s.timeBuffer(1), spike.foodTriggerTimes + s.feedingTrialTimeWindow + + s.timeBuffer(2));
        minTime = ephysTime(1);
        maxTime = ephysTime(end);
        validTrials = find(timeEpochs(:,1)>minTime & timeEpochs(:,2)<maxTime); % make sure trials aren't too long
        
        if length(s.specificTrials) == 0  % specific trials are not indicated 
            trialsToShow = validTrials(round(linspace(1,length(validTrials),s.trialsToShow)));
            timeEpochs = timeEpochs(trialsToShow, :);
        else
            trialNum = s.specificTrials;
            trialsToShow = trialNum(ismember(trialNum, validTrials));
            timeEpochs = timeEpochs(trialsToShow, :);
        end
        
%     case 'showRewardEvents'
%         timeEpochs = cat(2, rewardTimes - s.timeBuffer(1), rewardTimes + s.timeBuffer(2));
%         minTime = timeStamps(find(~isnan(spkRates(unitInd,:)),1,'first'));
%         maxTime = timeStamps(find(~isnan(spkRates(unitInd,:)),1,'last'));
%         validTrials = find(timeEpochs(:,1)>minTime & timeEpochs(:,2)<maxTime); 
%         % make sure trials aren't too long
%         if length(s.specificObsTrials) == 0 & length(s.specificRewardTrials) == 0  % specific trials are not indicated 
%             trialsToShow = validTrials(round(linspace(1,length(validTrials),s.trialsToShow)));
%             timeEpochs = timeEpochs(trialsToShow, :);
%         else
%             trialNum = s.specificRewardTrials;
%             trialsToShow = trialNum(ismember(trialNum, validTrials));
%             timeEpochs = timeEpochs(trialsToShow, :);
%         end
        
    case 'showSpecificTimeWindows'
        timeEpochs = [];
        timeWindows = s.specificTimeWindows;

        minTime = ephysTime(1);
        maxTime = ephysTime(end);
        for i = 1:size(timeWindows, 1)
            if timeWindows(i, 1) < minTime || timeWindows(i, 2) > maxTime
                disp(['WARNING: time windows ' num2str(i) ' you selected have exceeded the min/max unit time!!']);
            else
                timeEpochs = [timeEpochs; timeWindows(i, :)];       
            end            
        end
        trialsToShow = 1:size(timeEpochs, 1);        
end

          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create video
fprintf('\nwriting video... \n')

for i = 1:length(trialsToShow)
    fprintf('%d/%d\n', i, length(trialsToShow));
    % get trial frame indices
    trialInds = find(frameTimeStamps>timeEpochs(i,1) & ...
        frameTimeStamps<timeEpochs(i,2));
    
    % get voltage for entire trial
    voltageBins = timeStampsMapped>timeEpochs(i,1)-1 & ...
                  timeStampsMapped<timeEpochs(i,2)+1; % add and subtract 1 as a buffer
    timeStampsSub = timeStampsMapped(voltageBins);
    voltage = ephysVoltage(voltageBins);
    
    [~, spikeInds, ~] = crossdet(voltage, s.ephysThresh, 'thresholdDown');
    chunkEphysTime = ephysTime(voltageBins);
    trialSpkTimes = chunkEphysTime(spikeInds);
    audioTemp = false(1, length(voltage));
    audioTemp(spikeInds) = true;
    
    
    % update obstacle and whisker contact lines
    if strcmp(s.vidType, 'showFeedingTrials')
        foodTriggerString = 'dispenser triggered';
        updateTextAndLine(foodTriggerText, foodTriggerLine, spike.foodTriggerTimes, foodTriggerString);
    end
   
    
    % get frames for trials
    for j = trialInds'
        if j == trialInds(end)
            break;
        end
        frame = rgb2gray(read(vid, j));
        % add trial number onto frames
        position = [10, 10];
        if strcmp(s.vidType, 'showFeedingTrials')
            textString = ['trial ', num2str(trialsToShow(i))];
            RGB = insertText(frame, position, textString, 'TextColor','white');
            frame = rgb2gray(RGB);
        end
                
        % update frame
        set(im, 'CData', frame);
        
        % get voltage
        traceStartInd = find(timeStampsSub>(frameTimeStamps(j)-s.voltageWindow/2), 1, 'first');
        traceInds = traceStartInd:traceStartInd+traceLength/2-1;
        trace = voltage(traceInds);
        times = timeStampsSub(traceInds)';
        set(tracePlot, 'xdata', times, 'ydata', trace);
        set(gca, 'xlim', [times(1) times(end)])
        
        % scatter the detected spike!
        scatter(trialSpkTimes, ...  
            repmat(s.yLims(1),1,length(trialSpkTimes)) + range(s.yLims)*.1, ... % y values
            10, s.spkScatterColor, 'filled');
        currentTime = (times(1)+times(end))/2;
        updateTextAndLine('', currentTimeLine, currentTime);
        
        % get audio
        audioStartInd = find(timeStampsSub>=currentTime, 1, 'first');
        audio = audioTemp(audioStartInd:audioStartInd+audioSmpsPerFrame-1)';
        audio = audio * s.audioGain;
        
%         % add fade in/out for first and last sample
%         if i==1 && j==trialInds(1)
%             audio = int16(double(audio) .* linspace(0,1,length(audio))'); % fade in on first sample
%         elseif i==length(trialsToShow) && j==trialInds(end)
%             audio = int16(double(audio) .* linspace(1,0,length(audio))'); % fade out on last sample
%         end
   
        % write to file   
        frame = getframe(fig);
        vidWriter(frame.cdata, audio);
    end
end
release(vidWriter);
close(fig)


% compress video
if s.compressVideo
    baseDir = fileparts(fileName);
    [~,~] = system([fileName(1) ': & ffmpeg -i ' fileName ' -vcodec mpeg4 -vb 10M -y ' baseDir '\temp.avi']); % run ffmpeg to compress file
    %delete(fileName)
    %movefile([baseDir '\temp.avi'], fileName)
end

disp(' all done!')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% update position of lines and text marking events of interest
function updateTextAndLine(text, line, eventTimes, textString)
    eventInd = find(eventTimes>=timeEpochs(i,1)-s.voltageWindow & ...
                    eventTimes<=timeEpochs(i,2),1,'first');
    if ~isempty(eventInd)
        eventTime = eventTimes(eventInd);
        set(text, 'Position', [eventTime+s.voltageWindow*.01 s.yLims(2)])
        set(line, 'XData', [eventTime eventTime])
    end
    if exist('textString', 'var'); set(text, 'String', textString); end
end

end