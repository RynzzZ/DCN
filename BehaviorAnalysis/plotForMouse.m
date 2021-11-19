%% extract jaw movement and overlay it with DCN LFP

session = '20211111_002';

% initializations
rootFolder = 'Z:\Qianyun\DCN\';
gitFolder = 'D:\DCN_Project\Github\DCN\';
sessionFolder = fullfile(rootFolder, 'Data', session);

% load analyzed mat files
load(fullfile(sessionFolder, 'spikeAnalyzed.mat'));
load(fullfile(sessionFolder, 'cameraAnalyzed.mat'));
load(fullfile(sessionFolder, 'sessionEphysInfo.mat'));

% get food trigger times
randomSampleNum = 5;
randomFoodTriggerTimes = sort(randsample(spike.foodTriggerTimes, randomSampleNum));

for i = 1:randomSampleNum
    
    figure('Color', 'white', 'position', get(0,'ScreenSize')); clf; hold on;
    
    startTime = randomFoodTriggerTimes(i) - 0.3;
    stopTime = randomFoodTriggerTimes(i) + 10;
    
    % plot mic signal
    startInd = find(spike.micSignalTimes > startTime, 1, 'first');
    stopInd = find(spike.micSignalTimes < stopTime, 1, 'last');
    ind = startInd:stopInd;
    
    plot(spike.micSignalTimes(ind) - randomFoodTriggerTimes(i), spike.micSignal(ind)/max(spike.micSignal(ind))); 
    hold on
    
    % plot jaw movement
    yoffset = 2;
    
    startInd = find(videoTracking.frameTimestamps > startTime, 1, 'first');
    stopInd = find(videoTracking.frameTimestamps < stopTime, 1, 'last');
    ind = startInd:stopInd;
        
    plot(videoTracking.frameTimestamps(ind) - randomFoodTriggerTimes(i),...
        videoTracking.jawDistance(ind)/max(videoTracking.jawDistance(ind)) + yoffset);
    
    
    % plot tongue movement
    yoffset = 3;
    
    startInd = find(videoTracking.frameTimestamps > startTime, 1, 'first');
    stopInd = find(videoTracking.frameTimestamps < stopTime, 1, 'last');
    ind = startInd:stopInd;
    
    plot(videoTracking.frameTimestamps(ind) - randomFoodTriggerTimes(i), videoTracking.tongueExist(ind) + yoffset);
    plot(videoTracking.frameTimestamps(ind) - randomFoodTriggerTimes(i), videoTracking.tongue_confidence(ind) + yoffset);
    
    % plot DCN LFP
    yoffset = 1;
    
    LFPChannel = 30;
    trace = showChannelsTimeChunk(session, startTime, stopTime, 'selectedChannels', LFPChannel, 'dataOnly', true);
    plot(linspace(startTime, stopTime, length(trace)) - randomFoodTriggerTimes(i), ((trace - min(trace))/(max(trace) - min(trace))) + yoffset);
    
    xlim([startTime , stopTime] - randomFoodTriggerTimes(i));
    ylim([-1, 4]);
    
    xlabel('time (sec, 0 = food dispenser trigger)');
    legend('mic signal', 'jaw movement', 'tongue', 'tongue confidence', 'DCN LFP');
    ax = gca;
    ax.YAxis.Visible = 'off'; % remove y-axis
end






