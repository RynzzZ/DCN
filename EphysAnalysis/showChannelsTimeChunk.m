function showChannelsTimeChunk(session, startTimes, stopTimes, varargin) 

% given session name, shows traces of all channels on a probe as rows in a
% large column // does this for evenly spaced intervals over the duration
% of recording // use to assess drift!


% settings
s.yLims = [-500 500];             % (microvolts)
s.showSortedSpikes = false;       % whether to overlay sorted spikes
s.spkWindow = [-.5 1];            % (ms) pre and post spike time to plot for spike overlays
s.bestChannels = [];              % for each unit, the best channel for that unit // if provided, only overlay sorted spikes on this channel
s.selectedChannels = 'all';       % 'all' plots raw traces for all connected channles, otherwise specify the channel numbers to plot
s.figureName = fullfile('D:\DCN_Project\Figures', 'Ephys', [session '.png']);

if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs


% initializations
rootFolder = 'Z:\Qianyun\DCN\';
gitFolder = 'D:\DCN_Project\Github\';

load(fullfile(rootFolder, 'Data', session, 'sessionEphysInfo.mat'));
info = sessionEphysInfo;
load(fullfile('Z:\obstacleData\ephys\channelMaps\kilosort', [info.mapFile, '.mat']), ...
    'connected', 'channelNum_OpenEphys');

if length(startTimes) ~= length(stopTimes)
    error('ERROR: startTimes length does NOT match endTimes length!');
end

% properly determine number of channels to plot
if strcmp(s.selectedChannels, 'all')
    selectedChannels = channelNum_OpenEphys(connected);
else
    s.selectedChannels = sort(s.selectedChannels);
    selectedChannels = channelNum_OpenEphys(s.selectedChannels);
end

% properly determine the color scheme
if s.showSortedSpikes
    [spkInds, unit_ids] = getGoodSpkInds(session);
    colors = hsv(length(spkInds));
else
    colors = parula(selectedChannels);
end


% function to extract voltage from binary file
getVoltage = @(data, channel, inds) ...
    double(data.Data.Data(channel,inds))*info.bitVolts; % extract voltage from memmapfile, converting to voltage, highpassing, and only return specific channel

% load data
contFiles = dir(fullfile(sessionEphysInfo.ephysFolder, '*.continuous'));
data = memmapfile(fullfile(sessionEphysInfo.ephysFolder, [contFiles(1).name(1:end-12), 's.dat']), ...
    'Format', {'int16', [sessionEphysInfo.channelNum sessionEphysInfo.smps], 'Data'}, 'Writable', false);

% plot channels at specific times
traces = nan(length(selectedChannels), stopInd-startInd+1);
for i = 1:length(startTimes)    
    figure('color', 'white', 'Units', 'pixels', 'position', get(0,'ScreenSize'));
    
    startInd = find(info.convertedEphysTimestamps>startTimes(i),1,'first');
    stopInd = find(info.convertedEphysTimestamps<stopTimes(i),1,'last');
    traces = getVoltage(data, selectedChannels, startInd:stopInd);
    
    for j = 1:length(selectedChannels)
        offset = j*range(s.yLims);
        startTime = info.convertedEphysTimestamps(startInd);
        stopTime = info.convertedEphysTimestamps(stopInd);
        if s.showSortedSpikes; color=[.5 .5 .5]; else; color=colors(j,:); end
        plot(linspace(startTime, stopTime, stopInd-startInd+1), traces(j, :)+offset, 'Color', color);
        hold on; box off;
        if connected(i); color='black'; else; color='red'; end
        text(startTime, offset, ['(' num2str(channelNum_OpenEphys(j)) ')' num2str(selectedChannels(j))], 'Color', color);
    end
    
    
end

% % plot spikes on top of traces (if showSortedSpikes)
% if s.showSortedSpikes
%     lines = nan(1, length(unit_ids));
%     for i = 1:length(unit_ids)
%         for j = 1:length(timeStarts)
%             subplot(1,length(timeStarts),j); hold on
%             startInd = int64(find(info.timeStamps>timeStarts(j),1,'first'));
% 
%             % overlay spikes
%             traceSpkInds = spkInds{i}(info.timeStamps(spkInds{i})>=timeBinEdges(j) & info.timeStamps(spkInds{i})<timeBinEdges(j)+s.windowSize);
%             spkBins = repmat(spkWindowInds,length(traceSpkInds),1) + int64(traceSpkInds); % matrix where each row is inds for given spike
%             spkBins = reshape(spkBins',[],1);
%             spkBins = spkBins-startInd;
%             spkBins = spkBins(spkBins<size(traces,3) & spkBins>0);
% 
%             if ~isempty(s.bestChannels); channelsToShow=s.bestChannels(i); else; channelsToShow=1:info.channelNum; end
%             
%             for k = channelsToShow
%                 spikesTrace = nan(1,size(traces,3));
%                 spikesTrace(spkBins) = traces(k,j,spkBins);
%                 lines(i) = plot(timesSub, spikesTrace + offsets(k), ...
%                     'Color', [colors(i,:) .6], 'LineWidth', 1.5); % overlay spikes!
%             end
%         end
%     end
% end
% 
% 
% % fancify
% for i = 1:length(timeStarts)
%     subplot(1,length(timeStarts),i);
%     
%     set(gca, 'XLim', [0 timesSub(end)], ...
%         'YLim', [0 (info.channelNum+.5)*range(s.yLims)], ...
%         'Visible', 'off')
%     
%     text(range(timesSub)*.1, (info.channelNum+1)*(range(s.yLims)), ...
%         [num2str(round((timeStarts(i)-min(timeStarts))/60)) ' minutes'], ...
%         'FontWeight', 'bold')
% end

if s.showSortedSpikes; legend(lines, cellstr(num2str(unit_ids)), 'location', 'northeast'); end


% save that ish
saveas(gcf, s.figureName);
