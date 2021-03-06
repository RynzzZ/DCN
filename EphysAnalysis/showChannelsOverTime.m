function showChannelsOverTime(session, varargin) 

% given session name, shows traces of all channels on a probe as rows in a
% large column // does this for evenly spaced intervals over the duration
% of recording // use to assess drift!


% settings
s.yLims = [-500 500];             % (microvolts)
s.timeBinNum = 5;                 % number of time periods to show
s.showSortedSpikes = false;       % whether to overlay sorted spikes
s.windowSize = 2.0/s.timeBinNum;  % (s) width of each time window
s.spkWindow = [-.5 1];            % (ms) pre and post spike time to plot for spike overlays
s.bestChannels = [];              % for each unit, the best channel for that unit // if provided, only overlay sorted spikes on this channel
s.figureName = fullfile(getenv('OBSDATADIR'), 'figures', 'ephys', 'drift', [session '.png']);


% initializations
rootFolder = 'Z:\Qianyun\DCN\';
gitFolder = 'D:\DCN_Project\Github\';

if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs

load(fullfile(rootFolder, 'Data', session, 'sessionEphhysInfo.mat'));
info = sessionEphysInfo;
load(fullfile('Z:\obstacleData\ephys\channelMaps\kilosort', [info.mapFile, '.mat']), ...
    'connected', 'channelNum_OpenEphys')


% transfer the openephys channel order into its physical location order on the probe
% in order to determine the spacing offset for every channel for the drift plots.
[~, sortedInds] = sort(channelNum_OpenEphys);

if s.showSortedSpikes
    [spkInds, unit_ids] = getGoodSpkInds(session);
    spkWindowInds = int64((s.spkWindow(1)/1000*info.fs) : (s.spkWindow(2)/1000*info.fs));
    colors = hsv(length(spkInds));
else
    colors = parula(info.channelNum);
end


% function to extract voltage from binary file
getVoltage = @(data, channel, inds) ...
    double(data.Data.Data(channel,inds))*info.bitVolts; % extract voltage from memmapfile, converting to voltage, highpassing, and only return specific channel

% load data
contFiles = dir(fullfile(sessionEphysInfo.ephysFolder, '*.continuous'));
data = memmapfile(fullfile(sessionEphysInfo.ephysFolder, [contFiles(1).name(1:end-12), 's.dat']), ...
    'Format', {'int16', [sessionEphysInfo.channelNum sessionEphysInfo.smps], 'Data'}, 'Writable', false);


% plot all channels at dft time intervals
[~, name, ~] = fileparts(s.figureName);
figure('color', 'white', 'Units', 'pixels', 'position', get(0,'ScreenSize'));

timeBinEdges = linspace(min(info.timeStamps), max(info.timeStamps), s.timeBinNum+1);
% timeStarts = timeBinEdges(1:end-1);
timeStart = 143.6;
timeEnd = 147.7;
% traces = nan(info.channelNum, length(timeStarts), s.windowSize*info.fs);
traces = nan(info.channelNum, int32((timeEnd - timeStart)*info.fs));

% collect traces
for j = 1:length(timeStarts)
    startInd = find(info.timeStamps>timeStarts(j),1,'first');
    inds = startInd:startInd+s.windowSize*info.fs-1;
    traces(:,j,:) = getVoltage(data, 1:info.channelNum, inds);
end

startInd = find(info.convertedEphysTimestamps > timeStart, 1, 'first');
endInd = find(info.convertedEphysTimestamps > timeEnd, 1, 'first');
inds = startInd:endInd - 1;

traces = getVoltage(data, 1:info.channelNum, inds);

% plot traces
timesSub = linspace(0, s.windowSize, s.windowSize*info.fs); % xlim for plotting the traces
offsets = (info.channelNum - sortedInds + 1)*range(s.yLims); % set offsets for every channel based on its physical location on the probe
for i = 1:info.channelNum
    for j = 1:length(timeStart)
        subplot(1,length(timeStarts),j); hold on
        if s.showSortedSpikes; color=[.5 .5 .5]; else; color=colors(i,:); end        
        plot(timesSub, squeeze(traces(i,j,:)) + offsets(i), 'Color', color);
        if connected(i); color='black'; else; color='red'; end
        text(timesSub(1), offsets(i), ['(' num2str(i) ')' num2str(sortedInds(i))], 'Color', color)
    end
end

% plot traces
timesSub = linspace(0, s.windowSize, s.windowSize*info.fs); % xlim for plotting the traces
offsets = (info.channelNum - sortedInds + 1)*range(s.yLims); % set offsets for every channel based on its physical location on the probe
for i = 1:info.channelNum

        hold on
        if s.showSortedSpikes; color=[.5 .5 .5]; else; color=colors(i,:); end        
        plot(info.convertedEphysTimestamps(inds), traces(i, :) + offsets(i), 'Color', color);
        if connected(i); color='black'; else; color='red'; end
        text(info.convertedEphysTimestamps(inds(1)), offsets(i), ['(' num2str(i) ')' num2str(sortedInds(i))], 'Color', color)

end

% plot spikes on top of traces (if showSortedSpikes)
if s.showSortedSpikes
    lines = nan(1, length(unit_ids));
    for i = 1:length(unit_ids)
        for j = 1:length(timeStarts)
            subplot(1,length(timeStarts),j); hold on
            startInd = int64(find(info.timeStamps>timeStarts(j),1,'first'));

            % overlay spikes
            traceSpkInds = spkInds{i}(info.timeStamps(spkInds{i})>=timeBinEdges(j) & info.timeStamps(spkInds{i})<timeBinEdges(j)+s.windowSize);
            spkBins = repmat(spkWindowInds,length(traceSpkInds),1) + int64(traceSpkInds); % matrix where each row is inds for given spike
            spkBins = reshape(spkBins',[],1);
            spkBins = spkBins-startInd;
            spkBins = spkBins(spkBins<size(traces,3) & spkBins>0);

            if ~isempty(s.bestChannels); channelsToShow=s.bestChannels(i); else; channelsToShow=1:info.channelNum; end
            
            for k = channelsToShow
                spikesTrace = nan(1,size(traces,3));
                spikesTrace(spkBins) = traces(k,j,spkBins);
                lines(i) = plot(timesSub, spikesTrace + offsets(k), ...
                    'Color', [colors(i,:) .6], 'LineWidth', 1.5); % overlay spikes!
            end
        end
    end
end


% fancify
for i = 1:length(timeStarts)
    subplot(1,length(timeStarts),i);
    
    set(gca, 'XLim', [0 timesSub(end)], ...
        'YLim', [0 (info.channelNum+.5)*range(s.yLims)], ...
        'Visible', 'off')
    
    text(range(timesSub)*.1, (info.channelNum+1)*(range(s.yLims)), ...
        [num2str(round((timeStarts(i)-min(timeStarts))/60)) ' minutes'], ...
        'FontWeight', 'bold')
end

if s.showSortedSpikes; legend(lines, cellstr(num2str(unit_ids)), 'location', 'northeast'); end


% save that ish
saveas(gcf, s.figureName);
