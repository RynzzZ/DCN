function convertOpenEphysToSpikeTimes(session, varargin)

% This function syncs the time b/w the open ephys time system and the
% spike2 time system. It also converts the openEphys timestamps into spike2
% times and save the converted timestamps into the sessionEphysInfo.mat. 

% settings
s.forceAlignment = false;  % whether to run the alignement algorithm to find best matches between spike and ephys sync signals // if false, only runs the algorithm when there are different numbers of events in each channel
s.plot = true;

if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs

% initialization 
rootFolder = 'Z:\Qianyun\DCN\';
gitFolder = 'D:\DCN_Project\Github\DCN';


% getting ephys session info
sessionEphysInfo = getSessionEphysInfo(session);
ephysFolder = sessionEphysInfo.ephysFolder;

% get ephys event times
[channel, eventTimes, info] = load_open_ephys_data_faster(...
    fullfile(ephysFolder, 'all_channels.events'));  % load event data
eventChannel = unique(channel);  % assumes only one digital input is used!
ephysEventTimes = eventTimes(logical(info.eventId) & ismember(channel, eventChannel)); % only take rising edge of event channel

% get spike event times
spikeEvents = load(fullfile(rootFolder, 'Data', session, 'Spike.mat'), 'Keyboard', 'Food', 'AudioSyn');

spikeEventTimes = sort([spikeEvents.Food.times; spikeEvents.AudioSyn.times]);

% check if the length of event time in spike and in openephys are the same.
if length(spikeEventTimes) ~= length(ephysEventTimes)
    warning('ephysTimes do NOT match spikeTimes!');
end

% check if we need to force the alignment algorithm
runAlignment = length(ephysEventTimes)~=length(spikeEventTimes) || s.forceAlignment;

% simple linear mapping alignment
if ~runAlignment
    openEphysToSpikeMapping = polyfit(ephysEventTimes, spikeEventTimes, 1); % get linear mapping from open ephys to spike
    predictedEventSpikeTimes = polyval(openEphysToSpikeMapping, ephysEventTimes);
    if max(abs(predictedEventSpikeTimes - spikeEventTimes)) > .004
        fprintf('%s: WARNING! Same number of events in both channels, but linear mapping fails. Running alignment algorithm.\n', session)
        runAlignment = true;
    else
        fprintf('%s: OpenEphys Time successfully converted to Spike time. Max diff is %12.8f. \n', session, max(abs(predictedEventSpikeTimes - spikeEventTimes)));
    end
end

if s.plot
    figure; box off;
    plot(spikeEventTimes, ones(size(spikeEventTimes)), '.', 'MarkerSize', 10); hold on;
    plot(predictedEventSpikeTimes, zeros(size(predictedEventSpikeTimes)), '.', 'MarkerSize', 10);
    legend('spikeEventTimes', 'ephysEventTimesConvertedToSpikeTime');
    ylim([-1, 2]);
end

% save the converted ephys time and the conversion factors in the
% sessionEphysInfo.mat
sessionEphysInfo.ephysTimeConversionFactors = openEphysToSpikeMapping;
sessionEphysInfo.convertedEphysTimestamps = polyval(openEphysToSpikeMapping, sessionEphysInfo.timeStamps); 
save(fullfile(rootFolder, 'Data', session, 'sessionEphysInfo.mat'), 'sessionEphysInfo');
 
end