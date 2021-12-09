function formatEphysData(session, varargin)
% prepares spikes for matlab-land by getting spike times wrt spike clock,
% and maybe by computing instantaneous firing rates as well...


% settings
s.spkRateFs = 200;         % sampling frequency of instantaneous firing rate
s.kernelRise = .001;       % (s) rise for double exponential kernel
s.kernelFall = .02;        % (s) fall for double exponential kernel
s.kernelSig = .02;         % (s) if a gaussian kernel is used
s.kernel = 'doubleExp';    % 'gauss', or 'doubleExp'
s.forceAlignment = false;  % whether to run the alignement algorithm to find best matches between spike and ephys sync signals // if false, only runs the algorithm when there are different numbers of events in each channel
s.plot = true;             % whether to show plot when forcing alignment... 
s.outputFileName = 'neuralData.mat'; 


% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs

% some paths
rootFolder = 'Z:\Qianyun\DCN\';
gitFolder = 'D:\DCN_Project\Github\DCN\';
sessionFolder = fullfile(rootFolder, 'Data', session);
GITDIR = 'D:\Cerebellum_DCN_Project\Github';
outputFilePath = fullfile(sessionFolder, s.outputFileName);

addpath(genpath(fullfile(GITDIR, 'analysis-tools')));
addpath(genpath(fullfile(GITDIR, 'npy-matlab')));

load(fullfile(sessionFolder, 'sessionEphysInfo.mat'));
ephysInfo = sessionEphysInfo;
ephysFolder = ephysInfo.ephysFolder;


% get spike times for good units
[spkInds, unit_ids] = getGoodSpkInds(session);
[bestChannels, ~] = getBestChannels(session, 'returnPhysicalLayout', true);
cellData = readtable(fullfile(sessionFolder, 'cellData.csv'));
if ~isequal(cellData.unit_id(:), unit_ids(:)); disp('WARNING! cellData.csv unit_ids do not match those in ephysFolder'); keyboard; end
if ~isequal(length(unit_ids), length(bestChannels)); disp('WARNING! bestChannels do not match unit_ids'); keyboard; end


% restrict to good units
goodBins = cellData.include==1;
unit_ids = unit_ids(goodBins);
bestChannels = bestChannels(goodBins);
spkInds = spkInds(goodBins);
cellData = cellData(goodBins,:);

spkTimes = cell(1, length(unit_ids));
for i = 1:length(unit_ids)
    spkTimes{i} = ephysInfo.convertedEphysTimestamps(spkInds{i});
end

% convert to instantaneous firing rate
minTime = min(cellfun(@(x) x(1), spkTimes)); % latest spike across all neurons
maxTime = max(cellfun(@(x) x(end), spkTimes)); % latest spike across all neurons
[~, timeStamps] = getFiringRate(spkTimes{1}, 'tLims', [minTime maxTime], 'fs', s.spkRateFs, ...
    'kernel', s.kernel, 'kernelRise', s.kernelRise, 'kernelFall', s.kernelFall, 'sig', s.kernelSig);
spkRates = nan(length(spkTimes), length(timeStamps));

for i = 1:length(spkTimes)
    
    [spkRates(i,:), timeStamps] = getFiringRate(spkTimes{i}, 'tLims', [minTime maxTime], 'fs', s.spkRateFs, ...
        'kernel', s.kernel, 'kernelRise', s.kernelRise, 'kernelFall', s.kernelFall, 'sig', s.kernelSig);    
    
    % get min and max time for cell
    cellMinTime = polyval(ephysInfo.ephysTimeConversionFactors, cellData.timeStart(i)*60);
    if strcmp(cellData.timeEnd(i), 'max')
        cellMaxTime = timeStamps(end);
    else
        if iscell(cellData.timeEnd(i)); timeEnd = str2double(cellData.timeEnd(i)); else; timeEnd = cellData.timeEnd(i); end
        cellMaxTime = polyval(ephysInfo.ephysTimeConversionFactors, timeEnd*60);
    end
    
    % remove spikes that are out of min and max times
    spkRates(i, timeStamps<cellMinTime | timeStamps>cellMaxTime) = nan;
    spkTimes{i} = spkTimes{i}(spkTimes{i}>cellMinTime & spkTimes{i}<cellMaxTime);
end

settings = s;
openEphysToSpikeMapping = ephysInfo.ephysTimeConversionFactors;
save(outputFilePath, 'spkRates', 'spkTimes', 'timeStamps', 'unit_ids', 'bestChannels', 'openEphysToSpikeMapping', 'settings')




