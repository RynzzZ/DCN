function frameTimes = convertCameraToSpikeTimes(session, varargin)

% settings
s.ttlGap = [0.995, 1.005]; % the last frame is collected after a gap // look for a gap that is within these limits // this shojld match the setting in the Arduino script that controls the camTrigger.

if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs


% initializations
rootFolder = 'Z:\Qianyun\DCN\';
gitFolder = 'D:\DCN_Project\Github\DCN\';
sessionFolder = fullfile(rootFolder, 'Data', session);

% load both Bonsai and Spike camera timestamp data
if ~exist(fullfile(rootFolder, 'Data', session, 'cameraTimestamps.csv'), 'file') || ~exist(fullfile(rootFolder, 'Data', session, 'Spike.mat'), 'file')
    error('cameraTimestamps.csv or Spike.mat is missing! Can not proceed!');
else
    fprintf('%s: getting frame time stamps\n', session);
    load(fullfile(sessionFolder, 'Spike.mat'), 'CamTrig');
    spikeCamTimes = CamTrig.times;
    
    camMetadata = dlmread(fullfile(sessionFolder, 'cameraTimestamps.csv')); % columns: bonsai timestamps, point grey counter, point grey timestamps (uninterpretted)
    frameCounts = camMetadata(:,1);
    bonsaiCamTimes = timeStampDecoderFLIR(camMetadata(:,2));
end


% find ttlGaps in spike ttls camera metadata
ttlGap = s.ttlGap;
ttlGaps = find(diff(spikeCamTimes)>ttlGap(1) & diff(spikeCamTimes)<ttlGap(2));
camGaps = find(diff(bonsaiCamTimes)>ttlGap(1) & diff(bonsaiCamTimes)<ttlGap(2));

if length(ttlGaps)==1 && length(camGaps)==1  % if same number of ttl gaps, just use the first one
    ttlGap = 1;
    camGap = 1;
else
    ttlGapDurations = diff(spikeCamTimes(ttlGaps));
    camGapDurations = diff(bonsaiCamTimes(camGaps));
    
    if length(ttlGaps)==length(camGaps)  % if same number of gaps, use the most similar gap (to avoid using a gap where a frame may be lost at the beginning or the end)
        [~, minInd] = min(abs(ttlGapDurations-camGapDurations));
        ttlGap = minInd;
        camGap = minInd;
    else  % otherwise try to find the best match (note that the following approach will fail for evenly spaced ttl gaps or if two gaps happen to be of similar duration) // would be smarter to loop across all entries, only keeping those with closely match durations
        fprintf('%s: WARNING! %i gaps in camera and % gaps in spike!\n', session, length(camGapDurations), length(ttlGapDurations));
        [inds, diffs] = knnsearch(ttlGapDurations, camGapDurations);
        [~, minInd] = min(diffs);
        ttlGap = inds(minInd);
        camGap = minInd;
    end
end

camCountsShifted = frameCounts - frameCounts(camGaps(camGap)) + ttlGaps(ttlGap);

% if camera was started before spike, account for missing ttls at beginning
if any(camCountsShifted<1)
    beginningNans = nan(sum(camCountsShifted<1), 1);
    camCountsShifted = camCountsShifted(camCountsShifted>0);
    fprintf('%s: WARNING! %i camera frame(s) detected before spike recording turned on!\n', session, length(beginningNans));
end

% if camera was stopped after spike, account for extra ttls at end
if any(camCountsShifted>length(spikeCamTimes))
    endNans = nan(sum(camCountsShifted>length(spikeCamTimes)), 1);
    camCountsShifted = camCountsShifted(camCountsShifted<=length(spikeCamTimes));
    fprintf('%s: WARNING! %i camera frame(s) detected after last spike ttl!\n', session, length(endNans));
end

% get spike times for each frame
frameTimes = spikeCamTimes(camCountsShifted);
if exist('beginningNans', 'var'); frameTimes = [beginningNans; frameTimes]; end
if exist('endNans', 'var'); frameTimes = [frameTimes; endNans]; end

% display number of missed frames
missedFrames = length(spikeCamTimes)-length(bonsaiCamTimes);
if exist('beginningNans', 'var'); missedFrames = missedFrames + length(beginningNans); end
if exist('endNans', 'var'); missedFrames = missedFrames + length(endNans); end
fprintf('%s: %i missed frame(s), %i at the very beginning, %i elsewhere\n', ...
    session, missedFrames, camCountsShifted(1)-1, missedFrames - (camCountsShifted(1)-1));

end