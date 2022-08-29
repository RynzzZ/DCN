function [FR, optoFR, intervalFR, optoIntervalFR] = processAuditoryData(keyboardTime, audioSig, audioSigTime, ephysData, ephysTime, toneNumber, s, tag)


tempInd = int32(find(audioSigTime >= keyboardTime, 1, 'first'));

FR = nan(toneNumber, 1);
optoFR = nan(toneNumber, 1);
intervalFR = nan(toneNumber, 1);
optoIntervalFR = nan(toneNumber, 1);
optoFlag = false;

for j = 1:toneNumber
    
    if mod(j, 10) == 1
        disp([tag, ': ', num2str(round(j*100/toneNumber)), '%']);
    end
    
    % get tone start and stop times&inds
    spikeStartInds = int32(find(audioSig(tempInd:end) >= 0.003, 1, 'first') + tempInd);
    startTime = audioSigTime(spikeStartInds) + s.toneEphysLag;
    ephysStartInds = find(ephysTime >= startTime, 1, 'first');
    
    spikeStopInds = int32(spikeStartInds + s.toneDuration*s.audioFs);
    stopTime = audioSigTime(spikeStopInds) + s.toneEphysLag;
    ephysStopInds = find(ephysTime <= stopTime, 1, 'last');
    
    % get interval start and stop times&inds
    intervalEphysStartInds = ephysStopInds + 1;
    intervalStartTime = ephysTime(intervalEphysStartInds);
    intervalStopTime = stopTime + s.intervalDuration;
    intervalEphysStopInds = find(ephysTime <= intervalStopTime, 1, 'last');
 
    % decide if this tone is with opto
    if ~isempty(s.optoTrainTimes)
        if sum(s.optoTrainTimes >= startTime & s.optoTrainTimes <= stopTime) >= s.toneDuration*s.optoTrainFreq*0.8
            optoFlag = true;
        else
            optoFlag = false;
        end
    end
    
    % calculate firing rate for the tone
    chunkEphysData = ephysData(ephysStartInds:ephysStopInds);
    [~, ~, ncrs] = crossdet(chunkEphysData, s.threshold, s.thresholdType);
    if ~optoFlag
        FR(j) = ncrs/(stopTime - startTime);
    else
        optoFR(j) = ncrs/(stopTime - startTime);
    end    
    
    % calculate firing rate for the interval
    chunkEphysData = ephysData(intervalEphysStartInds:intervalEphysStopInds);
    [~, ~, ncrs] = crossdet(chunkEphysData, s.threshold, s.thresholdType);
    if ~optoFlag
        intervalFR(j) = ncrs/(intervalStopTime - intervalStartTime);
    else
        optoIntervalFR(j) = ncrs/(intervalStopTime - intervalStartTime);
    end
    
    tempInd = int32(spikeStopInds + s.intervalDuration/2*s.audioFs);
    
    % quality check
%     if mod(j, 8) == 0 && ~s.suppressFig
%         figure;
%         x = linspace(startTime, stopTime, ephysStopInds(j) - ephysStartInds(j) + 1);
%         plot(x, chunkEphysData); axis tight; box off; hold on;
%         y = ones(size(locs))*s.threshold;
%         scatter(x(locs), y, '.', 'r');
%         % title([num2str(j), '/', num2str(toneNumber)]);
%     end
    
%     if j == toneNumber && ~s.suppressFig
%         figure;
%         toneStartInd = int32(find(audioSigTime >= keyboardTime, 1, 'first'));
%         toneStopInd = spikeStopInds(end);
%         toneStartTime = audioSigTime(toneStartInd);
%         toneStopTime = audioSigTime(toneStopInd);
%         ephysStartInd = find(ephysTime >= toneStartTime, 1, 'first');
%         ephysStopInd = find(ephysTime <= toneStopTime, 1, 'last');
%         
%         subplot(2, 1, 1)
%         x = linspace(toneStartTime, toneStopTime, ephysStopInd - ephysStartInd + 1);
%         y = ephysData(1, ephysStartInd:ephysStopInd);
%         plot(x, y); 
%         axis tight; hold on;
%         for i = 1:toneNumber
%             plot([ephysTime(ephysStartInds(i)), ephysTime(ephysStartInds(i))], [min(y), max(y)], '-c');
%             plot([ephysTime(ephysStopInds(i)), ephysTime(ephysStopInds(i))], [min(y), max(y)], '-r');
%         end
%         title('Ephys Channel Data');
%         
%         subplot(2, 1, 2)
%         x = linspace(toneStartTime, toneStopTime, toneStopInd - toneStartInd + 1);
%         y = audioSig(toneStartInd:toneStopInd);
%         plot(x, y);
%         axis tight; hold on;
%         for i = 1:toneNumber
%             plot([audioSigTime(spikeStartInds(i)), audioSigTime(spikeStartInds(i))], [min(y), max(y)], '-c');
%             plot([audioSigTime(spikeStopInds(i)), audioSigTime(spikeStopInds(i))], [min(y), max(y)], '-r');
%         end
%         title('Audio Signal Data');
%         
%     end

end

if size(FR, 2) == 1; FR = FR'; end
if size(optoFR, 2) == 1; optoFR = optoFR'; end
if size(intervalFR, 2) == 1; intervalFR = intervalFR'; end
if size(optoIntervalFR, 2) == 1; optoIntervalFR = optoIntervalFR'; end 

end