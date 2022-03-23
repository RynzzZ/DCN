function FR = processAuditoryData(keyboardTime, audioSig, audioSigTime, ephysData, ephysTime, toneNumber, s, tag)


tempInd = int32(find(audioSigTime >= keyboardTime, 1, 'first'));

spikeStartInds = nan(toneNumber, 1);
ephysStartInds = nan(toneNumber, 1);
spikeStopInds = nan(toneNumber, 1);
ephysStopInds = nan(toneNumber, 1);
FR = nan(toneNumber, 1);

for j = 1:toneNumber
    if mod(j, 10) == 1
        disp([tag, ': ', num2str(round(j*100/toneNumber)), '%']);
    end
    
    spikeStartInds(j) = int32(find(audioSig(tempInd:end) >= 0.003, 1, 'first') + tempInd);
    startTime = audioSigTime(spikeStartInds(j));
    ephysStartInds(j) = find(ephysTime >= startTime, 1, 'first');
    
    spikeStopInds(j) = int32(spikeStartInds(j) + s.toneDuration*s.audioFs);
    stopTime = audioSigTime(spikeStopInds(j));
    ephysStopInds(j) = find(ephysTime <= stopTime, 1, 'last');
    
    
    chunkEphysData = ephysData(ephysStartInds(j):ephysStopInds(j));
    [~, locs, ncrs] = crossdet(chunkEphysData, s.threshold, s.thresholdType);
    FR(j) = ncrs/(stopTime - startTime);
    
    tempInd = int32(spikeStopInds(j) + s.intervalDuration/2*s.audioFs);
    
    % quality check
    if mod(j, 8) == 0 && ~s.suppressFig
        figure;
        x = linspace(startTime, stopTime, ephysStopInds(j) - ephysStartInds(j) + 1);
        plot(x, chunkEphysData); axis tight; box off; hold on;
        y = ones(size(locs))*s.threshold;
        scatter(x(locs), y, '.', 'r');
        title([num2str(j), '/', num2str(tonNumber)]);
    end
end

if size(FR, 2) == 1; FR = FR'; end

end