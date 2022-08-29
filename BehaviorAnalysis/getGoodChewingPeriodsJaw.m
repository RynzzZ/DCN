function [chewingStruct] = getGoodChewingPeriodsJaw(jawDistance, frameTimestamps, trialStartTimes, varargin)
%GETGOODCHEWINGPERIODS Summary of this function goes here
%   Detailed explanation goes here

% settings
s.chewingSearchLength = 15; % sec, detemine the time range (how long) of the chewing search 
s.chewingStartTimeShift = 4; % sec, how many seconds after the trial start should we strat the chewing search
s.chewingChunkLength = 4; % sec
s.chewingSearchStepLength = 0.1; % sec
s.audioFs = 150000; % sampling rate for mic signal
s.fs = 30000; % sampling rate for ephys
s.cameraFs = 100; % sampling rate for the camera

if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs



chewingStartStopTimes = nan(length(trialStartTimes), 2);
chewingFreq = nan(length(trialStartTimes), 1);
chewingJawDistance = nan(length(trialStartTimes), s.chewingChunkLength*s.cameraFs);
chewingJawTime = nan(length(trialStartTimes), s.chewingChunkLength*s.cameraFs);
chewingJawInds = nan(length(trialStartTimes), s.chewingChunkLength*s.cameraFs);
chewingJawDistanceLowPass = nan(length(trialStartTimes), s.chewingChunkLength*s.cameraFs);
jawPhase = nan(length(trialStartTimes), s.chewingChunkLength*s.cameraFs);
selectedTrialNum = 1;


for i = 1:length(trialStartTimes)
    
    trialStartTime = trialStartTimes(i); 
    trialEndTime = trialStartTime + s.chewingStartTimeShift + s.chewingSearchLength;
    
    % locate chewing period
    chewingChunkStartTime = trialStartTime + s.chewingStartTimeShift;
    chewingChunkEndTime = chewingChunkStartTime + s.chewingChunkLength;
    fitResults = [];
    
    while(chewingChunkEndTime <= trialEndTime)
        chewingChunkVideoStartInd = find(frameTimestamps >= chewingChunkStartTime, 1, 'first');
        chewingChunkVideoEndInd = find(frameTimestamps <= chewingChunkEndTime, 1, 'last');
        
        chewingChunkJaw = jawDistance(chewingChunkVideoStartInd:chewingChunkVideoEndInd);
        x = frameTimestamps(chewingChunkVideoStartInd:chewingChunkVideoEndInd)' - frameTimestamps(chewingChunkVideoStartInd);
        fitResult = sineFit(x, chewingChunkJaw', 0);
        fitResults = [fitResults; fitResult];
        
        chewingChunkStartTime = chewingChunkStartTime + s.chewingSearchStepLength;
        chewingChunkEndTime = chewingChunkEndTime + s.chewingSearchStepLength;
    end
    
    fitMSE = fitResults(:, 5);
    tempInd = find(fitMSE == min(fitMSE(fitMSE>0)));
    minFitMSE = min(fitMSE(fitMSE>0));
    disp(['minFitMSE = ', num2str(minFitMSE)]);
    
    if fitResults(tempInd, 2) > 1.5 && fitResults(tempInd, 3) > 0.03
        disp(['trial ' num2str(i) ' selected'])
        chewingFreq(selectedTrialNum,1) = fitResults(tempInd, 3);
        chewingStartStopTimes(selectedTrialNum, 1) = trialStartTime + s.chewingStartTimeShift + s.chewingSearchStepLength*(tempInd - 1);
        chewingStartStopTimes(selectedTrialNum, 2) = chewingStartStopTimes(selectedTrialNum, 1) + s.chewingChunkLength;
        chewingChunkVideoStartInd = find(frameTimestamps >= chewingStartStopTimes(selectedTrialNum, 1), 1, 'first');
        chewingChunkVideoEndInd = find(frameTimestamps <= chewingStartStopTimes(selectedTrialNum, 2), 1, 'last');
        chewingChunkJaw = jawDistance(chewingChunkVideoStartInd:chewingChunkVideoEndInd);
        % sanity check
        sineFit(1:length(chewingChunkJaw), chewingChunkJaw'); % sanity check
        
        % get jaw distance & time for good chewing periods
        chewingJawDistance(selectedTrialNum, 1:length(chewingChunkJaw)) = chewingChunkJaw;
        chewingJawTime(selectedTrialNum, 1:length(chewingChunkJaw)) = frameTimestamps(chewingChunkVideoStartInd:chewingChunkVideoEndInd);
        chewingJawInds(selectedTrialNum, 1:length(chewingChunkJaw)) = [chewingChunkVideoStartInd:chewingChunkVideoEndInd];
        
        
        % get jaw phase for good chewing periods
        y = chewingChunkJaw;
        x = linspace(chewingStartStopTimes(selectedTrialNum, 1), chewingStartStopTimes(selectedTrialNum, 2), length(chewingChunkJaw));
        y2 = lowpass(y, 5, 1/(x(2) - x(1)));
        chewingJawDistanceLowPass(selectedTrialNum, 1:length(y2)) = y2;
        
        jawPhaseTemp = hilbert(lowpass(y, 5, 1/(x(2) - x(1))));
        jawPhaseTemp = angle(jawPhaseTemp);
        jawPhase(selectedTrialNum, 1:length(jawPhaseTemp)) = jawPhaseTemp;
        
        % increment selectedTrialNum
        selectedTrialNum = selectedTrialNum + 1;
        
    end
    
      
end



% build the output variable - chewingStruct
chewingStruct.chewingJawDistance = chewingJawDistance(1:selectedTrialNum-1, :);
chewingStruct.chewingJawTime = chewingJawTime(1:selectedTrialNum-1, :);
chewingStruct.chewingJawCameraInds = chewingJawInds(1:selectedTrialNum-1, :);
chewingStruct.chewingStartStopTimes = chewingStartStopTimes(1:selectedTrialNum-1, :);
chewingStruct.chewingJawDistanceLowPass = chewingJawDistanceLowPass(1:selectedTrialNum-1, :);
chewingStruct.chewingjawPhase = jawPhase(1:selectedTrialNum-1, :);




end

