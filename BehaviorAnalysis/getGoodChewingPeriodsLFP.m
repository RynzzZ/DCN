function [chewingStruct] = getGoodChewingPeriodsLFP(LFPData, LFPTime, trialStartTimes, varargin)
%GETGOODCHEWINGPERIODSLFP Summary of this function goes here
%   Detailed explanation goes here

% settings
s.chewingSearchLength = 15; % sec, detemine the time range (how long) of the chewing search 
s.chewingStartTimeShift = 4; % sec, how many seconds after the trial start should we strat the chewing search
s.chewingChunkLength = 4; % sec
s.chewingSearchStepLength = 0.25; % sec
s.audioFs = 150000; % sampling rate for mic signal
s.fs = 30000; % sampling rate for ephys
s.cameraFs = 100; % sampling rate for the camera

if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs



chewingStartStopTimes = nan(length(trialStartTimes), 2);
chewingFreq = nan(length(trialStartTimes), 1);
chewingLFP = nan(length(trialStartTimes), s.chewingChunkLength*s.cameraFs);
chewingLFPTime = nan(length(trialStartTimes), s.chewingChunkLength*s.cameraFs);
chewingLFPInds = nan(length(trialStartTimes), s.chewingChunkLength*s.cameraFs);
chewingLFPLowPass = nan(length(trialStartTimes), s.chewingChunkLength*s.cameraFs);
LFPPhase = nan(length(trialStartTimes), s.chewingChunkLength*s.cameraFs);
selectedTrialNum = 1;


for i = 1:length(trialStartTimes)
    
    trialStartTime = trialStartTimes(i); 
    trialEndTime = trialStartTime + s.chewingStartTimeShift + s.chewingSearchLength;
    
    % locate chewing period
    chewingChunkStartTime = trialStartTime + s.chewingStartTimeShift;
    chewingChunkEndTime = chewingChunkStartTime + s.chewingChunkLength;
    fitResults = [];
    
    while(chewingChunkEndTime <= trialEndTime)
        chewingChunkEphysStartInd = find(LFPTime >= chewingChunkStartTime, 1, 'first');
        chewingChunkEphysEndInd = find(LFPTime <= chewingChunkEndTime, 1, 'last');
        
        chewingChunkTime = LFPTime(chewingChunkEphysStartInd:chewingChunkEphysEndInd);
        
        chewingChunkLFP = LFPData(chewingChunkEphysStartInd:chewingChunkEphysEndInd);
        chewingChunkLFP = smooth(chewingChunkLFP, 0.02);
        chewingChunkLowPass = lowpass(chewingChunkLFP, 5, 1/(chewingChunkTime(2) - chewingChunkTime(1)));
        chewingChunkLowSample = interp1(chewingChunkTime, chewingChunkLowPass, ...
            linspace(chewingChunkTime(1), chewingChunkTime(end), length(chewingChunkTime)/20));
        chewingChunkLowSample = smooth(chewingChunkLowSample, 0.009);
        
        chewingChunkLFPTime = LFPTime(chewingChunkEphysStartInd:chewingChunkEphysEndInd)' - LFPTime(chewingChunkEphysStartInd);
        chewingChunkTimeLowSample = interp1(1:length(chewingChunkLFPTime), chewingChunkLFPTime, ...
            linspace(1, length(chewingChunkLFPTime), length(chewingChunkLFPTime)/20));
        
        if size(chewingChunkLowSample, 1) ~= 1; chewingChunkLowSample = chewingChunkLowSample'; end
        if size(chewingChunkTimeLowSample, 1) ~= 1; chewingChunkTimeLowSample = chewingChunkTimeLowSample'; end
        
        fitResult = sineFit(chewingChunkTimeLowSample, chewingChunkLowSample, 0);
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
        chewingChunkEphysStartInd = find(LFPTime >= chewingStartStopTimes(selectedTrialNum, 1), 1, 'first');
        chewingChunkEphysEndInd = find(LFPTime <= chewingStartStopTimes(selectedTrialNum, 2), 1, 'last');
        
        chewingChunkLFP = LFPData(chewingChunkEphysStartInd:chewingChunkEphysEndInd);
        chewingChunkLFP = smooth(chewingChunkLFP, 0.02);
        chewingChunkLowPass = lowpass(chewingChunkLFP, 5, 1/(chewingChunkTime(2) - chewingChunkTime(1)));
        chewingChunkLowSample = interp1(chewingChunkTime, chewingChunkLowPass, ...
            linspace(chewingChunkTime(1), chewingChunkTime(end), length(chewingChunkTime)/20));
        chewingChunkLowSample = smooth(chewingChunkLowSample, 0.009);
        
        chewingChunkLFPTime = LFPTime(chewingChunkEphysStartInd:chewingChunkEphysEndInd)' - LFPTime(chewingChunkEphysStartInd);
        chewingChunkTimeLowSample = interp1(1:length(chewingChunkLFPTime), chewingChunkLFPTime, ...
            linspace(1, length(chewingChunkLFPTime), length(chewingChunkLFPTime)/20));
        
        
        % sanity check
        sineFit(chewingChunkTimeLowSample', chewingChunkLowSample); % sanity check
    
    
        % get LFP data & time for good chewing periods
        chewingLFP(selectedTrialNum, 1:length(chewingChunkLFP)) = chewingChunkLowSample;
        chewingLFPTime(selectedTrialNum, 1:length(chewingChunkLFP)) = chewingChunkTimeLowSample + LFPTime(chewingChunkEphysStartInd);

        
        
        % get jaw phase for good chewing periods
        y = chewingChunkLowSample;
        x = linspace(chewingStartStopTimes(selectedTrialNum, 1), chewingStartStopTimes(selectedTrialNum, 2), length(chewingChunkLowSample));
        % y2 = lowpass(y, 5, 1/(x(2) - x(1)));
        % chewingLFPLowPass(selectedTrialNum, 1:length(y2)) = y2;
        
        phaseTemp = hilbert(lowpass(y, 5, 1/(x(2) - x(1))));
        % phaseTemp = hilbert(y);
        phaseTemp = angle(phaseTemp);
        LFPPhase(selectedTrialNum, 1:length(phaseTemp)) = phaseTemp;
        
        % increment selectedTrialNum
        selectedTrialNum = selectedTrialNum + 1;
        
    end
    
      
end



% build the output variable - chewingStruct
chewingStruct.chewingLFPrmsv = chewingLFP;
chewingStruct.chewingLFPTime = chewingLFPTime;
chewingStruct.chewingLFPInds = chewingLFPInds;
chewingStruct.chewingStartStopTimes = chewingStartStopTimes;
chewingStruct.chewingLFPLowPass = chewingLFPLowPass;
chewingStruct.chewingLFPPhase = LFPPhase;




end

