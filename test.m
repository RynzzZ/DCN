%% Save the processed calibration data

% run the Sound_Analysis.m first, detrend the signal, normalize the signal,
% then plot the signal spectrogram.

fileName = 'MicCalibration_20210926.mat';
folderPath = 'D:\DCN_Project\Data\Calibration_20210926\';
fullPath = fullfile([folderPath, fileName]);
save(fullPath, 'F', 'T', 'STPS');


%% make microphone/speaker calibration curve

% load data
load(fullPath);

% settings 
minFreq = 3000;
maxFreq = 48000;
numFreq = 30;

% calculate the frequency for all the tones played.
freqLog = calculateFrequencies(minFreq, maxFreq, numFreq);

% get the power for all the frequencies played
STPS_freqInds = zeros(numFreq, 1);
powers = zeros(numFreq, 1);
for i = 1:length(freqLog)
    STPS_freqInds(i) = find(abs(F - freqLog(i)) == min(abs(F - freqLog(i)))); 
    powers(i) = max(STPS(STPS_freqInds(i), :));
end

% plot the calibration curve
figure('Color', 'white', 'position', get(0,'ScreenSize')); clf; hold on;

plot(freqLog/1000, powers, '-', 'LineWidth', 2);
ylim([-70, 0]);

xlabel('Frequency (kHz)');
ylabel('power (dBV^2)');
title('Mic Calibration Curve, 3KHz to 48KHz, 0.5V');

box off;

% % plot Mic1 and Mic2 together
% figure('Color', 'white', 'position', get(0,'ScreenSize')); clf; hold on;
% 
% plot(frequencyMic1, peakPowerMic1, '-', 'LineWidth', 2);
% xlabel('Frequency (kHz)');
% ylabel('power (dBV^2)');
% 
% ylim([-70, 0]);
% box off; hold on;
% 
% plot(frequencyMic2, peakPowerMic2, '-', 'LineWidth', 2);
% 
% legend('Mic1', 'Mic2');
% 
% title('Mic Test, Calibration Curve');
% 




%%
