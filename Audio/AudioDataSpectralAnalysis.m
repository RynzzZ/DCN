%% Load the audio file and plot spectrograms
close all
clear all
clc

fileFolder = 'D:\DCN_Project\AudioFiles\MouseUSVs\';
fileName = 'SocialContextComparison_B6D2F1J_2.WAV';
filePath = fullfile(fileFolder, fileName);

[audioData, Fs] = audioread(filePath);  
audioData = audioData';
n = length(audioData); 
t = 0:1/Fs:(n-1)/Fs;             % Time 

n1 = pow2(nextpow2(n));
audioDataFFT = fft(audioData, n1);
f = (0:n1-1)*(Fs/n1);            % frequency range
f = f/1000;
power = abs(audioDataFFT).^2/n;  % power of the DFT

% AudioData in Time Domain
figure('Color', 'white', 'position', get(0,'ScreenSize'));
subplot(2, 1, 1);
plot(t/60, audioData);
xlabel('time (min)');
ylabel('Amplitude');
title('Signals in Time Domain');
xlim([0, max(t/60)]);

% AudioData Spectral Analysis
subplot(2, 1, 2);
plot(f(1:floor(n1/2)),power(1:floor(n1/2)))
xlabel('Frequency (kHz)');
ylabel('Power');
title('Spectral Analysis');
xlim([0,100]);
xticks(0:5:120);

% spectrogram
figure('Color', 'white', 'position', get(0,'ScreenSize'));
% Window duration (in seconds):
dur = 0.5;
% Spectrogram settings (in samples):
winSize = round(Fs*dur);
overlap = round(winSize/2);
fftSize = winSize;
% Plot the spectrogram:
spectrogram(audioData,winSize,overlap,fftSize,Fs,'yaxis');
title('Spectrogram - Social Context B6D2F1 mice');
yticks(0:5:120);

%% Locate Start and Stop Position
hold on;

startXPosition = 1.85; 
stopXPosition = 1.91; 

line([startXPosition, startXPosition],[0, 120]);
line([stopXPosition, stopXPosition], [0, 120]);

%% Create Audio Clips
fileName = 'chirp_17Kto20K_100msRepeat_2min_2.wav';

startTime = startXPosition*60;
endTime = stopXPosition*60;

audioClip = test(startTime*Fs:endTime*Fs);
audiowrite(fileName, audioData2, Fs, 'BitsPerSample',24);

%% Make waveform .txt files

TXTFileName = [fileName(1:end-4), '.txt'];
fid = fopen(TXTFileName, 'w')
fprintf(fid, '%12.8f\r\n', audioData2)

%% test for highpass filtering

test = highpass(audioData, 3000, 250000);
% test = highpass(audioClip, 3000, 250000);

n = length(test); 
t = 0:1/Fs:(n-1)/Fs;             % Time 

n1 = pow2(nextpow2(n));
audioDataFFT = fft(test, n1);
f = (0:n1-1)*(Fs/n1);            % frequency range
f = f/1000;
power = abs(audioDataFFT).^2/n;  % power of the DFT

% AudioData in Time Domain
figure('Color', 'white', 'position', get(0,'ScreenSize'));
subplot(2, 1, 1);
plot(t/60, test);
xlabel('time (min)');
ylabel('Amplitude');
title('Signals in Time Domain');
xlim([0, max(t/60)]);

% AudioData Spectral Analysis
subplot(2, 1, 2);
plot(f(1:floor(n1/2)),power(1:floor(n1/2)))
xlabel('Frequency (kHz)');
ylabel('Power');
title('Spectral Analysis');
xlim([0,100]);
xticks(0:5:120);

% spectrogram
figure('Color', 'white', 'position', get(0,'ScreenSize'));
% Window duration (in seconds):
dur = 0.5;
% Spectrogram settings (in samples):
winSize = round(Fs*dur);
overlap = round(winSize/2);
fftSize = winSize;
% Plot the spectrogram:
spectrogram(test,winSize,overlap,fftSize,Fs,'yaxis');
title('Spectrogram - Social Context B6D2F1 mice');
yticks(0:5:120);






