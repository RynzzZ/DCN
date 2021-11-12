clear all;

minFreq = 7000;
maxFreq = 14000;
numFreq = 40;

% minFreq = 6000;
% maxFreq = 12000;
% numFreq = 11;
 
% minFreq = 12000;
% maxFreq = 24000;
% numFreq = 11;

% minFreq = 24000;
% maxFreq = 48000;
% numFreq = 11;

FreqLinear = [ minFreq:(maxFreq - minFreq)/(numFreq - 1):maxFreq ];        % Tone Frequency (linear)
FreqLog = logspace(log10(minFreq), log10(maxFreq), numFreq);               % Tone Frequency (log)

Tdur = 200/1000;
Tint = 50/1000;
Tram = 10/1000;
Amplitude = convertDBSPL2Voltage(20);                                      % Amplitude: volt.

Fs = 150*1000;                                                             % Default Sampling Frequency (Hz)
Ts = 1/Fs;                                                                 % Sampling Interval (s)
T1 = 0:Ts:(Ts*(Fs-1/Fs))*Tdur;                                             % Time duration for each pure tone (eg, 200ms)
S2 = zeros(1, Tint*Fs);                                                    % Silent Interval
T3 = Tram*Fs;                                                              % Time for ramp up and ramp down


% make the tone for playing
Tone = zeros(1, length(T1)*numFreq + length(S2)*(numFreq-1));
j = 1;
for i = 1:numFreq
    temp = sin(2*pi*FreqLog(i)*T1)*Amplitude;
    Tone(j : j + length(T1) - 1) = temp.*[linspace(0, 1, T3), ones(1, length(T1) - T3*2), linspace(1, 0, T3)];
    j = j + length(T1);
    disp(num2str(j))
    if i < numFreq
        Tone(j : j + length(S2) -1) = 0;
        j = j + length(S2);
        disp(num2str(j))
    end
end

% soundsc(Tone,Fs);                     % Play Sound

%% make wave.txt files

fileFolder = 'D:\DCN\AudioFiles\PureTones\BF_ToneSequence\';
fileName = 'DCN_ToneSequence_7Kto14K_20DBSPL_200msDuration_10msRamps_50msInterval.WAV';
filePath = fullfile(fileFolder, fileName);
audiowrite(filePath, Tone, Fs);


TXTFileName = [fileName(1:end-4), '.txt'];
FullPath = fullfile(fileFolder, TXTFileName);
fid = fopen(FullPath, 'w');
fprintf(fid, '%2.4f\r\n', Tone);
fclose(fid);



%% generate bbn noise

fs = 150000;
Length = 0.002; % sec
amp = [0.025, 0.05, 0.075, 0.1:0.1:1, 1.2:0.2:3]*2;
intervalDuration = 0.2; % sec

% make the noise signal
randomAmp = randsample(amp, length(amp));
noise = [];
noiseTemp = [];

for i = 1:length(amp)
    noiseTemp = rand(1, fs*Length) * randomAmp(i);
    noiseTemp = noiseTemp - (randomAmp(i)/2);
    
    if exist('highpassFreq', 'var')
        noise = highpass(noise, highpassFreq, fs);
    end
    
    noise = [noise, noiseTemp];
    
    interval = zeros(1, fs*intervalDuration);
    
    noise = [noise, interval];    
end

%% save bbn noise file

fileFolder = 'D:\DCN_Project\AudioFiles\Noise\';
fileName = 'DCN_NoiseRandomSequence6_2msDuration200msInterval_Amp0.025Vto3V.WAV';
filePath = fullfile(fileFolder, fileName);
audiowrite(filePath, noise, fs);


TXTFileName = [fileName(1:end-4), '.txt'];
FullPath = fullfile(fileFolder, TXTFileName);
fid = fopen(FullPath, 'w');
fprintf(fid, '%6.4f\r\n', noise);
fclose(fid);

MATFileName = [fileName(1:end-4), '.mat'];
FullPath = fullfile(fileFolder, MATFileName);
save(FullPath, 'noise', 'randomAmp');

%% generate and save rate level function audio files

BF = 10000;
startDBSPL = 10;
stopDBSPL = 70;

% this part generates rate level function file for BF
Sound = generateRateLevelAudioFile(BF, startDBSPL, stopDBSPL);
fileFolder = 'D:\DCN\AudioFiles\RateLevelFunction\';
fileName = 'BF10KHz_10to70DBSPL_1DBSPLInterval_RateLevelFunction.txt';

FullPath = fullfile(fileFolder, fileName);
fid = fopen(FullPath, 'w');
fprintf(fid, '%6.4f\r\n', Sound);
fclose(fid)

% this part generates rate level function file for BBN
BBN = generateRateLevelAudioFile(BF, startDBSPL, stopDBSPL, 'isBF', false);

fileFolder = 'D:\DCN\AudioFiles\RateLevelFunction\';
fileName = 'BBNNEW_10to70DBSPL_1DBSPLInterval_RateLevelFunction.txt';

FullPath = fullfile(fileFolder, fileName);
fid = fopen(FullPath, 'w');
fprintf(fid, '%6.4f\r\n', BBN);
fclose(fid);


%% generate and save BBN chewing sound mimic

% settings
samplingFreq = 150000; % Hz
chewingFreq = 4; % Hz 

Tdur = 50/1000; % duration for each chew

DBSPLAmplitude = 70;
voltageAmplitude = convertDBSPL2Voltage(DBSPLAmplitude); % loudness of the sound

repeat = 10; % how many times the chewing mimic reprats
highpassFreq = 4000; % for highpassing the bbn, unit: Hz

% make brief bbn sample for mimicing the chews
noiseTemp = rand(1, samplingFreq*Tdur) * voltageAmplitude;
noiseTemp = noiseTemp - (voltageAmplitude/2);
noiseTemp = noiseTemp*2;
noiseTemp = highpass(noiseTemp, highpassFreq, samplingFreq);

% make audio files for the chewing mimic
chewingMimic = zeros(1, repeat/chewingFreq*samplingFreq);
Ttotal = repeat/chewingFreq; % total time, unit: sec
startTimes = linspace(0, Ttotal, repeat+1); % unit: sec
startInds = startTimes*samplingFreq+1;
for i = 1:repeat
    startInd = startInds(i);
    chewingMimic(startInd : startInd + length(noiseTemp) - 1) = noiseTemp;    
end


fileFolder = 'D:\DCN\AudioFiles\Mimic_BBN_Noise\';
fileName = 'ChewMimics_50msDuration_10Repeats_70DBSPL.txt';

FullPath = fullfile(fileFolder, fileName);
fid = fopen(FullPath, 'w');
fprintf(fid, '%6.4f\r\n', chewingMimic);
fclose(fid);

%% generate BBN chewing sound mimic 2 

% settings
samplingFreq = 150000; % Hz
chewingFreq = 4; % Hz
DBSPLAmplitude = 70; % DBSPL
Tdur = 2/1000; % duration for fast BBN;
Tdur2 = 10/1000; % duration for each 

voltageAmplitude = convertDBSPL2Voltage(DBSPLAmplitude); % loudness of the sound

mimicRepeat = 10; % how many times the chewing mimic reprats
fastBBNRepeat = 3; % how many times the brief BBN (2ms) repeats
highpassFreq = 4000; % for highpassing the bbn, unit: Hz

% make brief BBN (2ms)
noiseTemp = rand(1, samplingFreq*Tdur) * voltageAmplitude;
noiseTemp = noiseTemp - (voltageAmplitude/2);
noiseTemp = noiseTemp*2;
noiseTemp = highpass(noiseTemp, highpassFreq, samplingFreq);

chewingMimic = zeros(1, mimicRepeat/chewingFreq*samplingFreq);
Ttotal = mimicRepeat/chewingFreq; % total time, unit: sec
startTimes = linspace(0, Ttotal, mimicRepeat+1); % unit: sec
startInds = startTimes*samplingFreq+1;
for i = 1:mimicRepeat
    % make BBN chewing sound mimic
    for j = 1:fastBBNRepeat
        startInd = startInds(i) + (j-1)*samplingFreq*Tdur2;
        chewingMimic(startInd : startInd + length(noiseTemp) - 1) = noiseTemp;
    end
end


fileFolder = 'D:\DCN\AudioFiles\Mimic_BBN_Noise\';
fileName = 'ChewMimics_3X2msFastBBN_10Repeats_70DBSPL.txt';

FullPath = fullfile(fileFolder, fileName);
fid = fopen(FullPath, 'w');
fprintf(fid, '%6.4f\r\n', chewingMimic);
fclose(fid);




















