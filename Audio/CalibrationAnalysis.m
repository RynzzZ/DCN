filePath = 'D:\DCN_Project\Data\Calibration_20210926\';
fileName = '20210926_002_AudioPlay_48K3K_0.5V.mat';

load(fullfile(filePath, fileName));

micData = Mic.values;
audioData = AudioSig.values;
sampleFreq = 150000;
freqCount = 30;
freqDuration = 0.2; % second
intervalDuration = 0.1; % second

%% detrend the signal
% organize question dialog menu about the detrending

x = micData;

quest = 'Do you want to detrend the signal?';
dlgtitle = 'Detrending';
btn1 = 'Yes, detrend the signal';
btn2 = 'No, do not detrend the signal';
defbtn = btn1;
answer = questdlg(quest, dlgtitle, btn1, btn2, defbtn);

% detrend the signal 
switch answer
    case btn1
    % detrend the signal    
    x = detrend(x);                             
    case btn2
    % do not detrend the signal
end

%% normalize the signal
% organize question dialog menu about the normalization
quest = 'What type of normalization do you want?';
dlgtitle = 'Normalization';
btn1 = 'Normalize the signal to unity peak';
btn2 = 'Normalize the signal to unity RMS-value';
btn3 = 'Do not normalize the signal';
defbtn = btn1;
answer = questdlg(quest, dlgtitle, btn1, btn2, btn3, defbtn);

% normalize the signal
switch answer
    case btn1
    % normalize to unity peak
    x = x/max(abs(x));
    case btn2
    % normalize to unity RMS-value
    x = x/std(x); 
    case btn3
    % do not normalize the signal
end

%% find the period of the audio signal and calculate the power

figure('Color', 'white', 'position', get(0,'ScreenSize')); 
h1 = subplot(2, 1, 1);
micData = x;
plot(micData); hold on; box off
h2 = subplot(2, 1, 2);
plot(audioData); hold on; box off

threshold = max(audioData*0.5);
tempInds = zeros(1, freqCount);
onsetInds = zeros(1, freqCount);
offsetInds = zeros(1, freqCount);
audioGain = zeros(1, freqCount);
inds = zeros(freqCount, 2);
temp = 1; % the start ind in every loop.
for i = 1:freqCount
    tempInds(i) = find(audioData(temp:end) > threshold, 1, 'first');
    
    plot(h1, [tempInds(i)+temp, tempInds(i)+temp], [-max(micData) - 0.1, max(micData) + 0.1], '-r');
    plot(h1, [tempInds(i)+temp+sampleFreq*freqDuration, tempInds(i)+temp+sampleFreq*freqDuration], [-max(micData) - 0.1, max(micData) + 0.1], '-c');
    ylim(h1, [-max(micData) - 0.1, max(micData) + 0.1]);
    
    
    plot(h2, [tempInds(i)+temp, tempInds(i)+temp], [-1, 1], '-r');
    plot(h2, [tempInds(i)+temp+sampleFreq, tempInds(i)+temp+sampleFreq], [-1, 1], '-c');
    
    % save the audio onset and offset index
    onsetInds(i) = tempInds(i)+temp;
    offsetInds(i) = tempInds(i) + temp + sampleFreq*freqDuration;
    
    % calculate the gain
    audioGain(i) = range(micData(onsetInds(i):offsetInds(i)))/2;
    
    % get to the next frequency
    inds(i, 1) = tempInds(i) + temp;
    inds(i, 2) = tempInds(i) + temp + sampleFreq*(freqDuration + 0.5*intervalDuration);
    temp = temp + tempInds(i) + 150000 + 50000;
end

% plot the gain for each frequency on the micData
xpos = onsetInds + sampleFreq*(freqDuration/2);
plot(h1, xpos, audioGain, 'o-', 'MarkerSize', 5, 'LineWidth', 1.5);

% audioNormalizeFactor = min(audioGain)./audioGain;
audioNormalizeFactor = mean(audioGain)./audioGain;




