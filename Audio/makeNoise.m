function makeNoise(fs, Length, ampPeakPeak, highpassFreq)

% make the noise signal
noise = rand(1, fs*Length) * ampPeakPeak;
noise = noise - (ampPeakPeak/2);

%     dlmwrite(['C:\Users\sawtell\Desktop\github\obstacleRig\matlab\noiseFs' num2str(fs)...
%         'Length' num2str(length) '.txt'], noise', 'delimiter', '\t');

if exist('highpassFreq', 'var')
    noise = highpass(noise, highpassFreq, fs);
end

% make the file name
if exist('highpassFreq', 'var')
    fid = fopen(['noiseFs' num2str(fs) 'Length' num2str(Length) 'Peak' num2str(ampPeakPeak/2) 'V' ...
        'highpass' num2str(highpassFreq/1000) 'K.txt'], 'w');
else
    fid = fopen(['noiseFs' num2str(fs) 'Length' num2str(Length) 'Peak' num2str(ampPeakPeak/2) 'V.txt'], 'w');
end

% save the file
fprintf(fid, '%12.8f\r\n', noise);
fclose(fid);

end