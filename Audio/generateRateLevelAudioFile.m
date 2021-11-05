function Sound = generateRateLevelAudioFile(BF, startDBSPL, stopDBSPL, varargin)

% settings
s.isBF = true;  % Ture for generating rate level audio file for BF, False for broadbandnoise.
s.highpass = 4000;  % only for BBN
s.DBSPLStepLength = 1;  % step length for DBSPL level
s.Tdur = 200/1000;  % tone/bbn duration (unit: sec)
s.Tint = 200/1000;  % interval duration
s.Tram = 10/1000;  % ramps duration
s.Fs = 150000;  % sampling frequency

if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs

Fs = s.Fs;
Ts = 1/Fs;                                                                 % Sampling Interval (s)
T1 = 0:Ts:(Ts*(Fs-1/Fs))*s.Tdur;                                           % Time duration for each pure tone (eg, 200ms)
S2 = zeros(1, s.Tint*Fs);                                                  % Silent Interval
T3 = s.Tram*Fs;                                                            % Time for ramp up and ramp down

DBSPLLevel = startDBSPL:s.DBSPLStepLength:stopDBSPL;

Sound = zeros(1, length(T1)*length(DBSPLLevel) + length(S2)*(length(DBSPLLevel)-1));
j = 1;
for i = 1:length(DBSPLLevel)
    Amplitude = convertDBSPL2Voltage(DBSPLLevel(i));
    if s.isBF
        % This part generates the rate level sound file for BF tones.
        temp = sin(2*pi*BF*T1)*Amplitude;
        Sound(j : j + length(T1) - 1) = temp.*[linspace(0, 1, T3), ones(1, length(T1) - T3*2), linspace(1, 0, T3)];
        j = j + length(T1);
        if i < length(DBSPLLevel)
            Sound(j : j + length(S2) -1) = 0;
            j = j + length(S2);
        end
    else
        % This part generates the rate level sound file for BBN.
        noiseTemp = rand(1, Fs*s.Tdur) * Amplitude;
        noiseTemp = noiseTemp - (Amplitude/2);
        noiseTemp = noiseTemp*2;
        
        if s.highpass ~= 0
            noiseTemp = highpass(noiseTemp, s.highpass, Fs);
        end
      
        Sound(j : j + length(T1) - 1) = noiseTemp;
        j = j + length(T1);
        if i < length(DBSPLLevel)
            Sound(j : j + length(S2) -1) = 0;
            j = j + length(S2);
        end     
    end
end

end