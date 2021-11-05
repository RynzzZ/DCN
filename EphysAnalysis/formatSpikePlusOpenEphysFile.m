function formatSpikePlusOpenEphysFile(session, varargin)

% settings

if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs

% initialization 
rootFolder = 'Z:\Qianyun\DCN\';
gitFolder = 'D:\DCN_Project\Github\';
sessionFolder = fullfile(rootFolder, 'Data', session);

% add path to CED code

setenv('CEDS64ML', 'D:\DCN_Project\Github\Spike2\Spike2_MATLAB_Interface\CEDS64ML');

cedpath = getenv('CEDS64ML');
addpath(cedpath);
% load ceds64int.dll
CEDS64LoadLib( cedpath );

% Open the .smr spike file in the session folder and copy it to a new .smr
% file
fhand1 = CEDS64Open( fullfile(sessionFolder, [session(1:8) '_000.smr']) );
if (fhand1 <= 0);  CEDS64ErrorMessage(fhand1); unloadlibrary ceds64int; return; end
maxchans = CEDS64MaxChan( fhand1 );

spikeFile = load(fullfile(rootFolder, 'Data', session, 'Spike.mat'));
spikeChanNumber = length(fieldnames(spikeFile));
load(fullfile(rootFolder, 'Data', session, 'sessionEphysInfo.mat'));
ephysChanNumber = sessionEphysInfo.channelNum;


% Create the new .smr file - Spike + OpenEphys daya together
fhand2 = CEDS64Create( fullfile(sessionFolder, 'SpikePlusOpenEphys.smr'), maxchans + ephysChanNumber, 2 );
if (fhand2 <= 0);  CEDS64ErrorMessage(fhand2); unloadlibrary ceds64int; return; end
% Set timebase in new file
timebase = CEDS64TimeBase( fhand1 );
if timebase < 0, CEDS64ErrorMessage(timebase), return; end
CEDS64TimeBase( fhand2, timebase );
maxTimeTicks = CEDS64MaxTime( fhand1 )+2;

% set file comment
[ iOk, sCom1 ] = CEDS64FileComment( fhand1, 1 );
CEDS64FileComment( fhand2, 1, sCom1 );
[ iOk, sCom2 ] = CEDS64FileComment( fhand1, 2 );
CEDS64FileComment( fhand2, 2, sCom2 );
[ iOk, sCom3 ] = CEDS64FileComment( fhand1, 3 );
CEDS64FileComment( fhand2, 3, sCom3 );
[ iOk, sCom4 ] = CEDS64FileComment( fhand1, 4 );
CEDS64FileComment( fhand2, 4, sCom4 );
[ iOk, sCom5 ] = CEDS64FileComment( fhand1, 5 );
CEDS64FileComment( fhand2, 5, sCom5 );

% the maximum number of items to copy from each channel
maxpoints = maxTimeTicks;

% loop through all channels in fhand1 and copy them to fhand2
for m = 1:maxchans
    
    disp(num2str(m)); % display which channel is processing now
    
    chan = CEDS64ChanType( fhand1, m );
    chan2 = CEDS64GetFreeChan( fhand2 );
    
    disp([' ch', num2str(chan2), ' in fhand2.']);
    
    if (chan > 0) % is there a channel m?
        chandiv = CEDS64ChanDiv( fhand1, m );
        rate = CEDS64IdealRate( fhand1, m );
    end
    
    switch(chan)    
        case 0 % there is no channel with this number
        case 1 % ADC channel
            [shortread, shortvals, shorttime] = CEDS64ReadWaveS( fhand1, m, maxpoints, 0 );
            CEDS64SetWaveChan( fhand2, chan2, chandiv, 1, rate );
            CEDS64WriteWave( fhand2, chan2, shortvals, shorttime );
        case 2 % Event Fall
            [evread, evtimes] = CEDS64ReadEvents( fhand1, m, maxpoints, 0 );
            CEDS64SetEventChan( fhand2, chan2, rate, 2 );
            CEDS64WriteEvents( fhand2, chan2, evtimes );
        case 3 % Event Rise
            [evread, evtimes] = CEDS64ReadEvents( fhand1, m, maxpoints, 0 );
            CEDS64SetEventChan( fhand2, chan2, rate, 3 );
            CEDS64WriteEvents( fhand2, chan2, evtimes );
        case 4 % Event Both.
            [levread, levtimes, levinit] = CEDS64ReadLevels(fhand1, m, maxpoints, 0);
            CEDS64SetLevelChan( fhand2, chan2, rate );
            CEDS64SetInitLevel( fhand2, chan2, levinit );
            CEDS64WriteLevels( fhand2, chan2, levtimes );
        case 5 % Marker
            [markerread, markervals] = CEDS64ReadMarkers( fhand1, m, 100, 0 );
            CEDS64SetMarkerChan( fhand2, chan2, rate, 5 );
            CEDS64WriteMarkers( fhand2, chan2, markervals );
        case 6 % Wave Mark
            [ iOk, Rows, Cols ] = CEDS64GetExtMarkInfo( fhand1, m );
            [wmarkerread, wmarkervals] = CEDS64ReadExtMarks( fhand1, m, maxpoints, 0 );
            CEDS64SetExtMarkChan(fhand2, chan2, rate, 6, Rows, Cols, chandiv);
            CEDS64WriteExtMarks( fhand2, chan2, wmarkervals);
        case 7 % Real Mark
            [ iOk, Rows, Cols ] = CEDS64GetExtMarkInfo( fhand1, m );
            [rmarkerread, rmarkervals] = CEDS64ReadExtMarks( fhand1, m, maxpoints, 0 );
            CEDS64SetExtMarkChan( fhand2, chan2, rate, 7, Rows, Cols, -1);
            CEDS64WriteExtMarks( fhand2, chan2, rmarkervals);
        case 8 % Text Mark
            [ iOk, Rows, Cols ] = CEDS64GetExtMarkInfo( fhand1, m );
            [tmarkerread, tmarkervals] = CEDS64ReadExtMarks( fhand1, m, 100, 0 );
            CEDS64SetExtMarkChan( fhand2, chan2, rate, 8, Rows, Cols, -1 );
            CEDS64WriteExtMarks( fhand2, chan2, tmarkervals);
        case 9 % Realwave
            [floatread, floatvals, floattime] = CEDS64ReadWaveF( fhand1, m, maxpoints, 0 );
            CEDS64SetWaveChan( fhand2, chan2, chandiv, 9, rate );
            CEDS64WriteWave( fhand2, chan2, floatvals, floattime );
    end
    % copy units, comments, offsets etc.
    if (chan > 0)
        [ iOk, dOffset ] = CEDS64ChanOffset( fhand1, m );
        [ iOk ] = CEDS64ChanOffset( fhand2, chan2, dOffset );
        [ iOk, dScale ] = CEDS64ChanScale( fhand1, m );
        [ iOk ] = CEDS64ChanScale( fhand2, chan2, dScale );
        [ iOk, dTitle ] = CEDS64ChanTitle( fhand1, m );
        [ iOk ] = CEDS64ChanTitle( fhand2, chan2, dTitle );
        [ iOk, sUnits ] = CEDS64ChanUnits( fhand1, m );
        [ iOk ] = CEDS64ChanUnits( fhand2, chan2, sUnits );
    end
end

% write open ephys data to new channels in the SpikePlusOpenEphys.smr file
% get channel maps 
mapFile = sessionEphysInfo.mapFile;
load(fullfile('Z:\obstacleData\ephys\channelMaps\kilosort', [mapFile, '.mat']), 'channelNum_OpenEphys');
if ephysChanNumber ~= length(channelNum_OpenEphys)
    warning('channelNum_OpenEphys does NOT match ephysChanNumber!');
end


% function to extract voltage from binary file
getVoltage = @(data) ...
    double(data)*sessionEphysInfo.bitVolts; % extract voltage from memmapfile, converting to votlage, highpassing, and only return specific channel

% load data
contFiles = dir(fullfile(sessionEphysInfo.ephysFolder, '*.continuous'));
data = memmapfile(fullfile(sessionEphysInfo.ephysFolder, [contFiles(1).name(1:end-12), 's.dat']), ...
    'Format', {'int16', [sessionEphysInfo.channelNum sessionEphysInfo.smps], 'Data'}, 'Writable', false);

% calculate ephys start and stop time (in spike timestamps) and channel
% divide rate.
ephysStartTime = sessionEphysInfo.convertedEphysTimestamps(1);
ephysEndTime = sessionEphysInfo.convertedEphysTimestamps(end);
ephysSamps = length(sessionEphysInfo.convertedEphysTimestamps);
ephysChanDiv = (ephysEndTime - ephysStartTime)/(ephysSamps*timebase);
sTime = CEDS64SecsToTicks( fhand2, ephysStartTime );

ephysOriginalFs = 1/(ephysChanDiv*timebase);
ephysDesiredFs = 1/(round(ephysChanDiv)*timebase);
[p, q] = rat(ephysDesiredFs/ephysOriginalFs);

for i = 1:length(channelNum_OpenEphys)    
    
    disp(num2str(i));
    
    % create adc channel
    wavechan = CEDS64GetFreeChan( fhand2 );
    disp(['  wavechan = ', num2str(wavechan)]);
    % fs = 1/mean(diff(sessionEphysInfo.convertedEphysTimestamps));
    createret = CEDS64SetWaveChan( fhand2, wavechan, round(ephysChanDiv), 1, ephysDesiredFs);
    if createret ~= 0, warning('waveform channel not created correctly'); end
    CEDS64ChanTitle( fhand2, wavechan, ['ch ' num2str(i)]); % 1 being the most top recording site on the probe, 32/64 being the most bottom one.
    CEDS64ChanComment( fhand2, wavechan, '');
    CEDS64ChanUnits( fhand2, wavechan, 'Volt' );
    
    % prepare ephys data for writing into the file
    waveDataVoltage = getVoltage(data.Data.Data(channelNum_OpenEphys(i), :));
    waveDataResampled = resample(waveDataVoltage,p,q); % resample the data to fit into the timebase of the spike file
    waveDataNew = int16(waveDataResampled/sessionEphysInfo.bitVolts); % convert back to int16 type
    
    
    % write ephys data into the file
    fillret = CEDS64WriteWave( fhand2, wavechan, waveDataNew, sTime );
    if fillret < 0, warning(['Wave Channel ', num2str(wavechan), ' not filled correctly']); end
    
%      [iOk, ~, ~] = CEDS64ChanYRange( fhand2, wavechan, -0.5, 0.5);
%      if iOk < 0, warning('Chan Y range not set correctly'); end
%     
%     if strcmp(class(waveData), 'int16')
%         % write to the channels
%         ephysTimeOffset = sessionEphysInfo.convertedEphysTimestamps(1);
%         sTime = CEDS64SecsToTicks( fhand2, ephysTimeOffset ); %
%         % fillret = CEDS64WriteWave( fhand2, wavechan, waveData, sTime );
%         fillret = calllib('ceds64int', 'S64WriteWaveS', fhand2, wavechan, waveData, length(waveData), sTime);
%         if fillret < 0, warning('waveform channel not filled correctly'); end;
%         clear waveData;
%     else
%         warning('waveform channel not filled correctly, waveData is not int16 class!');
%     end
    
end

%
% close the file
CEDS64CloseAll();
% unload ceds64int.dll
unloadlibrary ceds64int;

end
