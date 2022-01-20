function addOneChannelToSpike(folderPath, oldFileName, newFileName, data, timestamps, varargin)

% settings
s.channelTitle = 'Jaw';
s.channelUnit = 'Pixel';

if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs

% initialization 
rootFolder = 'Z:\Qianyun\DCN\';
gitFolder = 'D:\DCN_Project\Github\DCN';
filePath = fullfile(folderPath, oldFileName);
newFilePath = fullfile(folderPath, newFileName);

% add path to CED code
setenv('CEDS64ML', 'D:\DCN_Project\Github\DCN\Spike2\Spike2_MATLAB_Interface\CEDS64ML');

cedpath = getenv('CEDS64ML');
addpath(cedpath);
% load ceds64int.dll
CEDS64LoadLib( cedpath );

% Open the old .smr spike file in the session folder and copy it to a new
% .smr file.
fhand1 = CEDS64Open(filePath);
if (fhand1 <= 0);  CEDS64ErrorMessage(fhand1); unloadlibrary ceds64int; return; end
maxchans = CEDS64MaxChan( fhand1 );

% Create the new .smr file - Spike + OpenEphys daya together
fhand2 = CEDS64Create(newFilePath, maxchans + 1, 2 );
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


% calculate channel divide rate for the new channel
dataStartTime = timestamps(find(~isnan(timestamps), 1, 'first'));
dataEndTime = timestamps(find(~isnan(timestamps), 1, 'last'));
dataSamps = length(timestamps);
dataChanDiv = (dataEndTime - dataStartTime)/(dataSamps*timebase);
sTime = CEDS64SecsToTicks( fhand2, dataStartTime );

dataOriginalFs = 1/(dataChanDiv*timebase);
dataDesiredFs = 1/(round(dataChanDiv)*timebase);
[p, q] = rat(dataDesiredFs/dataOriginalFs);

% add the data channel to the original spike .smr file
wavechan = CEDS64GetFreeChan( fhand2 );
disp(['  newchan = ', num2str(wavechan)]);
createret = CEDS64SetWaveChan( fhand2, wavechan, round(dataChanDiv), 1, dataDesiredFs);
if createret ~= 0, warning('waveform channel not created correctly'); end
CEDS64ChanTitle( fhand2, wavechan, s.channelTitle);
CEDS64ChanComment( fhand2, wavechan, '');
CEDS64ChanUnits( fhand2, wavechan, s.channelUnit );

% prepare data for writing into the file
waveDataResampled = resample(data,p,q); % resample the data to fit into the timebase of the spike file
waveDataNew = int16(waveDataResampled); % convert back to int16 type

% write the data into the file
fillret = CEDS64WriteWave( fhand2, wavechan, waveDataNew, sTime );
if fillret < 0, warning(['Wave Channel ', num2str(wavechan), ' not filled correctly']); end

% close the file
CEDS64CloseAll();
% unload ceds64int.dll
unloadlibrary ceds64int;

end