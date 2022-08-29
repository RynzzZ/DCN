function [optoOnOffTimes] = getOptoSpFROnOffTimes(spike, varargin)
% output: optoOnOffTimes - column1: optoOnTimes;
%                          column2: optoOffTimes;   

% only take into consideration opto trains given to affect spontaneous
% firing rate and get its onset and offset times. Remove any light on periods 
% happening during auditory stimuli or feeding behavior (s.expType =
% 'awake').

s.expType = 'anes'; % set to 'awake' to also remove any opto given during chewing.

if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs

optoTriggerTimes = spike.optoTrainTimes;
audioTriggerTimes = spike.audioTriggerTimes;

backup = optoTriggerTimes;

for i = 1:length(optoTriggerTimes)
    
    if any(audioTriggerTimes >= optoTriggerTimes(i) - 5 & audioTriggerTimes <= optoTriggerTimes(i) + 5)
        optoTriggerTimes(i) = 0;
    end
end

optoTriggerTimes(optoTriggerTimes == 0) = [];

inds = find(diff(optoTriggerTimes)>1) + 1;    
optoOnTimes = [optoTriggerTimes(1); optoTriggerTimes(inds)];
optoOffTimes = [optoTriggerTimes(inds-1); optoTriggerTimes(end)];

if strcmp(s.expType, 'awake')
    inds = strfind(convertCharsToStrings(spike.keyboardInput), 'F');
    Ftimes = spike.keyboardTimes(inds);
    
    for i = 1:length(optoOnTimes)
        if any(abs(Ftimes - optoOnTimes(i)) < 15)
            optoOnTimes(i) = 0;
            optoOffTimes(i) = 0;
        end
    end
    optoOnTimes(optoOnTimes == 0) = [];
    optoOffTimes(optoOffTimes == 0) = [];
end

optoOnOffTimes = [optoOnTimes, optoOffTimes];
       
% quality check
figure('Color', 'white', 'position', get(0,'ScreenSize')); clf;
plot(backup, ones(size(backup)), '.'); 
hold on; 

plot(optoOnTimes, ones(size(optoOnTimes))*0.99, '.g'); 
plot(optoOffTimes, ones(size(optoOffTimes))*0.99, '.r'); 
ylim([0.8, 1.1])
end

