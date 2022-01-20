function [cros,locs,ncrs] = crossdet(sig, th, type)
%CROSSDET  Detect crossing a threshold into a signal
%   returns cross values and its location that occurred in a input signal.
%   It can detect next first points crossing up through or down through the
%   threshold. Signal data requires a row or column vector with real-valued 
%   elements. If there are no local crossing returns empty.
%
% Syntax:
%   [cros,locs,ncrs] = crossdet(sig, th, type)
%
% Inputs:
%   sig - Signal column vector
%    th - threshold scalar absolute value
%  type - Detection crossing 'thresholdUp' or 'thresholdDown'
%
% Outputs:
%  cros - signal values after each crossing
%  locs - signal indices after each crossing
%  ncrs - number of crossings
%
% Example:
% t = [0:.01:1]';                           % time vector
% s = -sin(2*pi*5*t);                       % signal vector
% s = s + .6*randn(length(s),1);            % add noise
% th = 0.5;                                 % choose th
% [crosUP,locsUP] = crossdet(s, th);        % perform crossing-up detection
% [crosDOWN,locsDOWN] = crossdet(-s, th);   % perform crossing-down detection
% plot(t,s,'k'); hold on;                   % plot signal
% stem(t(locsUP),crosUP,'g')                % plot up crossings
% stem(t(locsDOWN),-crosDOWN,'r')           % plot down crossings
% line([t(1) t(end)],[th th], ...           % plot thresholds
%     'LineStyle','--','color','m');
% line([t(1) t(end)],[-th -th],'LineStyle','--','color','m');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: peakdet;

% Author: Marco Borges, Ph.D. Student, Computer/Biomedical Engineer
% UFMG, PPGEE, Neurodinamica Lab, Brazil
% email address: marcoafborges@gmail.com
% Website: http://www.cpdee.ufmg.br/
% Jul 2014; Version: v2; Last revision: 2014-07-19
% References:
%
% Changelog:
%   v2 - bugfix handle only one cross and return number of crossings
%
% Edit:
%   v3 - split 'threshold' into 'thresholdUp' and 'thresholdDown'; make this
%        function able to work with either type of threshold.
%------------------------------- BEGIN CODE -------------------------------
cros = [];
locs = [];

[m,n] = size(sig);
if m == 1 && n > 1
    sig = sig';
    [m,n] = size(sig);
end

if nargin < 2
    th = 1;
end

if length(th) > 1
    error('Input argument th (threshold) must be a scalar');
end

if ~ischar(type)
    error('Input argument type must be string "upthreshold" or "downthreshold".');
end

if strcmp(type,'thresholdUp') || strcmp(type, 'thresholdDown')
    cross = th;
else
    error('Input argument type must be string "upthreshold" or "downthreshold".');
end

switch type
    case 'thresholdUp'        
        idx = find(sig > th,1,'first');
        if ~isempty(idx) && idx < m
            idd = find(sig(idx+1:end) < cross,1,'first');
            idd = idd + idx;
            ncrs = 1;
        else
            ncrs = 0;
            idd = [];
        end
        
        while ~isempty(idx) && ~isempty(idd) || ncrs == 1
            cros = [cros; sig(idx)];
            locs = [locs; idx];
            idx = find(sig(idd+1:end) > th,1,'first');
            idx = idx + idd;
            if idx < m
                idd = find(sig(idx+1:end) < cross,1,'first');
                idd = idd + idx;
                ncrs = ncrs +1;
            else
                idd = [];
            end
        end
        
    case 'thresholdDown'
        idx = find(sig < th,1,'first');
        if ~isempty(idx) && idx < m
            idd = find(sig(idx+1:end) > cross,1,'first');
            idd = idd + idx;
            ncrs = 1;
            round = 1;
        else
            ncrs = 0;
            idd = [];
            idx = [];
            round = 0;
        end
        
        while ~isempty(idx) && ~isempty(idd) || round == 1
            cros = [cros; sig(idx)];
            locs = [locs; idx];
            idx = find(sig(idd+1:end) < th,1,'first');
            idx = idx + idd;
            if idx < m
                idd = find(sig(idx+1:end) > cross,1,'first');
                idd = idd + idx;
                ncrs = ncrs +1;
            else
                idd = [];         
            end
            round = round + 1;
        end
end

end
%-------------------------------- END CODE --------------------------------