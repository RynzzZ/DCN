function [result] = removeOptoArtifacts(data, bottomRange, topRange)

% to get rid of the artifacts caused by giving opto stimuli.
% data: raw ephys data;
% bottomRange: range below zero where opto artifacts exist, eg, [-400,
% -100];
% topRange: range above zero where opto artifacts exist, eg, [80, 200].

if any(bottomRange)
    data(data >= bottomRange(1) & data <= bottomRange(2)) = nan;
    
    inds = find(~isnan(data));
    result = interp1(inds, data(inds), 1:length(data), 'nearest');
end

if any(topRange)
    result(result >= topRange(1) & result <= topRange(2)) = nan;
    inds = find(~isnan(result));
    result = interp1(inds, result(inds), 1:length(result), 'nearest');
end


end

