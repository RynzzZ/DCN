function preprocessSessionData(session, varargin)

s.node = '109';
s.analyze = 'spike';
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs


convertOpenEphysToSpikeTimes(session, 'node', s.node);
formatSpikePlusOpenEphysFile(session);
analyzeSession(session, 'analyze', s.analyze);


end

