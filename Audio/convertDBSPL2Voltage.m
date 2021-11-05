function Voltage = convertDBSPL2Voltage(DBSPL, varargin)

% settings
s.Vref = 2.8117*2;  % Vpp of the reference
s.DBSPLref = 70;  % dB SPL level of the reference 

if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs


Vrefrms = s.Vref/sqrt(2);
Vrms = 10^((DBSPL - s.DBSPLref)/20)*Vrefrms;
Vpp = Vrms*sqrt(2);
Voltage = Vpp/2;


end