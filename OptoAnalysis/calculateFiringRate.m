function [FR, FRTimes] = calculateFiringRate(ephysSpkTimes, ephysTimes, fs, varargin)


s.ephysFs = 30000;

if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs


nbins = round(ephysTimes(end) - ephysTimes(1))*fs;
binWidth = 1/fs;
[counts, centers] = hist(ephysSpkTimes, nbins);

h = 1;
r=ksr(centers, counts/binWidth, h);
FR = r.f;
FRTimes = r.x;

figure('Color', 'white','WindowState','maximized'); clf;
plot(centers, counts/binWidth, 'co', centers, counts/binWidth, 'b-', r.x,r.f,'r--','linewidth',2);
legend('FR Hist', 'FR Line', 'FR Regression','location','northwest');
title('FR Smoothed By Gaussian kernel regression');


end

