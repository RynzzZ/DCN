function resizeVideo(vidPath, oldVidName, newVidName, newWidth, newHeight, topOffset, leftOffset)

% load the video
addpath(fullfile(vidPath));
v = VideoReader(fullfile(vidPath, oldVidName));

% check if the new dimensions and offsets work
if topOffset + newHeight > v.Height
    error('topOffset + newHeight + 1 exceeds the original height!');
end

if leftOffset + newWidth > v.Width
    error('leftOffset + newWidth + 1 exceeds the original width!');
end

% convert new dimension and offsets into start and stop index to resize
% each frame
widStartInd = 1 + leftOffset;
widStopInd = leftOffset + newWidth;
hgtStartInd = 1 + topOffset;
hgtStopInd = topOffset + newHeight;

% open video writter
vidWriter = VideoWriter(fullfile(vidPath, newVidName));
vidWriter.FrameRate = v.FrameRate;
open(vidWriter);

% resize the vid
progressReport = 0.1;
for i = 1:v.NumFrames
    frame = read(v, i);  
    frameNew = frame(hgtStartInd:hgtStopInd, widStartInd:widStopInd, :);
    writeVideo(vidWriter, frameNew);
    
    if i/v.NumFrames > progressReport
        disp([num2str(100*progressReport) '%']);
        progressReport = progressReport + 0.1;
    end
end
disp('100%');


end