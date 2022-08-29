function generateVideoClips(vidPath, oldVidName, newVidName, startTime, stopTime)


% load the video
addpath(fullfile(vidPath));
v = VideoReader(fullfile(vidPath, oldVidName));

% open video writter
vidWriter = VideoWriter(fullfile(vidPath, newVidName), 'MPEG-4');
vidWriter.FrameRate = v.FrameRate;
open(vidWriter);

% clip the video
startFrame = int32((floor(startTime)*60 + (startTime - floor(startTime))*100)*v.FrameRate);
stopFrame = int32((floor(stopTime)*60 + (stopTime - floor(stopTime))*100)*v.FrameRate);

disp('start clipping...')
for i = startFrame:stopFrame
    frame = read(v, i);  
    writeVideo(vidWriter, frame);
end
disp('done!');


end