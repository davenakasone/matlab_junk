clc;
close all;
clear;
%dash video
vehicleDetector = load ('ssdVehicleDetector.mat', 'detector');
detector_ssd = vehicleDetector.detector;

% load from url
url = 'http://www.ee.unlv.edu/~b1morris/labs/cv_intro/dash02.mp4';
outfilename = websave('dash02.mp4', url);

videoSource = VideoReader('dash02.mp4');

vidout = uint8(zeros(240, 360, 3, 9));
t = zeros(1 ,9); % processing time vector

% try out different set of frames to see performance .
index = [1, 9];
frames = read(videoSource, index);

for i = 1:9
    frame = frames(:, :, :, i);

    % crop image to center
    frame_crop = frame(240-120:240+120-1 , 427-180:427+180-1 , :);

    %[ TODO be sure to try both the SSD and Fast RCNN models ]
    tic
    [bboxes, scores, labels] = detect(detector_ssd, frame_crop);
    t(i) = toc;

    if ~ isempty(bboxes)
        detectedI = insertObjectAnnotation(frame_crop, 'rectangle', bboxes, cellstr(labels));
    else
        detectedI = insertText (frame_crop, [10 10], 'No Detections');
    end
    
    vidout(:, :, :, i) = detectedI;
end

tmean = mean(t);
fps = 1 / tmean;
fprintf (1, 'Average FPS = %2.2f (%2.4f s)\n', fps ,tmean);

% display vid output as frames
for i = 1:9
    figure
    imshow (vidout(:, :, :, i) ,[]);
    title(sprintf('frame = %d', index(1)+i-1));
end