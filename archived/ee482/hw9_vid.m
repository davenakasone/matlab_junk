clc;
close all;
clear;
%{
zz = 1;
vehicleDetector = load('ssdVehicleDetector.mat', 'detector');
detector_ = vehicleDetector.detector;
%}
%{
zz = 2;
data = load('fasterRCNNVehicleTrainingData.mat', 'detector');
detector_ = data.detector;
%}
%
zz = 3;
name = 'tiny-yolov3-coco';
detector_ = yolov3ObjectDetector(name);
%}

hits = 0;
% load from url
%url = 'http://www.ee.unlv.edu/~b1morris/labs/cv_intro/dash02.mp4';
%outfilename = websave('dash02.mp4', url);
videoSource = VideoReader('dash02.mp4');
vidout = uint8(zeros(240, 360, 3, 9));
t = zeros(1 ,9); % processing time vector

% try out different set of frames to see performance.
index_start = 1;
index_stop = 9;
index = [index_start, index_stop];
frames = read(videoSource, index);
%
for iFrame = index_start:1:index_stop
  img_out = read(videoSource, iFrame);
  imwrite(img_out, sprintf('dash_%d.jpg', iFrame));
end 
%}
for i = index_start:1:index_stop
    frame = frames(:, :, :, i);

    % crop image to center
    frame_crop = frame(240-120:240+120-1 , 427-180:427+180-1 , :);

    tic
    [bboxes, scores, labels] = detect(detector_, frame_crop);
    t(i) = toc;
  
    if ~ isempty(bboxes)
        detectedI = insertObjectAnnotation(frame_crop, 'rectangle', bboxes, cellstr(labels));
        temp = size(bboxes);
        hits = hits + temp(1);
    else
        detectedI = insertText (frame_crop, [10 10], 'No Detections');
    end
    
    vidout(:, :, :, i) = detectedI;
end

% display vid output as frames
for i = index_start:1:index_stop
    figure
    imshow(vidout(:, :, :, i) ,[]);
    title(sprintf('frame = %d', index(1)+i-1));
end

tmean = mean(t);
fps = 1 / tmean;
fprintf (1, 'Average FPS = %2.2f (%2.4f s)\n\tdetections:  %d\n', fps ,tmean, hits);
if (zz==1)
    fprintf("using SSD\n\n");
end
if (zz==2)
    fprintf("using RCNN\n\n");
end
if (zz==3)
    fprintf("using YOLOv3\n\n");
end


%%%%~~~~END>