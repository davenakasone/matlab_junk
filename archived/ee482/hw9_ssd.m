% SSD
clc;
close all;
clear;
img_input = 'visionteam1.jpg'; % change here
%img_input = 'highway.jpg';

vehicleDetector = load ('ssdVehicleDetector.mat', 'detector');
detector_ssd = vehicleDetector.detector;

I = imread (img_input);
tic
[bboxes, scores, labels] = detect(detector_ssd ,I);
p_time = toc;

if ~ isempty(bboxes)
    detectedI = insertObjectAnnotation(I, 'rectangle', bboxes, cellstr(labels));
else
    detectedI = insertText (I ,[10 10] , 'No Detections');
end

figure;
imshow(detectedI);

temp = size(bboxes);
fprintf("\nSSD on [ %s ]  ,  detections:  %d\n", img_input, temp(1));
fprintf("processing time:  %0.3f seconds\n\n", p_time);


%%%%~~~~END>
