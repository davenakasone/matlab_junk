% SSD
clc;
close all;
clear all;
img_input = 'visionteam1.jpg'; % change here

%vehicleDetector = load ('ssdVehicleDetector.mat', 'detector');
%detector_ssd = vehicleDetector.detector;

vehicleDetector = load ('ssdObjectDetector.mat', 'detector');
detector_ssd = objectDetector.detector;



I = imread (img_input);
tic
[bboxes, scores, labels] = detect(detector_ssd ,I);
t_ssd = toc;

if ~ isempty(bboxes)
    detectedI = insertObjectAnnotation(I, 'rectangle', bboxes, cellstr(labels));
else
    detectedI = insertText (I ,[10 10] , 'No Detections');
end
fprintf("\nSSD on [ %s ]  ,  processing time:  %0.3f seconds\n\n", img_input, t_ssd);
figure
imshow(detectedI)
