clc;
close all;
clear;
% SSD
vehicleDetector = load ('ssdVehicleDetector.mat', 'detector');
detector_ssd = vehicleDetector . detector ;
I = imread ('highway.png');

tic
[bboxes, scores, labels] = detect(detector_ssd ,I);
t_ssd = toc

if ~ isempty(bboxes)
    detectedI = insertObjectAnnotation(I, 'rectangle', bboxes ,cellstr(labels));
else
    detectedI = insertText (I ,[10 10] , 'No Detections');
end

figure
imshow(detectedI)
