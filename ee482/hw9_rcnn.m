clc;
close all;
clear;
% fast RCNN
data = load('fasterRCNNVehicleTrainingData.mat', 'detector');
detector_fastrcnn = data.detector ;

% Display and inspect the properties of object detector .
disp(detector_fastrcnn)

% use analyzeNetwork to display the network architecture and get
% information about the network layers . The network has two
% detection heads attached to the feature extraction network .
analyzeNetwork(detector_fastrcnn.Network)

I = imread('highway.png');
tic
[bboxes, scores, labels] = detect(detector_fastrcnn ,I);
t_fastrcnn = toc
if ~ isempty (bboxes)
    detectedI = insertObjectAnnotation (I,'rectangle',bboxes, cellstr(labels));
else
    detectedI = insertText (I ,[10 10] ,'No Detections');
end

figure
imshow(detectedI)