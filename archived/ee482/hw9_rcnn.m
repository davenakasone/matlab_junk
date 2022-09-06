% fast RCNN
clc;
close all;
clear;
img_input = 'visionteam1.jpg'; % change here
%img_input = 'highway.jpg';

data = load('fasterRCNNVehicleTrainingData.mat', 'detector');
detector_fastrcnn = data.detector;

% Display and inspect the properties of object detector .
%disp(detector_fastrcnn)

% use analyzeNetwork to display the network architecture and get
% information about the network layers . The network has two
% detection heads attached to the feature extraction network .
%analyzeNetwork(detector_fastrcnn.Network)

I = imread(img_input);
tic;
[bboxes, scores, labels] = detect(detector_fastrcnn ,I);
p_time = toc;
if ~ isempty (bboxes)
    detectedI = insertObjectAnnotation (I,'rectangle',bboxes, cellstr(labels));
else
    detectedI = insertText (I ,[10 10] ,'No Detections');
end
figure;
imshow(detectedI);

temp = size(bboxes);
fprintf("\nRCNN on [ %s ]  ,  detections:  %d\n", img_input, temp(1));
fprintf("processing time:  %0.3f seconds\n\n", p_time);


%%%%~~~~END>