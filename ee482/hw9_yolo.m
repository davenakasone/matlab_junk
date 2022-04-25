clc;
close all;
clear;
%YOLOv3
% Specify the name of a pretrained YOLO v3 deep learning network .
name = 'tiny-yolov3-coco'; %or bigger model 'darknet53 -coco '

% Create YOLO v3 object detector by using the pretrained YOLO v3 network .
YOLO = yolov3ObjectDetector(name);

% Display and inspect the properties of the YOLO v3 object detector
disp(YOLO)

% use analyzeNetwork to display the YOLO v3 network architecture
% and get information about the network layers . The network has
% two detection heads attached to the feature extraction network.
analyzeNetwork(YOLO.Network)

I = imread('barcelona-team.jpg');
% put image into YOLO format
Iyolo = preprocess(YOLO, I);
% do detection
tic
[bboxes, scores, labels] = detect(YOLO, Iyolo , 'DetectionPreprocessing', 'none');
t_yolo = toc

detectedImg = insertObjectAnnotation(Iyolo, 'Rectangle', bboxes, labels);
figure
imshow(detectedImg)