%YOLOv3
clc;
close all;
clear;
img_input = 'visionteam1.jpg'; % change here


% Specify the name of a pretrained YOLO v3 deep learning network .
name = 'tiny-yolov3-coco'; % big: 'darknet53-coco' <> small: 'tiny-yolov3-coco'

% Create YOLO v3 object detector by using the pretrained YOLO v3 network .
YOLO = yolov3ObjectDetector(name);

% Display and inspect the properties of the YOLO v3 object detector
%disp(YOLO)

% use analyzeNetwork to display the YOLO v3 network architecture
% and get information about the network layers . The network has
% two detection heads attached to the feature extraction network.
%analyzeNetwork(YOLO.Network)

I = imread(img_input);
% put image into YOLO format
Iyolo = preprocess(YOLO, I);
Iyolo = im2single(Iyolo);
% do detection
tic;
[bboxes, scores, labels] = detect(YOLO, Iyolo, 'DetectionPreprocessing', 'none');
t_yolo = toc;
fprintf("\nyolov3 on [ %s ]  ,  processing time:  %0.3f seconds\n\n", img_input, t_yolo);

detectedImg = insertObjectAnnotation(Iyolo, 'Rectangle', bboxes, labels);
figure;
imshow(detectedImg)