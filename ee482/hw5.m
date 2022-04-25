%{
    1  : problem 1a, comapring merge threshold 4 to 0, 'FrontalFaceCART'(Default)
    2  : problem 1b, "                              ", 'FrontalFaceLBP'
    3  : problem 1c, using problem 1b on "barcelona-team.jpeg"
    4  : problem 1c, using problem 1b on "SJEarthquakesteampic.jpeg"
    5  : problem 2a, HOG, vary parameters, get "best" with default model
    6  : problem 2b, HOG, use best parameters, use "best", change WindowStride
    7  : problem 2c, HOG, "best", changing window stride, "barcelona-team.jpeg"
    8  : problem 2c, HOG, "best", changing window stride, "SJEarthquakesteampic.jpeg"
    9  : problem 3b, train model, stop signs, 5-stage, mt=1
    10 : problem 3c, Haar
    11 : problem 3c, HOG
    12 : problem 3c, LBP
    13 : problem 3d, YOLOv2
    14 : problem 3d, ACF
    15 : problem 3d, fast RCNN
    16 : problem 3e, train custom

    97 : original problem 1 code
        more parameters https ://www.mathworks.com/help/vision/ref/vision.cascadeobjectdetector-system-object.html
    98 : original problem 2 code
        more parameters https://www.mathworks.com/help/vision/ref/vision .peopledetector-system-object.html
%}
format compact;
clc;
close all;
clear all;
select = 66;


%------------------------------------------------------------------------------------------
if select == 1
    I = imread ('visionteam.jpg');
    MergeThreshold = 0;
    VJ_0 = vision.CascadeObjectDetector('MergeThreshold', MergeThreshold);
    tic
    bboxes = VJ_0(I);
    t_VJ_0 = toc
    Ifaces = insertObjectAnnotation (I, 'rectangle', bboxes , 'Face');
    figure(1);
    hold on;
    imshow (Ifaces , []);
    title ('Viola Jones Face Detector, threshold=0');
    hold off;
    release(VJ_0);
    %%%%~~~~
    I = imread ('visionteam.jpg');
    MergeThreshold = 4;
    VJ_4 = vision.CascadeObjectDetector('MergeThreshold', MergeThreshold);
    tic
    bboxes = VJ_4(I);
    t_VJ_4 = toc
    Ifaces = insertObjectAnnotation (I, 'rectangle', bboxes , 'Face');
    figure(2);
    hold on;
    imshow (Ifaces , []);
    title ('Viola Jones Face Detector, threshold=4');
    hold off;
    release(VJ_4);
end


%------------------------------------------------------------------------------------------
if select == 2
    I = imread ('visionteam.jpg');
    MergeThreshold = 0;
    VJ_0 = vision.CascadeObjectDetector('FrontalFaceLBP', 'MergeThreshold', MergeThreshold);
    tic
    bboxes = VJ_0(I);
    t_VJ_0 = toc
    Ifaces = insertObjectAnnotation (I, 'rectangle', bboxes , 'Face');
    figure(1);
    hold on;
    imshow (Ifaces , []);
    title ('Viola Jones Face Detector, threshold=0');
    hold off;
    release(VJ_0);
    %%%%~~~~
    I = imread ('visionteam.jpg');
    MergeThreshold = 4;
    VJ_4 = vision.CascadeObjectDetector('FrontalFaceLBP', 'MergeThreshold', MergeThreshold);
    tic
    bboxes = VJ_4(I);
    t_VJ_4 = toc
    Ifaces = insertObjectAnnotation (I, 'rectangle', bboxes , 'Face');
    figure(2);
    hold on;
    imshow (Ifaces , []);
    title ('Viola Jones Face Detector, threshold=4');
    hold off;
    release(VJ_4);
end


%------------------------------------------------------------------------------------------
if select == 3
    I = imread ('barcelona-team.jpeg');
    MergeThreshold = 0;
    VJ_0 = vision.CascadeObjectDetector('FrontalFaceLBP', 'MergeThreshold', MergeThreshold);
    tic;
    bboxes = VJ_0(I);
    t_VJ_0 = toc;
    fprintf("\nthe threshold=0 processing time:  %f seconds\n", t_VJ_0);
    Ifaces = insertObjectAnnotation (I, 'rectangle', bboxes , 'Face');
    figure(1);
    hold on;
    imshow (Ifaces , []);
    title ('Viola Jones Face Detector, threshold=0');
    hold off;
    release(VJ_0);
    %%%%~~~~
    I = imread ('barcelona-team.jpeg');
    MergeThreshold = 4;
    VJ_4 = vision.CascadeObjectDetector('FrontalFaceLBP', 'MergeThreshold', MergeThreshold);
    tic;
    bboxes = VJ_4(I);
    t_VJ_4 = toc;
    fprintf("\nthe threshold=4 processing time:  %f seconds\n", t_VJ_4);
    Ifaces = insertObjectAnnotation (I, 'rectangle', bboxes , 'Face');
    figure(2);
    hold on;
    imshow (Ifaces , []);
    title ('Viola Jones Face Detector, threshold=4');
    hold off;
    release(VJ_4);
end


%------------------------------------------------------------------------------------------
if select == 4
    I = imread ('SJEarthquakesteampic.jpeg');
    MergeThreshold = 0;
    VJ_0 = vision.CascadeObjectDetector('FrontalFaceLBP', 'MergeThreshold', MergeThreshold);
    tic;
    bboxes = VJ_0(I);
    t_VJ_0 = toc;
    fprintf("\nthe threshold=0 processing time:  %f seconds\n", t_VJ_0);
    Ifaces = insertObjectAnnotation (I, 'rectangle', bboxes , 'Face');
    figure(1);
    hold on;
    imshow (Ifaces , []);
    title ('Viola Jones Face Detector, threshold=0');
    hold off;
    release(VJ_0);
    %%%%~~~~
    I = imread ('SJEarthquakesteampic.jpeg');
    MergeThreshold = 4;
    VJ_4 = vision.CascadeObjectDetector('FrontalFaceLBP', 'MergeThreshold', MergeThreshold);
    tic;
    bboxes = VJ_4(I);
    t_VJ_4 = toc;
    fprintf("\nthe threshold=4 processing time:  %f seconds\n", t_VJ_4);
    Ifaces = insertObjectAnnotation (I, 'rectangle', bboxes , 'Face');
    figure(2);
    hold on;
    imshow (Ifaces , []);
    title ('Viola Jones Face Detector, threshold=4');
    hold off;
    release(VJ_4);
end


%------------------------------------------------------------------------------------------
if select == 5
    I_bad = imread('visionteam.jpg');
    ClassificationThreshold = 0;
    MergeDetections = false ;
    HOG_bad = vision.PeopleDetector( 'ClassificationThreshold',...
        ClassificationThreshold , 'MergeDetections', MergeDetections);
    tic;
    [bboxes_bad , scores_bad] = HOG_bad(I_bad);
    t_HOG_bad = toc;
    fprintf("\nprocessing time, bad parameters:  %f\n", t_HOG_bad);
    Iped_bad = insertObjectAnnotation (I_bad,'rectangle', bboxes_bad , scores_bad);
    figure(1);
    hold on;
    imshow (Iped_bad , []);
    title ('HOG Pedestrian Detector, bad parameters');
    hold off;
    release (HOG_bad);
    %%%%~~~~
    I_good = imread('visionteam.jpg');
    ClassificationThreshold = 0.7;
    MergeDetections = true;
    HOG_good = vision.PeopleDetector( 'ClassificationThreshold',...
        ClassificationThreshold , 'MergeDetections', MergeDetections);
    tic;
    [bboxes_good , scores_good] = HOG_good(I_good);
    [selected_bbox, selected_score] = selectStrongestBbox(bboxes_good, scores_good);
    t_HOG_good = toc;
    fprintf("\nprocessing time, good parameters:  %f\n", t_HOG_good);
    Iped_good = insertObjectAnnotation (I_good,'rectangle', selected_bbox , selected_score);
    figure(2);
    hold on;
    imshow (Iped_good , []);
    title ('HOG Pedestrian Detector, good parameters');
    hold off;
    release (HOG_good);
end


%------------------------------------------------------------------------------------------
if select == 6
    for ii = 2:2:10
        I = imread('visionteam.jpg');
        ClassificationThreshold = 0;
        MergeDetections = true;
        HOG = vision.PeopleDetector(... 
            'ClassificationThreshold', ClassificationThreshold,...
            'MergeDetections', MergeDetections,...
            'WindowStride', [ii,ii]);
        tic;
        [bboxes , scores] = HOG(I);
        t_HOG = toc;
        fprintf("\nprocessing time, WindowStride[%d,%d]:  %f\n", ii, ii, t_HOG);
        Iped = insertObjectAnnotation (I,'rectangle', bboxes , scores);
        figure(ii);
        hold on;
        imshow (Iped , []);
        temp = sprintf("HOG Pedestrian Detector, WindowStride[%d, %d]", ii, ii);
        title (temp, 'fontsize', 18);
        hold off;
        release (HOG);
    end
end


%------------------------------------------------------------------------------------------
if select == 7
    I = imread('barcelona-team.jpg');
    ClassificationThreshold = 0.2;
    MergeDetections = false;
    for ii = 3:1:9
        HOG = vision.PeopleDetector(... 
                'ClassificationThreshold', ClassificationThreshold,...
                'MergeDetections', MergeDetections,...
                'WindowStride', ii);
        tic;
        [bboxes , scores] = HOG(I);
        t_HOG = toc;
        Iped = insertObjectAnnotation (I,'rectangle', bboxes , scores);
        fprintf("\nprocessing time, WindowStride[%d,%d]:  %f\n", ii, ii, t_HOG);
        temp = sprintf("WindowStride[%d, %d]", ii, ii);
        figure('name', temp);
        hold on;
        imshow (Iped , []);
        title (temp, 'fontsize', 18);
        hold off;
        release (HOG);
        pause(1);
    end
end


%------------------------------------------------------------------------------------------
if select == 8
    I = imread('SJEarthquakesteampic.jpeg');
    ClassificationThreshold = 0.7;
    MergeDetections = false;
    for ii = 3:1:9
        HOG = vision.PeopleDetector(... 
                'ClassificationThreshold', ClassificationThreshold,...
                'MergeDetections', MergeDetections,...
                'WindowStride', ii);
        tic;
        [bboxes , scores] = HOG(I);
        t_HOG = toc;
        Iped = insertObjectAnnotation (I,'rectangle', bboxes , scores);
        fprintf("\nprocessing time, WindowStride[%d,%d]:  %f\n", ii, ii, t_HOG);
        temp = sprintf("WindowStride[%d, %d]", ii, ii);
        figure('name', temp);
        hold on;
        imshow (Iped , []);
        title (temp, 'fontsize', 18);
        hold off;
        release (HOG);
        pause(1);
    end
end


%------------------------------------------------------------------------------------------
if select == 97
%Viola Jones Face Detector
    % read image
    I = imread ('visionteam.jpg');
    MergeThreshold = 4;
    
    % create Viola Jones detector
    VJ = vision.CascadeObjectDetector ('MergeThreshold', MergeThreshold);

    % detect faces
    tic
    bboxes = VJ(I);
    % time detector
    t_VJ = toc

    % annotate detected faces
    Ifaces = insertObjectAnnotation (I, 'rectangle', bboxes , 'Face');
    figure;
    imshow (Ifaces , []);
    title ('Viola Jones Face Detector');
    
    % clean up
    release(VJ);
end


%------------------------------------------------------------------------------------------
if select == 98
%HOG detector
    % read image
    I = imread('visionteam.jpg');
    
    % increase value when there are too many false positives
    ClassificationThreshold = 0;
    
    % adjust to group multiple responses
    MergeDetections = false ;

    HOG = vision.PeopleDetector( 'ClassificationThreshold',...
        ClassificationThreshold , 'MergeDetections', MergeDetections);

    tic
    [bboxes , scores] = HOG(I);
    t_HOG = toc

    Iped = insertObjectAnnotation (I,'rectangle', bboxes , scores);

    figure ;
    imshow (Iped , []);
    title ('HOG Pedestrian Detector');
    % notice detections contain area around the person (bounding boxes are not tight)
    
    % cleanup
    release (HOG);
end


%------------------------------------------------------------------------------------------
if select == 9
    load('stopSignsAndCars.mat');
    stopSigns = fullfile(toolboxdir('vision'),'visiondata',stopSignsAndCars{:,1});
    imds = imageDatastore(stopSigns);
    blds = boxLabelDatastore(stopSignsAndCars(:,2));
    positiveInstances = combine(imds,blds);
    imDir = fullfile(matlabroot,'toolbox','vision','visiondata','stopSignImages');
    addpath(imDir);
    negativeFolder = fullfile(matlabroot,'toolbox','vision','visiondata','nonStopSigns');
    negativeImages = imageDatastore(negativeFolder);
    trainCascadeObjectDetector('stopSignDetector.xml', positiveInstances,...
        negativeFolder, FalseAlarmRate=0.1, NumCascadeStages=5);
    detector = vision.CascadeObjectDetector('stopSignDetector.xml', MergeThreshold=1);
    img = imread('stopSignTest.jpg');
    bbox = step(detector,img);
    detectedImg = insertObjectAnnotation(img,'rectangle',bbox,'stop sign');
    figure;
    imshow(detectedImg);
    rmpath(imDir);
end


%------------------------------------------------------------------------------------------
if select == 10
    load('stopSignsAndCars.mat');
    stopSigns = fullfile(toolboxdir('vision'),'visiondata',stopSignsAndCars{:,1});
    imds = imageDatastore(stopSigns);
    blds = boxLabelDatastore(stopSignsAndCars(:,2));
    positiveInstances = combine(imds,blds);
    imDir = fullfile(matlabroot,'toolbox','vision','visiondata','stopSignImages');
    addpath(imDir);
    negativeFolder = fullfile(matlabroot,'toolbox','vision','visiondata','nonStopSigns');
    negativeImages = imageDatastore(negativeFolder);
    trainCascadeObjectDetector('stopSignDetector.xml', positiveInstances, negativeFolder,...
        FalseAlarmRate=0.1, NumCascadeStages=5, FeatureType='Haar');
    detector = vision.CascadeObjectDetector('stopSignDetector.xml', MergeThreshold=4);

    ptime = zeros(1, 4);
    for ii = 1:1:4
        if ii == 1
            img = imread('s1_dashcam-screen-shot1.jpeg');
        end
        if ii == 2
            img = imread('s2_pedstop.jpeg');
        end
        if ii == 3
            img = imread('s3_jogger.jpeg');
        end
        if ii == 4
            img = imread('s4_four-way-stop-sign-one-way-american-suburb.jpeg');
        end
        tic;
        bbox = step(detector,img);
        detectedImg = insertObjectAnnotation(img,'rectangle',bbox,'stop sign');
        ptime(1, ii) = toc;
        figure(ii);
        hold on;
        title('Haar', 'fontsize', 18);
        imshow(detectedImg);
        hold off;
        rmpath(imDir);
    end
    for ii = 1:1:4
        fprintf("processing time, img_%d  :  %0.3f\n", ii, ptime(1,ii));
    end
end


%------------------------------------------------------------------------------------------
if select == 11
    load('stopSignsAndCars.mat');
    stopSigns = fullfile(toolboxdir('vision'),'visiondata',stopSignsAndCars{:,1});
    imds = imageDatastore(stopSigns);
    blds = boxLabelDatastore(stopSignsAndCars(:,2));
    positiveInstances = combine(imds,blds);
    imDir = fullfile(matlabroot,'toolbox','vision','visiondata','stopSignImages');
    addpath(imDir);
    negativeFolder = fullfile(matlabroot,'toolbox','vision','visiondata','nonStopSigns');
    negativeImages = imageDatastore(negativeFolder);
    trainCascadeObjectDetector('stopSignDetector.xml', positiveInstances, negativeFolder,...
        FalseAlarmRate=0.1, NumCascadeStages=30, FeatureType='HOG');
    detector = vision.CascadeObjectDetector('stopSignDetector.xml', MergeThreshold=8);

    ptime = zeros(1, 4);
    for ii = 1:1:4
        if ii == 1
            img = imread('s1_dashcam-screen-shot1.jpeg');
        end
        if ii == 2
            img = imread('s2_pedstop.jpeg');
        end
        if ii == 3
            img = imread('s3_jogger.jpeg');
        end
        if ii == 4
            img = imread('s4_four-way-stop-sign-one-way-american-suburb.jpeg');
        end
        tic;
        bbox = step(detector,img);
        detectedImg = insertObjectAnnotation(img,'rectangle',bbox,'stop sign');
        ptime(1, ii) = toc;
        figure(ii);
        hold on;
        title('HOG', 'fontsize', 18);
        imshow(detectedImg);
        hold off;
        rmpath(imDir);
    end
    for ii = 1:1:4
        fprintf("processing time, img_%d  :  %0.3f\n", ii, ptime(1,ii));
    end
end


%------------------------------------------------------------------------------------------
if select == 12
    load('stopSignsAndCars.mat');
    stopSigns = fullfile(toolboxdir('vision'),'visiondata',stopSignsAndCars{:,1});
    imds = imageDatastore(stopSigns);
    blds = boxLabelDatastore(stopSignsAndCars(:,2));
    positiveInstances = combine(imds,blds);
    imDir = fullfile(matlabroot,'toolbox','vision','visiondata','stopSignImages');
    addpath(imDir);
    negativeFolder = fullfile(matlabroot,'toolbox','vision','visiondata','nonStopSigns');
    negativeImages = imageDatastore(negativeFolder);
    trainCascadeObjectDetector('stopSignDetector.xml', positiveInstances, negativeFolder,...
        FalseAlarmRate=0.1, NumCascadeStages=22, FeatureType='LBP');
    detector = vision.CascadeObjectDetector('stopSignDetector.xml', MergeThreshold=4);

    ptime = zeros(1, 4);
    for ii = 1:1:4
        if ii == 1
            img = imread('s1_dashcam-screen-shot1.jpeg');
        end
        if ii == 2
            img = imread('s2_pedstop.jpeg');
        end
        if ii == 3
            img = imread('s3_jogger.jpeg');
        end
        if ii == 4
            img = imread('s4_four-way-stop-sign-one-way-american-suburb.jpeg');
        end
        tic;
        bbox = step(detector,img);
        detectedImg = insertObjectAnnotation(img,'rectangle',bbox,'stop sign');
        ptime(1, ii) = toc;
        figure(ii);
        hold on;
        title('LBP', 'fontsize', 18);
        imshow(detectedImg);
        hold off;
        rmpath(imDir);
    end
    for ii = 1:1:4
        fprintf("processing time, img_%d  :  %0.3f\n", ii, ptime(1,ii));
    end
end


%------------------------------------------------------------------------------------------
if select == 13
    detector = vehicleDetectorYOLOv2();
    ptime = zeros(1, 4);
    for ii = 1:1:4
        if ii == 1
            img = imread('s1_dashcam-screen-shot1.jpeg');
        end
        if ii == 2
            img = imread('s2_pedstop.jpeg');
        end
        if ii == 3
            img = imread('s3_jogger.jpeg');
        end
        if ii == 4
            img = imread('s4_four-way-stop-sign-one-way-american-suburb.jpeg');
        end
        tic;
        [bboxes, scores] = detect(detector, img);
        detectedImg =... 
            insertObjectAnnotation(img,'rectangle',bboxes,'vix',...
            'fontsize', 50,...
            'LineWidth', 5,...
            'color', 'yellow');
        ptime(1, ii) = toc;
        figure(ii);
        hold on;
        title('YOLOv2', 'fontsize', 18);
        imshow(detectedImg);
        hold off;
    end
    for ii = 1:1:4
        fprintf("processing time, img_%d  :  %0.3f\n", ii, ptime(1,ii));
    end
end


%------------------------------------------------------------------------------------------
if select == 14
    detector = vehicleDetectorACF();
    ptime = zeros(1, 4);
    for ii = 1:1:4
        if ii == 1
            img = imread('s1_dashcam-screen-shot1.jpeg');
        end
        if ii == 2
            img = imread('s2_pedstop.jpeg');
        end
        if ii == 3
            img = imread('s3_jogger.jpeg');
        end
        if ii == 4
            img = imread('s4_four-way-stop-sign-one-way-american-suburb.jpeg');
        end
        tic;
        [bboxes, scores] = detect(detector, img);
        detectedImg =... 
            insertObjectAnnotation(img,'rectangle',bboxes,'vix',...
            'fontsize', 50,...
            'LineWidth', 5,...
            'color', 'cyan');
        ptime(1, ii) = toc;
        figure(ii);
        hold on;
        title('ACF', 'fontsize', 18);
        imshow(detectedImg);
        hold off;
    end
    for ii = 1:1:4
        fprintf("processing time, img_%d  :  %0.3f\n", ii, ptime(1,ii));
    end
end


%------------------------------------------------------------------------------------------
if select == 15
    detector = vehicleDetectorFasterRCNN();
    ptime = zeros(1, 4);
    for ii = 1:1:4
        if ii == 1
            img = imread('s1_dashcam-screen-shot1.jpeg');
        end
        if ii == 2
            img = imread('s2_pedstop.jpeg');
        end
        if ii == 3
            img = imread('s3_jogger.jpeg');
        end
        if ii == 4
            img = imread('s4_four-way-stop-sign-one-way-american-suburb.jpeg');
        end
        tic;
        [bboxes, scores] = detect(detector, img);
        detectedImg =... 
            insertObjectAnnotation(img,'rectangle',bboxes,'vix',...
            'fontsize', 50,...
            'LineWidth', 5,...
            'color', 'green');
        ptime(1, ii) = toc;
        figure(ii);
        hold on;
        title('Fast RCNN', 'fontsize', 18);
        imshow(detectedImg);
        hold off;
    end
    for ii = 1:1:4
        fprintf("processing time, img_%d  :  %0.3f\n", ii, ptime(1,ii));
    end
end


%------------------------------------------------------------------------------------------
if select == 16
    csv_a = readtable("./db_lisa_tiny/annotations.csv");
    blds = boxLabelDatastore(csv_a);
    imDir = fullfile("./db_lisa_tiny");
    addpath(imDir);
    negativeFolder = fullfile(matlabroot,'toolbox','vision','visiondata','nonStopSigns');
    negativeImages = imageDatastore(negativeFolder);
    trainCascadeObjectDetector('stopSignDetector.xml', positiveInstances, negativeFolder,...
        FalseAlarmRate=0.1, NumCascadeStages=30, FeatureType='HOG');
    detector = vision.CascadeObjectDetector('stopSignDetector.xml', MergeThreshold=9);

    ptime = zeros(1, 4);
    for ii = 1:1:4
        if ii == 1
            img = imread('s1_dashcam-screen-shot1.jpeg');
        end
        if ii == 2
            img = imread('s2_pedstop.jpeg');
        end
        if ii == 3
            img = imread('s3_jogger.jpeg');
        end
        if ii == 4
            img = imread('s4_four-way-stop-sign-one-way-american-suburb.jpeg');
        end
        tic;
        bbox = step(detector,img);
        detectedImg = insertObjectAnnotation(img,'rectangle',bbox,'stop sign');
        ptime(1, ii) = toc;
        figure(ii);
        hold on;
        title('HOG', 'fontsize', 18);
        imshow(detectedImg);
        hold off;
        rmpath(imDir);
    end
    for ii = 1:1:4
        fprintf("processing time, img_%d  :  %0.3f\n", ii, ptime(1,ii));
    end
    %
end


%------------------------------------------------------------------------------------------
if select == 97
%Viola Jones Face Detector
    % read image
    I = imread ('visionteam.jpg');
    MergeThreshold = 4;
    
    % create Viola Jones detector
    VJ = vision.CascadeObjectDetector ('MergeThreshold', MergeThreshold);

    % detect faces
    tic
    bboxes = VJ(I);
    % time detector
    t_VJ = toc

    % annotate detected faces
    Ifaces = insertObjectAnnotation (I, 'rectangle', bboxes , 'Face');
    figure;
    imshow (Ifaces , []);
    title ('Viola Jones Face Detector');
    
    % clean up
    release(VJ);
end


%------------------------------------------------------------------------------------------
if select == 98
%HOG detector
    % read image
    I = imread('visionteam.jpg');
    
    % increase value when there are too many false positives
    ClassificationThreshold = 0;
    
    % adjust to group multiple responses
    MergeDetections = false ;

    HOG = vision.PeopleDetector( 'ClassificationThreshold',...
        ClassificationThreshold , 'MergeDetections', MergeDetections);

    tic
    [bboxes , scores] = HOG(I);
    t_HOG = toc

    Iped = insertObjectAnnotation (I,'rectangle', bboxes , scores);

    figure ;
    imshow (Iped , []);
    title ('HOG Pedestrian Detector');
    % notice detections contain area around the person (bounding boxes are not tight)
    
    % cleanup
    release (HOG);
end


%------------------------------------------------------------------------------------------
if select == 99
    
end
%%%%%%%%~~~~~~~~END>