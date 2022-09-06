%{
    chapter 11, image processing
    
    1  : ex 11.1, RGB->Y
    2  : ex 11.2, white balance
    3  : ex 11.3, gamma table
    4  : ex 11.4, histogram and equalize
    5  : ex 11.5, make 2D filters in units
    6  : ex 11.6, 2D filter in frequency domain ...needs adjustment
    7  : ex 11.7, 1D DCT
    8  : ex 11.8, DWT
    9  : hw1, klt11.4, imhist() on gamma corrected image from ex11.3
    10 : hw2a, use ex11.4 and hisogram equalize
    11 : hw2b, use blockproc() with 32x32 for equalization
    12 : hw3a, add white noise and filter
    13 : hw3b, add salt/pepper and filter
    14 : hw3c, filter a noisy image
    15 : hw4a, 7x7 and gauss
    16 : hw4b, fspecial()
    17 : hw5a/b/c/d/e
    18 : hw6, DCT using ex11.7
    
    bilinear()
    blockproc()
    butter()
    buttord()
    cat()
    colormap()
    cfirpm()
    decimate()
    dfilt.dffir
    ellip()
    fdatool
    fclose()
    fft()
    fft2()
    fftshift()
    filter()
    filterDesigner
    filtfilt()
    fir1()
    fir1s()
    fir2()
    fircls()
    fircls1()
    firpm()
    firpmord()
    fopen()
    fread()
    freqz()
    fspecial()
    fvtool
    gcd()
    grpdelay()
    hamming()
    histeq()
    ifft()
    ifft2()
    iirpeak()
    image() , use with figure
    imagesc()
    imfilter()
    imhist()
    imread()
    imshow()
    interp()
    kaiserord()
    mean()
    mean2()
    padarray()
    rand()
    remez()
    reshape()
    residuez()
    rgb2ycbcr()
    sosfilt()
    std()
    std2()
    tic
    toc
    tf2zp()
    uint8()
    upfirdn() 
    wintool
    wvtool
    ycbcr2rgb()
    zplane()
    
%}
close all;
clear all;
clc;
select = 99;  %        <-----CHANGE


%------------------------------------------------------------------------------------------
if (select == 1)
    % Read in R,G, and B data
    fid = fopen('Sunset800x600R8.RGB','r');
    R=fread(fid);
    fclose(fid);
    fid = fopen('Sunset800x600G8.RGB','r');
    G=fread(fid);
    fclose(fid);
    fid = fopen('Sunset800x600B8.RGB','r');
    B=fread(fid);
    fclose(fid);

    % Form the RGB color space
    RGB = cat(3, reshape(R, 800, 600)', ...
                 reshape(G, 800, 600)', ...
                 reshape(B, 800, 600)');

    % Convert double to uint8
    RGB = uint8(RGB);

    % Show original RGB image
    figure, imshow(RGB),  title('Origianl RGB Image');

    % RGB to YCbCr conversion
    YCbCr = rgb2ycbcr(RGB);

    % Show Y component of YCbCr data
    figure, imshow(YCbCr(:,:,1)),  title('Y of the YCbCr Image');

    % Write 8-bit data out to files
    fid=fopen('Sunset800x600Y8.dat', 'wb');
    fwrite(fid, YCbCr(:,:,1)', 'uint8');  % Write Y
    fclose(fid);
    fid=fopen('Sunset800x600Cb8.dat', 'wb');
    fwrite(fid, YCbCr(:,:,2)', 'uint8');  % Write Cb
    fclose(fid);
    fid=fopen('Sunset800x600Cr8.dat', 'wb');
    fwrite(fid, YCbCr(:,:,3)', 'uint8');  % Write Cr
    fclose(fid);

    % Read back 8-bit YCbCr data 
    fid=fopen('Sunset800x600Y8.dat', 'rb');
    Y = fread(fid, inf, 'uint8');		% Read Y
    fclose(fid);
    fid=fopen('Sunset800x600Cb8.dat', 'rb');
    Cb = fread(fid, inf, 'uint8');		% Read Cb
    fclose(fid);
    fid=fopen('Sunset800x600Cr8.dat', 'rb');
    Cr = fread(fid, inf, 'uint8');		% Read Cr
    fclose(fid);

    Y = Y/255;
    Cb = Cb/255;
    Cr = Cr/255;
    YCbCr1 = cat(3, reshape(Y, 800, 600)', ...
                    reshape(Cb, 800, 600)', ...
                    reshape(Cr, 800, 600)');
    RGB1 = ycbcr2rgb(YCbCr1);
    RGB1 = uint8(RGB1*255);

    % Show images
    figure, imshow(RGB1), title('The RGB->YCbCr->RGB Converted Image');
end


%------------------------------------------------------------------------------------------
if (select == 2)
    RGB = imread('example11.2.jpg');
    [height, width, color] = size(RGB);
    figure, imshow(RGB); title('Origianl JPEG Image');
    R = sum(sum(RGB(:,:,1)));         % Compute sum of R
    G = sum(sum(RGB(:,:,2)));         % Compute sum of G
    B = sum(sum(RGB(:,:,3)));         % Compute sum of B
    gr = 1;
    gg = 1;
    gb = 1;
    if G<R
        if G<B   
            gr = G/R;                 % Normalized gain factor for R 
            gb = G/B;                 % Normalized gain factor for B 
        end
    end
    if R<G
        if R<B   
            gg = R/G;                 % Normalized gain factor for G 
            gb = R/B;                 % Normalized gain factor for B 
        end
    end
    if B<R
        if B<G
            gr = B/R;                 % Normalized gain factor for R 
            gg = B/G;                 % Normalized gain factor for G 
        end
    end
    Rw=RGB(:,:,1)*double(gr);		 % Apply normalized gain factor to R
    Gw=RGB(:,:,2)*double(gg);		 % Apply normalized gain factor to G 
    Bw=RGB(:,:,3)*double(gb);		 % Apply normalized gain factor to B
    RGBw(:,:,1) = Rw;
    RGBw(:,:,2) = Gw;
    RGBw(:,:,3) = Bw;
    figure, imshow(RGBw), title('White Balanced RGB Image');
end


%------------------------------------------------------------------------------------------
if (select == 3)
    RGB = imread('example11.3.bmp');   % Read the image in
    figure(1), image(RGB);             % Display the original image
    title('Origianl BMP Image');       

    [height, width, color]=size(RGB);
    gamma=2.2;                         % Default to 2.20 for Window PC

    x=0:255;                           % 256 linear input data
    y=power(x, 1/gamma);               % Generate x^1/gamma
    gammaTable = floor(255*y/y(256));  % Normalize and make 8-bit gamma table

    figure(2), plot(gammaTable); grid; % Display gamma curve
    title('8-bit gamma curve');

    fid=fopen('gammaTable.txt', 'w');  % Open file to write

    fprintf(fid, 'char gammaTable[256]={%4d,\n\r',gammaTable(1));
    for i=2:255
       fprintf(fid,'%4d,\n\r',gammaTable(i));
    end
    fprintf(fid,'%4d};\n\r',gammaTable(256));

    fclose(fid);                      % Close file


    R = RGB(:,:,1);                   % Get the R, G, andB from RGB array
    G = RGB(:,:,2);
    B = RGB(:,:,3);

    for i=1:height                    % Apply gamma correction to R, G, and B
        for j=1:width
          if R(i,j) > 0
             Rgamma(i,j) = gammaTable(R(i,j));
          else
             Rgamma(i,j) = 0;
          end

          if G(i,j) > 0
             Ggamma(i,j) = gammaTable(G(i,j));
          else
             Ggamma(i,j) = 0;
          end

          if B(i,j) > 0
             Bgamma(i,j) = gammaTable(B(i,j));
          else
             Bgamma(i,j) = 0;
          end
        end
    end

    RGBgamma(:,:,1) = Rgamma;
    RGBgamma(:,:,2) = Ggamma;
    RGBgamma(:,:,3) = Bgamma;
    RGBgamma = uint8(RGBgamma);              
    figure(3), image(RGBgamma);       % Display gamma corrected image 
    title('Gamma Corrected RGB Image');
end


%------------------------------------------------------------------------------------------
if (select == 4)
    RGB = imread('example11.4.jpg');    % Read the image in
    [height, width, color] = size(RGB); % Get the width and height
    len = 256;

    YCbCr = rgb2ycbcr(RGB);             % RGB to YCbCr conversion

    Y = double(YCbCr(:,:,1));           % Get Y only

    hist = zeros(1,len);
    hist2 = imhist(uint8(Y));           % Use MATLAB built-in histogram func

    for i=1:height                      % Compute histogram of image
        for j=1:width              
            index = uint8((Y(i,j))+1);
            hist(index) = hist(index) + 1;
        end
    end

    imageSize = width * height;         % Create equalize look up table
    eqFactor = 255 / imageSize;         
    sum = 0;
    for i=1:len
        sum = sum + hist(i);
        eqValue = sum * eqFactor;
        if eqValue > 255;               % Limit upper data
            eqValue = 255;
        end    
        eqTable(i) = uint8(eqValue);
    end

    for i=1:height                      % Image equalization, new Y       
        for j=1:width              
            index = uint8((Y(i,j))+1);
            newY(i,j) = eqTable(index);
        end
    end

    newY2 = histeq(uint8(Y),256);       % Use MATLAB built-in histeq funciton

    newYCbCr(:,:,1) = newY;
    newYCbCr(:,:,2) = YCbCr(:,:,2);
    newYCbCr(:,:,3) = YCbCr(:,:,3);
    newYCbCr2(:,:,1) = newY2;
    newYCbCr2(:,:,2) = YCbCr(:,:,2);
    newYCbCr2(:,:,3) = YCbCr(:,:,3);

    newRGB  = ycbcr2rgb(newYCbCr);
    newRGB2 = ycbcr2rgb(newYCbCr2);

    Y   = double(newYCbCr(:,:,1));     % Example result
    Y2  = double(newYCbCr2(:,:,1));    % Result from MATLAB function

    newHist  = zeros(1,len);
    newHist2  = zeros(1,len);

    for i=1:height                     % Compute new image histogram
        for j=1:width              
            index = uint8((Y(i,j))+1);
            newHist(index) = newHist(index) + 1;
            index = uint8((Y2(i,j))+1);
            newHist2(index) = newHist2(index) + 1;        
        end
    end

    figure(1), image(RGB);             % Display the original image
    title('Origianl JPEG Image');      

    figure(2), image(newRGB);          % Display equalized image 
    title('Equalized Image');

    figure(3), image(newRGB2);         % Display MATLAB equalized image 
    title('Equalized Image Using MATLAB histeq');


    x = 0:len-1;

    figure(4), subplot(2,1,1);
    imhist(YCbCr(:,:,1));
    title('Histograms of original image');

    figure(4), subplot(2,1,2);
    imhist(newY); 
    title('Histograms equalized image');

    figure(5), subplot(2,1,1);
    imhist(newY2); 
    title('MATLAB histeq equalized');
end


%------------------------------------------------------------------------------------------
if (select == 5)
%           'delta'     for a delat filter has no effect to the image
%         'average'   for an averaging filter
%         'sharp'     for a sharp contrast enhancement filter
%         'sobel'     for a Sobel horizontal edge-emphasizing filter 
%         'prewitt'   for a Prewitt vertical edge-emphasizing filter 
%         'laplacian' for a filter approximating 2-D Laplacian operator
%         'emboss'    for a filter approximating 3-D effect
%

    delta_filt = [0  0  0; ...         % Delat filter
                  0  1  0; ...
                  0  0  0];            

    lowpass_filt = [1  1  1; ...       % Low pass (average) filter
                    1  1  1; ...
                    1  1  1]/9;          

    highpass_filt = [-1  -4 -1; ...    % High pass (sharp) filter
                     -4  26 -4; ...
                     -1  -4 -1]/6; 

    sobel_filt = [1   2   1; ...       % Sobel filter (horizontal arrangement)
                  0   0   0; ...
                 -1  -2  -1];

    prewitt_filt = [1  0 -1; ...       % Prewitt filter (vertical arragenment)
                    1  0 -1; ...
                    1  0 -1];          

    laplacian_filt = [1   4   1; ...   % Laplacian filter
                      4 -20   4; ...
                      1   4   1];

    emboss_filt = [-4 -4   0; ...      % Modified diagnal emboss filter
                   -4  1   4; ...
                    0  4   4];

    blur_filt = [0   1   0; ...        % Blur filter
                 1   4   1; ...
                 0   1   0]/8;

    edge_filt = [-1  -1  -1; ...       % Edge filter (edge in all direction)
                 -1   8  -1; ...
                 -1  -1  -1];

    engrave_filt = [-1   0   0; ...    % Emboss and engrave filter
                     0   2   0; ...
                     0   0   0];

    RGB = imread('example11.5.jpg');   % Read the image in

    [imHeight, imWidth, color] = size(RGB); % Get the width and height

    figure(1), image(RGB);             % Display the original image
    title('Origianl JPEG Image');      

    for i=2:11
        if i == 2                      % Load filter
            filt = delta_filt;
        elseif i == 3
            filt = lowpass_filt;
        elseif i == 4
            filt = highpass_filt;
        elseif i == 5
            filt = sobel_filt;
        elseif i == 6
            filt = prewitt_filt;
        elseif i == 7
            filt = laplacian_filt;             
        elseif i == 8
            filt = emboss_filt;    
        elseif i == 9
            filt = blur_filt;
        elseif i == 10
            filt = edge_filt;      
        else
            filt = engrave_filt;
        end

    if 0 % Use 3 2-D filters
        R = filter2(filt, RGB(:,:,1)); % Filter the image
        G = filter2(filt, RGB(:,:,2));
        B = filter2(filt, RGB(:,:,3));
        newRGB(:,:,1) = uint8(R);
        newRGB(:,:,2) = uint8(G);
        newRGB(:,:,3) = uint8(B);
    else % Use image-filter
        newRGB = imfilter(RGB, filt);  % Filter the image
    end

        figure(i), image(newRGB);      % Display filtered image 

        if i == 2                      % Find which filter was used
            title('Delta Filtered Image');
        elseif i == 3
            title('Lowpass Filtered Image');
        elseif i == 4
            title('Highpass Filtered Image');
        elseif i == 5
            title('Sobel H-Filtered Image');
        elseif i == 6
            title('Prewitt V-Filtered Image');
        elseif i == 7
            title('Laplacian Filtered Image');
        elseif i == 8
            title('Emboss Filtered Image');  
        elseif i == 9
            title('Blur Filtered Image');      
        elseif i == 10
            title('Edge Filtered Image');      
        else        
            title('Engrave Filtered Image');                 
        end
    end
end


%------------------------------------------------------------------------------------------
if (select == 6)
    delta_filt = [0  0  0  0  0  0  0  0  0; ...  % Delat filter
                  0  0  0  0  0  0  0  0  0; ...
                  0  0  0  0  0  0  0  0  0; ...
                  0  0  0  0  0  0  0  0  0; ...
                  0  0  0  0  1  0  0  0  0; ...
                  0  0  0  0  0  0  0  0  0; ...
                  0  0  0  0  0  0  0  0  0; ...
                  0  0  0  0  0  0  0  0  0; ...
                  0  0  0  0  0  0  0  0  0];            

    edge_filt = [-1 -1 -1 -1 -1 -1 -1 -1 -1; ... % Edge filter
                 -1 -1 -1 -1 -1 -1 -1 -1 -1; ...
                 -1 -1  1  1  1  1  1 -1 -1; ...
                 -1 -1  1  1  1  1  1 -1 -1; ...
                 -1 -1  1  1 31  1  1 -1 -1; ...
                 -1 -1  1  1  1  1  1 -1 -1; ...
                 -1 -1  1  1  1  1  1 -1 -1; ...
                 -1 -1 -1 -1 -1 -1 -1 -1 -1; ...
                 -1 -1 -1 -1 -1 -1 -1 -1 -1];                          

    motion_filt = [1  0  0  0  0  0  0  0  0; ... % Motion filter
                   0  1  0  0  0  0  0  0  0; ...
                   0  0  1  0  0  0  0  0  0; ...
                   0  0  0  1  0  0  0  0  0; ...
                   0  0  0  0  1  0  0  0  0; ...
                   0  0  0  0  0  1  0  0  0; ...
                   0  0  0  0  0  0  1  0  0; ...
                   0  0  0  0  0  0  0  1  0; ...
                   0  0  0  0  0  0  0  0  1]/9;                          

    gaussian_filt = [0  0  0  0  0  0  0  0  0; ... % Gaussian filter
                     0  0  0  0  0  0  0  0  0; ...
                     0  0  1  4  6  4  1  0  0; ...
                     0  0  4 16 24 16  4  0  0; ...
                     0  0  6 24 36 24  6  0  0; ...
                     0  0  4 16 24 16  4  0  0; ...
                     0  0  1  4  6  4  1  0  0; ...
                     0  0  0  0  0  0  0  0  0; ...
                     0  0  0  0  0  0  0  0  0]/256;                          


    RGB1 = imread('example11.6.jpg');   % Read the image in
    RGB=padarray(RGB1,[8,8],0, 'post'); % Padding 8 rows and 8 columns zeros at the end
    fft2R = fft2(double(RGB(:,:,1)));  % 2-D FFT of R component
    fft2G = fft2(double(RGB(:,:,2)));  % 2-D FFT of G component
    fft2B = fft2(double(RGB(:,:,3)));  % 2-D FFT of B component

    [imHeight imWidth] = size(fft2R);  % Find image width and height

    figure(1), image(RGB);             % Display the original image
    title('Origianl JPEG Image');      

    for i=2:5
        if i == 2                      % Load filter
            filt = delta_filt;
        elseif i == 3
            filt = edge_filt;
        elseif i == 4
            filt = motion_filt;
        else
            filt = gaussian_filt;
        end

        fft2Filt = fft2(filt, imHeight, imWidth);  % 2-D FFT of filter kernel

        fft2FiltR = fft2Filt .* fft2R;             % 2-D fast convolution
        fft2FiltG = fft2Filt .* fft2G;
        fft2FiltB = fft2Filt .* fft2B;

        newRGB(:,:,1) = uint8(ifft2(fft2FiltR));   % Inverse 2-D FFT of R component
        newRGB(:,:,2) = uint8(ifft2(fft2FiltG));   % Inverse 2-D FFT of G component
        newRGB(:,:,3) = uint8(ifft2(fft2FiltB));   % Inverse 2-D FFT of B component

        figure(i), image(newRGB);                  % Display filtered image 

        if i == 2                                  % Find which filter was used
            title('Delta Filtered Image');
        elseif i == 3
            title('edge Filtered Image');
        elseif i == 4
            title('Motion Filtered Image');
        else
            title('Gaussian Filtered Image');
        end
    end
end


%------------------------------------------------------------------------------------------
if (select == 7)
    RGB = imread('example11.7.jpg');            % Read the image in
    fid = fopen('dctFile.dat', 'w');
    [imHeight imwidth color] = size(RGB);       % Get image dimensions

    YCbCr = rgb2ycbcr(RGB);                     % Converts to YCbCr color space
    Y = double(YCbCr(:,:,1));                   % Change to YUV space 
    U = double(YCbCr(:,:,2))-128;
    V = double(YCbCr(:,:,3))-128;

    for n=1:8:imHeight                         % Each block is 8 pixels in height
      for m=1:8:imwidth                        % Each block is 8 pixels in width
        for i=0:7
          for j=0:7
            mbY(i+1,j+1) = Y(n+i,m+j,1);       % Form 8x8 micro-block for Y          
            mbU(i+1,j+1) = U(n+i,m+j,1);       % Form 8x8 micro-block for U     
            mbV(i+1,j+1) = V(n+i,m+j,1);       % Form 8x8 micro-block for V     
          end
        end

        mbY = dct(double(mbY));                % Perform 1-D DCT horizontally          
        mbY = dct(double(mbY'));               % Perform 1-D DCT vertically          
        mbU = dct(double(mbU));                % Perform 1-D DCT horizontally          
        mbU = dct(double(mbU'));               % Perform 1-D DCT vertically          
        mbV = dct(double(mbV));                % Perform 1-D DCT horizontally          
        mbV = dct(double(mbV'));               % Perform 1-D DCT vertically          
        oY = int16(mbY);  
        oU = int16(mbU);  
        oV = int16(mbV);      
        fwrite(fid, oY, 'int16');  
        fwrite(fid, oU, 'int16');  
        fwrite(fid, oV, 'int16');  

      end
    end

    fclose(fid);                              % Close DCT data file

    %%%%%%%%%%%%%%%%%%%%%%%%%%%% iDCT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fid = fopen('dctFile.dat', 'r');          % Open DCT data file to read
    dctData = fread(fid, inf, 'int16');       % Read DCT data in
    dctData = int16(dctData);

    k = 1;
    for n=1:8:imHeight
      for m=1:8:imwidth

        for j=1:8
          for i=1:8
             imbY(i,j) = dctData(k);
             k = k+1;
          end
        end
        for j=1:8
          for i=1:8
            imbU(i,j) = dctData(k);
            k = k+1;
          end
        end
        for j=1:8
          for i=1:8
           imbV(i,j) = dctData(k);
           k = k+1;
          end
        end

        imbY = idct(double(imbY));            % 1-D iDCT horizontally
        imbY = idct(double(imbY'));           % 1-D iDCT vertically
        imbU = idct(double(imbU));            % 1-D iDCT horizontally
        imbU = idct(double(imbU'));           % 1-D iDCT vertically
        imbV = idct(double(imbV));            % 1-D iDCT horizontally
        imbV = idct(double(imbV'));           % 1-D iDCT vertically

        for i=0:7
          for j=0:7
            iY(n+i,m+j) = imbY(i+1,j+1);
            iU(n+i,m+j) = imbU(i+1,j+1);
            iV(n+i,m+j) = imbV(i+1,j+1);
          end
        end

      end
    end

    newYCbCr(:,:,1) = uint8(round(iY));      % Form YCbCr data format  
    newYCbCr(:,:,2) = uint8(round(iU) + 128);
    newYCbCr(:,:,3) = uint8(round(iV) + 128);

    newRGB = ycbcr2rgb(newYCbCr);            % Convert YCbCr to RGB color space

    err = RGB - newRGB;

    figure (1), imshow(RGB);                 % Show the original image
    title('The original image');

    figure (2), imshow(newRGB);              % Show the reconstructed image
    title('The DCT/IDCT reconstructed image');

    figure (3), plot(err(:,:,1), 'red');	% Show the error in R
    title('The error of R');
    figure (4), plot(err(:,:,2), 'green');	% Show the error in G
    title('The error of G');
    figure (5), plot(err(:,:,3), 'blue');	% Show the error in B
    title('The error of B');

    fclose(fid);							% Close DCT data file
    delete .\\dctFile.dat;
end


%------------------------------------------------------------------------------------------
if (select == 8)
    [X,map] = imread('example11.8.bmp'); % Read the image in

    % Map contains the loaded colormap. 
    nbcol = size(map,1);

    % generate Decomposition and Reconstruction filter coefficients
    [L_D,H_D,L_R,H_R] = wfilters('db2'); % Lowpass decomposite, highpass decomposite
                                         %   lowpass reconstruction, and highpass reconstruction

    % Perform single-level decomposition 
    [LL,LH,HL,HH] = dwt2(X,L_D,H_D);	% 2-D DWT

    % Images coding
    dec2d = [... 
            LL,     LH;     ... 
            HL,     HH      ... 
            ];
    % Construction filter
    LL = zeros(size(LL));              % Zero out the LL

    % Reconstruct the image
    Y = idwt2(LL,LH,HL,HH,L_R,H_R);    % Inverse 2-D DWT

    % Make the image B&W
    [height, width] = size(X);
    for i=1:height
      for j=1:width
        if Y(i,j) < 42
          Y(i,j) = 10;
        else
          Y(i,j) = 240;
        end
      end
    end

    % Plotting
    figure;
    clf
    imagesc(X);
    axis equal
    axis off
    colormap(gray(128));
    title('The original image');

    figure;
    imagesc(dec2d) ;
    axis equal
    axis off
    colormap(gray(128));
    title('Decomposed image');

    figure;
    imagesc(Y); 
    axis equal
    axis off
    colormap(gray(128));
    title('Reconstructed image without the LL subband');
end


%------------------------------------------------------------------------------------------
if (select == 9)  
    RGB = imread('example11.3.bmp');   % Read the image in
    figure(1), image(RGB);             % Display the original image
    title('Original BMP Image');       
    [height, width, color]=size(RGB);
    
    figure(2), imhist(RGB);            % histogram of original image
    title("from original image"); 
    gamma=2.2;                         % Default to 2.20 for Window PC

    x=0:255;                           % 256 linear input data
    y=power(x, 1/gamma);               % Generate x^1/gamma
    gammaTable = floor(255*y/y(256));  % Normalize and make 8-bit gamma table

    fid=fopen('gammaTable.txt', 'w');  % Open file to write

    fprintf(fid, 'char gammaTable[256]={%4d,\n\r',gammaTable(1));
    for i=2:255
       fprintf(fid,'%4d,\n\r',gammaTable(i));
    end
    fprintf(fid,'%4d};\n\r',gammaTable(256));

    fclose(fid);                      % Close file

    R = RGB(:,:,1);                   % Get the R, G, andB from RGB array
    G = RGB(:,:,2);
    B = RGB(:,:,3);

    for i=1:height                    % Apply gamma correction to R, G, and B
        for j=1:width
          if R(i,j) > 0
             Rgamma(i,j) = gammaTable(R(i,j));
          else
             Rgamma(i,j) = 0;
          end

          if G(i,j) > 0
             Ggamma(i,j) = gammaTable(G(i,j));
          else
             Ggamma(i,j) = 0;
          end

          if B(i,j) > 0
             Bgamma(i,j) = gammaTable(B(i,j));
          else
             Bgamma(i,j) = 0;
          end
        end
    end

    RGBgamma(:,:,1) = Rgamma;
    RGBgamma(:,:,2) = Ggamma;
    RGBgamma(:,:,3) = Bgamma;
    RGBgamma = uint8(RGBgamma);              
    figure(3), image(RGBgamma);       % Display gamma corrected image 
    title('Gamma Corrected RGB Image');
    figure(4), imhist(RGBgamma);       % histogram of gamma corrected image
    title("from gamma corrected image");
end


%------------------------------------------------------------------------------------------
if (select == 10)
    RGB = imread('jetplane.png');    % Read the image in
    [height, width, color] = size(RGB); % Get the width and height
    %len = 256;
    %len = 128;
    len = 64;

    YCbCr = rgb2ycbcr(RGB);           % RGB to YCbCr to isolate intensity
    Y = double(YCbCr(:,:,1));         % just need intensity
    hist2 = imhist(uint8(Y));         % make a histogram
    newY2 = histeq(uint8(Y), len);    % equalize it

    newYCbCr2(:,:,1) = newY2;           % equalized intensity
    newYCbCr2(:,:,2) = YCbCr(:,:,2);    % old values
    newYCbCr2(:,:,3) = YCbCr(:,:,3);    % old values

    newRGB2 = ycbcr2rgb(newYCbCr2);    % prepare image for display 

    Y2  = double(newYCbCr2(:,:,1));    % Result from MATLAB function

    figure(1), image(RGB);             % Display the original image
    title('Origianl JPEG Image');      

    figure(2), image(newRGB2);         % Display MATLAB equalized image 
    title('Equalized Image Using MATLAB histeq');

    figure(3), subplot(1,2,1);
    imhist(YCbCr(:,:,1));
    title('Histograms of original image');

    subplot(1,2,2);
    imhist(newY2); 
    title('MATLAB histeq equalized');
end


%------------------------------------------------------------------------------------------
if (select == 11)
    RGB = imread('jetplane.png');    % Read the image in
    figure(1), image(RGB);           % Display the original image
    title('Origianl JPEG Image');   
    
    my_fun = @(block_struct) histeq(block_struct.data);
    new_image = blockproc(RGB, [32, 32], my_fun);
    
    figure(2); image(new_image);
    title('after blockproc()');
    
    %{
    RGB = imread('jetplane.png');    % Read the image in
    figure(1), image(RGB);           % Display the original image
    title('Origianl JPEG Image');   
    
    YCbCr = rgb2ycbcr(RGB);           % RGB to YCbCr to isolate intensity
    
    my_fun = @(block_struct) histeq(block_struct.data);
    new_Y = blockproc(YCbCr, [32, 32], my_fun);
    
    new_image = ycbcr2rgb(new_Y);       % prepare image for display 
    new_image(:,:,2) = YCbCr(:,:,2);    % old values
    new_image(:,:,3) = YCbCr(:,:,3);    % old values
    figure(2); image(new_image);
    title('after blockproc()');
    %}
end


%------------------------------------------------------------------------------------------
if (select == 12)
    RGB = imread("DSCN0479-001.jpeg");
    [height, width, color] = size(RGB);
    figure(1); image(RGB); title("original image");
    
    YCbCr = rgb2ycbcr(RGB);           % RGB to YCbCr to isolate intensity
    Y = double(YCbCr(:,:,1));         % just need intensity
    
    dirty_Y = zeros(height, width);
    for ii = 1:1:height
        for jj = 1:1:width
            temp = Y(ii, jj) + normrnd(0, sqrt(0.005));
            if temp > 200
                dirty_Y(ii,jj) = 200;
            elseif temp < 0
                dirty_Y(ii,jj) = 0;
            else
                dirty_Y(ii,jj) = temp;
            end
        end
    end

    dirty_YCbCr(:,:,1) = uint8(dirty_Y);           % add the noise
    dirty_YCbCr(:,:,2) = YCbCr(:,:,2);    % old values
    dirty_YCbCr(:,:,3) = YCbCr(:,:,3);    % old values
    dirty_RGB = ycbcr2rgb(dirty_YCbCr);
    figure(2); image(dirty_RGB); title("noisy image");
    
    fprintf("\n3x3 smooth filter running time:\n");
    tic
    h_33 = (1/9) .* ones(3, 3);
    f33_RGB = imfilter(dirty_RGB, h_33);
    figure(3); image(f33_RGB); title("3x3 smooth");
    toc
    YCbCr_var = rgb2ycbcr(f33_RGB);
    Y_var = double(YCbCr_var(:, :, 1));
    mse = 0;
    for ii = 1:1:height
        for jj = 1:1:width
            mse =  mse + (Y(ii, jj) - Y_var(ii, jj))^2;
        end
    end
    mse = mse / (height * width);
    fprintf("MSE of 3x3 smooth filter:  %0.4f\n", mse);
    fprintf("PSNR =  %0.4f\n", 20 * log10(255/sqrt(mse)));
    
    fprintf("\n7x7 smooth filter running time:\n");
    tic
    h_77 = (1/49) .* ones(7, 7);
    f77_RGB = imfilter(dirty_RGB, h_77);
    figure(4); image(f77_RGB); title("7x7 smooth");
    toc
    YCbCr_var = rgb2ycbcr(f77_RGB);
    Y_var = double(YCbCr_var(:, :, 1));
    mse = 0;
    for ii = 1:1:height
        for jj = 1:1:width
            mse =  mse + (Y(ii, jj) - Y_var(ii, jj))^2;
        end
    end
    mse = mse / (height * width);
    fprintf("MSE of 7x7 smooth filter:  %0.4f\n", mse);
    fprintf("PSNR =  %0.4f\n", 20 * log10(255/sqrt(mse)));
    
    fprintf("\nmedian filter running time:\n");
    tic
    fmed = medfilt2(dirty_Y);
    fmed_YCbCr(:,:,1) = uint8(fmed);
    fmed_YCbCr(:,:,2) = YCbCr(:,:,2);    % old values
    fmed_YCbCr(:,:,3) = YCbCr(:,:,3);    % old values
    fmed_RGB = ycbcr2rgb(dirty_YCbCr);
    figure(5); image(fmed_RGB); title("median filter");
    toc
    YCbCr_var = rgb2ycbcr(fmed_RGB);
    Y_var = double(YCbCr_var(:, :, 1));
    mse = 0;
    for ii = 1:1:height
        for jj = 1:1:width
            mse =  mse + (Y(ii, jj) - Y_var(ii, jj))^2;
        end
    end
    mse = mse / (height * width);
    fprintf("MSE of median filter:  %0.4f\n", mse);
    fprintf("PSNR =  %0.4f\n", 20 * log10(255/sqrt(mse)));
end


%------------------------------------------------------------------------------------------
if (select == 13)
    RGB = imread("DSCN0479-001.jpeg");
    [height, width, color] = size(RGB);
    figure(1); image(RGB); title("original image");
    
    YCbCr = rgb2ycbcr(RGB);           % RGB to YCbCr to isolate intensity
    Y = double(YCbCr(:,:,1));         % just need intensity
 
    dirty_Y = zeros(height, width);
    for ii = 1:1:height
        for jj = 1:1:width
            temp = randi([0,19]); % density
            if temp == 1
                dirty_Y(ii, jj) = 2; % add pepper
            elseif temp == 2
                dirty_Y(ii,jj) = 200; % add salt
            else
                dirty_Y(ii,jj) = Y(ii,jj); % keep original
            end
        end
    end

    dirty_YCbCr(:,:,1) = uint8(dirty_Y);           % add the noise
    dirty_YCbCr(:,:,2) = YCbCr(:,:,2);    % old values
    dirty_YCbCr(:,:,3) = YCbCr(:,:,3);    % old values
    dirty_RGB = ycbcr2rgb(dirty_YCbCr);
    figure(2); image(dirty_RGB); title("noisy image");
    
    fprintf("\n3x3 smooth filter running time:\n");
    tic
    h_33 = (1/9) .* ones(3, 3);
    f33_RGB = imfilter(dirty_RGB, h_33);
    figure(3); image(f33_RGB); title("3x3 smooth");
    toc
    YCbCr_var = rgb2ycbcr(f33_RGB);
    Y_var = double(YCbCr_var(:, :, 1));
    mse = 0;
    for ii = 1:1:height
        for jj = 1:1:width
            mse =  mse + (Y(ii, jj) - Y_var(ii, jj))^2;
        end
    end
    mse = mse / (height * width);
    fprintf("MSE of 3x3 smooth filter:  %0.4f\n", mse);
    fprintf("PSNR =  %0.4f\n", 20 * log10(255/sqrt(mse)));
    
    fprintf("\n7x7 smooth filter running time:\n");
    tic
    h_77 = (1/49) .* ones(7, 7);
    f77_RGB = imfilter(dirty_RGB, h_77);
    figure(4); image(f77_RGB); title("7x7 smooth");
    toc
    YCbCr_var = rgb2ycbcr(f77_RGB);
    Y_var = double(YCbCr_var(:, :, 1));
    mse = 0;
    for ii = 1:1:height
        for jj = 1:1:width
            mse =  mse + (Y(ii, jj) - Y_var(ii, jj))^2;
        end
    end
    mse = mse / (height * width);
    fprintf("MSE of 7x7 smooth filter:  %0.4f\n", mse);
    fprintf("PSNR =  %0.4f\n", 20 * log10(255/sqrt(mse)));
    
    fprintf("\nmedian filter running time:\n");
    tic
    fmed = medfilt2(dirty_Y);
    fmed_YCbCr(:,:,1) = uint8(fmed);
    fmed_YCbCr(:,:,2) = YCbCr(:,:,2);    % old values
    fmed_YCbCr(:,:,3) = YCbCr(:,:,3);    % old values
    fmed_RGB = ycbcr2rgb(dirty_YCbCr);
    figure(5); image(fmed_RGB); title("median filter");
    toc
    YCbCr_var = rgb2ycbcr(fmed_RGB);
    Y_var = double(YCbCr_var(:, :, 1));
    mse = 0;
    for ii = 1:1:height
        for jj = 1:1:width
            mse =  mse + (Y(ii, jj) - Y_var(ii, jj))^2;
        end
    end
    mse = mse / (height * width);
    fprintf("MSE of median filter:  %0.4f\n", mse);
    fprintf("PSNR =  %0.4f\n", 20 * log10(255/sqrt(mse)));
end


%------------------------------------------------------------------------------------------
if (select == 14)
    RGB = imread("DSCN0482-001.jpeg");
    [height, width, color] = size(RGB);
    figure(1); image(RGB); title("original image");
    
    YCbCr = rgb2ycbcr(RGB);           % RGB to YCbCr to isolate intensity
    Y = double(YCbCr(:,:,1));         % just need intensity
    
    fprintf("\n3x3 smooth filter running time:\n");
    tic
    h_33 = (1/9) .* ones(3, 3);
    f33_RGB = imfilter(RGB, h_33);
    figure(3); image(f33_RGB); title("3x3 smooth");
    toc
    
    fprintf("\n7x7 smooth filter running time:\n");
    tic
    h_77 = (1/49) .* ones(7, 7);
    f77_RGB = imfilter(RGB, h_77);
    figure(4); image(f77_RGB); title("7x7 smooth");
    toc
    
    fprintf("\nmedian filter running time:\n");
    tic
    fmed = medfilt2(Y);
    fmed_YCbCr(:,:,1) = uint8(fmed);
    fmed_YCbCr(:,:,2) = YCbCr(:,:,2);    % old values
    fmed_YCbCr(:,:,3) = YCbCr(:,:,3);    % old values
    fmed_RGB = ycbcr2rgb(YCbCr);
    figure(5); image(fmed_RGB); title("median filter");
    toc
    YCbCr_var = rgb2ycbcr(fmed_RGB);
    Y_var = double(YCbCr_var(:, :, 1));
end


%------------------------------------------------------------------------------------------
if (select == 15)
    RGB = imread("city.jpeg");
    [height, width, color] = size(RGB);
    YCbCr = rgb2ycbcr(RGB);           
    Y = double(YCbCr(:,:,1));         
    figure(1); image(RGB); title("original image");
    figure(2); image(RGB); title("original image");
    
    h_77 = (1/49) .* ones(7, 7);
    f77_RGB = imfilter(RGB, h_77);
    figure(3); image(f77_RGB); title("7x7 smooth");
    
    fgaus = imgaussfilt(Y, 0.5);
    fgaus_YCbCr(:, :, 1) = uint8(fgaus);
    fgaus_YCbCr(:, :, 2) = YCbCr(:, :, 2);
    fgaus_YCbCr(:, :, 3) = YCbCr(:, :, 3);
    fguas_RGB = ycbcr2rgb(fgaus_YCbCr);
    figure(4); image(fguas_RGB); title("gauss, sig = 0.5");
end


%------------------------------------------------------------------------------------------
if (select == 16)
    %{
    RGB = imread("city.jpeg");
    [height, width, color] = size(RGB);
    YCbCr = rgb2ycbcr(RGB);           
    Y = double(YCbCr(:,:,1));         
    figure(1); image(RGB); title("original image");
    figure(2); image(RGB); title("original image");
    
    h_sobel = fspecial('sobel');
    sobel_RGB = imfilter(RGB, h_sobel);
    figure(3); image(sobel_RGB); title("sobel filtered");
    
    h_lap = fspecial('laplacian', .1);
    lap_RGB = imfilter(RGB, h_lap);
    figure(4); image(lap_RGB); title("laplacian filtered");
    %}
    %
    RGB = imread("city.jpeg");
    [height, width, color] = size(RGB);
    YCbCr = rgb2ycbcr(RGB);           
    Y = double(YCbCr(:,:,1));         
    figure(1); image(RGB); title("original image");
    figure(2); image(RGB); title("original image");
    
    h_sobel = fspecial('sobel');
    Y_sobel = imfilter(Y, h_sobel);
    YCbCr_sobel(:, :, 1) = uint8(Y_sobel);
    YCbCr_sobel(:, :, 2) = YCbCr(:, :, 2);
    YCbCr_sobel(:, :, 3) = YCbCr(:, :, 3);
    RGB_sobel = ycbcr2rgb(YCbCr_sobel);
    figure(3); image(RGB_sobel); title("sobel filtered");
    
    h_lap = fspecial('laplacian', .2);
    Y_lap = imfilter(Y, h_lap);
    YCbCr_lap(:, :, 1) = uint8(Y_lap);
    YCbCr_lap(:, :, 2) = YCbCr(:, :, 2);
    YCbCr_lap(:, :, 3) = YCbCr(:, :, 3);
    RGB_lap = ycbcr2rgb(YCbCr_lap);
    figure(4); image(RGB_lap); title("laplacian filtered");
    %}
end


%------------------------------------------------------------------------------------------
if (select == 17)
    hitter = 6;
    RGB1 = imread('city.jpeg');
    [height, width, color] = size(RGB1);
    figure(1); image(RGB1); title('Origianl JPEG Image');
    RGB = padarray(RGB1,[height, width],0, 'post'); % double pad
    YCbCr = rgb2ycbcr(RGB);           
    Y = double(YCbCr(:,:,1));
    
    fft2R = fft2(double(RGB(:,:,1)));
    fft2G = fft2(double(RGB(:,:,2)));  
    fft2B = fft2(double(RGB(:,:,3)));  
    fft2Y = fft2(Y);                   
    
    if hitter == 1
        h_77 = (1/49) .* ones(7, 7);  %2 
        fft2Filt_77 = fft2(h_77, 2*height, 2*width);     % kernel is ready
        fft2FiltR_77 = fft2Filt_77 .* fft2R;             % 2-D fast convolution
        fft2FiltG_77 = fft2Filt_77 .* fft2G;
        fft2FiltB_77 = fft2Filt_77 .* fft2B;
        newRGB_77(:, :,1) = uint8(ifft2(fft2FiltR_77));     % Inverse 2-D FFT of  components
        newRGB_77(:, :,2) = uint8(ifft2(fft2FiltG_77));
        newRGB_77(:, :,3) = uint8(ifft2(fft2FiltB_77));   
        figure(2), image(newRGB_77(1:height, 1:width,:));                 
        title('7x7 smooth');
        figure(3); freqz2(h_77); title("filter profile");
    end
    
    if hitter == 2
        h_33 = (1/9) .* ones(3, 3);   %3
        fft2Filt_33 = fft2(h_33, 2*height, 2*width);     % kernel is ready
        fft2FiltR_33 = fft2Filt_33 .* fft2R;             % 2-D fast convolution
        fft2FiltG_33 = fft2Filt_33 .* fft2G;
        fft2FiltB_33 = fft2Filt_33 .* fft2B;
        newRGB_33(:, :,1) = uint8(ifft2(fft2FiltR_33));     % Inverse 2-D FFT of  components
        newRGB_33(:, :,2) = uint8(ifft2(fft2FiltG_33));
        newRGB_33(:, :,3) = uint8(ifft2(fft2FiltB_33));   
        figure(2), image(newRGB_33(1:height, 1:width,:));                 
        title('3x3 smooth');
        figure(3); freqz2(h_33); title("filter profile");
    end
    
    if hitter == 3
        h_med = fspecial('average');
        fft2Filt_med = fft2(h_med, 2*height, 2*width);
        fft2FiltR_med = fft2Filt_med .* fft2R;            
        fft2FiltG_med = fft2Filt_med .* fft2G;
        fft2FiltB_med = fft2Filt_med .* fft2B;
        newRGB_med(:, :,1) = uint8(ifft2(fft2FiltR_med));   
        newRGB_med(:, :,2) = uint8(ifft2(fft2FiltG_med));
        newRGB_med(:, :,3) = uint8(ifft2(fft2FiltB_med));   
        figure(2), image(newRGB_med(1:height, 1:width,:));                 
        title('median');
        figure(3); freqz2(h_med); title("filter profile");
    end
    
    if hitter == 4
        h_sob = fspecial('sobel');
        fft2Filt_sob = fft2(h_sob, 2*height, 2*width);
        fft2FiltR_sob = fft2Filt_sob .* fft2R;            
        fft2FiltG_sob = fft2Filt_sob .* fft2G;
        fft2FiltB_sob = fft2Filt_sob .* fft2B;
        newRGB_sob(:, :,1) = uint8(ifft2(fft2FiltR_sob));   
        newRGB_sob(:, :,2) = uint8(ifft2(fft2FiltG_sob));
        newRGB_sob(:, :,3) = uint8(ifft2(fft2FiltB_sob));   
        figure(2), image(newRGB_sob(1:height, 1:width,:));                 
        title('sobel');
        figure(3); freqz2(h_sob); title("filter profile");
    end
    
    if hitter == 5
        h_lp = lpfilter('ideal', height, width, 50);
        temp = fft2(h_lp, 2*height, 2*width);
        fft2Filt_lp = fftshift(temp);
        fft2FiltR_lp = fft2Filt_lp .* fft2R;            
        fft2FiltG_lp = fft2Filt_lp .* fft2G;
        fft2FiltB_lp = fft2Filt_lp .* fft2B;
        newRGB_lp(:, :, 1) = uint8(ifft2(fft2FiltR_lp));   
        newRGB_lp(:, :, 2) = uint8(ifft2(fft2FiltG_lp));
        newRGB_lp(:, :, 3) = uint8(ifft2(fft2FiltB_lp));   
        figure(2), image(newRGB_lp(1:height, 1:width,:));                 
        title('LP 1/4');
        figure(3); freqz2(h_lp); title("filter profile");
    end
    
    if hitter == 6
        h_rec = [0  0  0  0  0  0  0  0  0; ...  
                  0  0  0  0  0  0  0  0  0; ...
                  0  0  0  0  0  0  0  0  0; ...
                  0  0  0  .1  .1  .1  0  0  0; ...
                  0  0  0  .1  1  .1  0  0  0; ...
                  0  0  0  .1  .1  .1  0  0  0; ...
                  0  0  0  0  0  0  0  0  0; ...
                  0  0  0  0  0  0  0  0  0; ...
                  0  0  0  0  0  0  0  0  0]; 
        fft2Filt_rec = fft2(h_rec, 2*height, 2*width);
        fft2FiltR_rec = fft2Filt_rec .* fft2R;            
        fft2FiltG_rec = fft2Filt_rec .* fft2G;
        fft2FiltB_rec = fft2Filt_rec .* fft2B;
        newRGB_rec(:, :,1) = uint8(ifft2(fft2FiltR_rec));   
        newRGB_rec(:, :,2) = uint8(ifft2(fft2FiltG_rec));
        newRGB_rec(:, :,3) = uint8(ifft2(fft2FiltB_rec));   
        figure(2), image(newRGB_rec(1:height, 1:width,:));                 
        title('rectangle');
        figure(3); freqz2(h_rec); title("filter profile");
    end
    
end


%------------------------------------------------------------------------------------------
if (select == 18)
    RGB = imread('jetplane.png');            % Read the image in
    fid = fopen('dctFile.dat', 'w');
    [imHeight imwidth color] = size(RGB);       % Get image dimensions

    YCbCr = rgb2ycbcr(RGB);                     % Converts to YCbCr color space
    Y = double(YCbCr(:,:,1));                   % Change to YUV space 
    U = double(YCbCr(:,:,2))-128;
    V = double(YCbCr(:,:,3))-128;
    kk = 1;
    for n=1:8:imHeight                         % Each block is 8 pixels in height
      for m=1:8:imwidth                        % Each block is 8 pixels in width
        for i=0:7
          for j=0:7
            mbY(i+1,j+1) = Y(n+i,m+j,1);       % Form 8x8 micro-block for Y          
            mbU(i+1,j+1) = U(n+i,m+j,1);       % Form 8x8 micro-block for U     
            mbV(i+1,j+1) = V(n+i,m+j,1);       % Form 8x8 micro-block for V     
          end
        end

        mbY = dct(double(mbY));                % Perform 1-D DCT horizontally          
        mbY = dct(double(mbY'));               % Perform 1-D DCT vertically          
        mbU = dct(double(mbU));                % Perform 1-D DCT horizontally          
        mbU = dct(double(mbU'));               % Perform 1-D DCT vertically          
        mbV = dct(double(mbV));                % Perform 1-D DCT horizontally          
        mbV = dct(double(mbV'));               % Perform 1-D DCT vertically          
        oY = int16(mbY);  
        oU = int16(mbU);  
        oV = int16(mbV);      
        fwrite(fid, oY, 'int16');  
        fwrite(fid, oU, 'int16');  
        fwrite(fid, oV, 'int16');  
        
        hold_3523 = zeros(8, 8);
        hold_2637 = zeros(8, 8);
        for ii = 1:8:imHeight
            for jj = 1:8:imwidth
                kk = kk + 1;
                if kk == 3523
                    fprintf("\nblock:  %d\n", kk);
                    for xx = 0:1:7
                        for yy = 0:1:7
                            hold_3523(xx+1,yy+1) = oY(xx+1,yy+1);
                            fprintf("%9d ", hold_3523(xx+1,yy+1));
                        end
                        fprintf("\n");
                    end
                end
                if kk == 2637
                    fprintf("\nblock:  %d\n", kk);
                    for xx = 0:1:7
                        for yy = 0:1:7
                            hold_2637(xx+1,yy+1) = oY(xx+1,yy+1);
                            fprintf("%9d ", hold_2637(xx+1,yy+1));
                        end
                        fprintf("\n");
                    end
                end
            end
        end

      end
    end

    fclose(fid);                              % Close DCT data file

    %%%%%%%%%%%%%%%%%%%%%%%%%%%% iDCT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fid = fopen('dctFile.dat', 'r');          % Open DCT data file to read
    dctData = fread(fid, inf, 'int16');       % Read DCT data in
    dctData = int16(dctData);

    k = 1;
    for n=1:8:imHeight
        for m=1:8:imwidth

            for j=1:8
                for i=1:8
                    imbY(i,j) = dctData(k);
                    k = k+1;
                end
            end
            for j=1:8
                for i=1:8
                    imbU(i,j) = dctData(k);
                    k = k+1;
                end
            end
            for j=1:8
                for i=1:8
                    imbV(i,j) = dctData(k);
                    k = k+1;
                end
            end

            imbY = idct(double(imbY));            % 1-D iDCT horizontally
            imbY = idct(double(imbY'));           % 1-D iDCT vertically
            imbU = idct(double(imbU));            % 1-D iDCT horizontally
            imbU = idct(double(imbU'));           % 1-D iDCT vertically
            imbV = idct(double(imbV));            % 1-D iDCT horizontally
            imbV = idct(double(imbV'));           % 1-D iDCT vertically
            for i=0:7
                for j=0:7
                    iY(n+i,m+j) = imbY(i+1,j+1);
                    iU(n+i,m+j) = imbU(i+1,j+1);
                    iV(n+i,m+j) = imbV(i+1,j+1);
                end
            end
        end
    end

    newYCbCr(:,:,1) = uint8(round(iY));      % Form YCbCr data format  
    newYCbCr(:,:,2) = uint8(round(iU) + 128);
    newYCbCr(:,:,3) = uint8(round(iV) + 128);

    newRGB = ycbcr2rgb(newYCbCr);            % Convert YCbCr to RGB color space

    figure (1), imshow(RGB);                 % Show the original image
    title('The original image');

    figure (2), imshow(newRGB);              % Show the reconstructed image
    title('The DCT/IDCT reconstructed image');

    fclose(fid);							% Close DCT data file
    delete dctFile.dat;
    
end


%------------------------------------------------------------------------------------------
if (select == 99)
    
end
%%%%%%%%~~~~~~~~END>
