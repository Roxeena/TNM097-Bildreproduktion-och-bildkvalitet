%% Read in the images and convert to XYZ and LAB
clear 
close all
clc

original_RGB = im2double(imread('OnePlus6 test2jpg.tif'));
OnePlus6_RGB = im2double(imread('600D test2jpg.tif'));

%Convert to XYZ and LAB
original_XYZ = rgb2xyz(original_RGB);
original_LAB = rgb2lab(original_RGB);

OnePlus6_XYZ = rgb2xyz(OnePlus6_RGB);
OnePlus6_LAB = rgb2lab(OnePlus6_RGB);

figure;
imshow(OnePlus6_RGB);
figure;
imshow(original_RGB);

%% Calculate SNR, Requires original
OnePlus6SNR = snr(original_RGB, original_RGB - OnePlus6_RGB)    % 17.6207

%% Calculate difference, Requires original
LAB_diff_1Plus6 = zeros(size(OnePlus6_RGB(:,:,1)));
for i = 1:size(OnePlus6_RGB(:, 1, 1))
    for j = 1:size(OnePlus6_RGB(1, :, 1))
        LAB_diff_1Plus6(i,j) = sqrt( (original_LAB(i,j, 1) - OnePlus6_LAB(i,j, 1))^2 + (original_LAB(i,j, 2) - OnePlus6_LAB(i,j, 2))^2 + (original_LAB(i,j, 3) - OnePlus6_LAB(i,j, 3))^2 );
    end
end

% Calculate mean and max values in difference
LAB_mean_1Plus6 = mean(mean(LAB_diff_1Plus6)) % 0.0056
LAB_max_1Plus6 = max(max(LAB_diff_1Plus6))    % 28.0498

%% Calculate SNR with filter on, Requires origianl
% DOES NOT WORK AND DONT KNOW WHY!?
OnePlus6FilterSNR = snr_filter(original_RGB, original_RGB - OnePlus6_RGB)   %

%% Apply HVS to RGB before calculating color differense, Requires original
distance = 500; % mm
dpi = 300;
inch = 25.4; 
dotSize = 1/dpi * inch; % mm

f = MFTsp(15,dotSize,distance);
% Denna funktion returnerar ett lågpass filter som representerar ögat. I
% detta fall betraktningsavståndet har satts till 500 mm och punkternas
% storlek till 0.0847. Observera att punktens storlek motsvarar ett tryck
% i 300 dpi. (0.0847 = 1/300 * 25.4 mm) 15?

sR = conv2(original_RGB(:,:, 1),f,'same');
sG = conv2(original_RGB(:,:, 2),f,'same');
sB = conv2(original_RGB(:,:, 3),f,'same');
% Ögats filter är applicerat till signalen (originalbilden)

nR = conv2(original_RGB(:,:, 1) - OnePlus6_RGB(:,:,1),f,'same');
nG = conv2(original_RGB(:,:, 2) - OnePlus6_RGB(:,:,2),f,'same');
nB = conv2(original_RGB(:,:, 3) - OnePlus6_RGB(:,:,3),f,'same');
% Ögats filter är applicerat till "noise"-en (skillnaden mellan originalbilden och rasterbilden)

% make sure no zeros
sR=(sR>0).*sR;
sG=(sG>0).*sG;
sB=(sB>0).*sB;

nR=(nR>0).*nR;
nG=(nG>0).*nG;
nB=(nB>0).*nB;

% Combine the channels to one image
original_RGB_S = zeros(size(original_RGB));
original_RGB_S(:,:, 1) = sR;
original_RGB_S(:,:, 2) = sG;
original_RGB_S(:,:, 3) = sB;

OnePlus6_RGB_N = zeros(size(OnePlus6_RGB));
OnePlus6_RGB_N(:,:, 1) = nR;
OnePlus6_RGB_N(:,:, 2) = nG;
OnePlus6_RGB_N(:,:, 3) = nB;

% Calculate color differense
original_LAB_S = rgb2lab(original_RGB_S);
OnePlus6_LAB_N = rgb2lab(OnePlus6_RGB_N);

LAB_diff_1Plus6_SN = zeros(size(OnePlus6_LAB_N(:,:,1)));
for i = 1:size(OnePlus6_LAB_N(:, 1, 1))
    for j = 1:size(OnePlus6_LAB_N(1, :, 1))
        LAB_diff_1Plus6_SN(i,j) = sqrt( (original_LAB_S(i,j, 1) - OnePlus6_LAB_N(i,j, 1))^2 + (original_LAB_S(i,j, 2) - OnePlus6_LAB_N(i,j, 2))^2 + (original_LAB_S(i,j, 3) - OnePlus6_LAB_N(i,j, 3))^2 );
    end
end

% Calculate mean and max values in difference
LAB_mean_OnePlus6_N = mean(mean(LAB_diff_1Plus6_SN))    % 0.0345
LAB_max_OnePlus6_N = max(max(LAB_diff_1Plus6_SN))       % 83.6210

%% Calculate SCIE Lab differense, Requires original
% Create lightsource
CIED65 = [95.05, 100, 108.9];
CIEA = [109.85, 100.00, 35.58];

ppi = 120;
d = 5; % inch
sampPerDeg = ppi * d * tan(pi/180);

OnePlus6_SCIELab = scielab(sampPerDeg, original_XYZ, OnePlus6_XYZ, CIEA, 'xyz');

SCIE_mean_nearest = mean(mean(OnePlus6_SCIELab))   % 1.9047
SCIE_max_nearest = max(max(OnePlus6_SCIELab))      % 24.1767

%% Calculate graininess, Requires SCIELAB
OnePlus6_SCIELab_1 = scielab(sampPerDeg, OnePlus6_XYZ);

OnePlus6_L1 = OnePlus6_SCIELab_1(:,:, 1);
OnePlus6_a1 = OnePlus6_SCIELab_1(:,:, 2);
OnePlus6_b1 = OnePlus6_SCIELab_1(:,:, 3);

graniness_OnePlus6 = sum(std(OnePlus6_L1(:)) + std(OnePlus6_a1(:)) + std(OnePlus6_b1(:)))    % 0.1562

%% Calculate SSIM, Requires original

[OnePlus6_SSIMVal, OnePlus6_SSIMMap] = ssim(OnePlus6_RGB, original_RGB);    % 0.9561

figure;
imshow(OnePlus6_SSIMMap);
