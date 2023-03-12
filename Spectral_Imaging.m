close all
clear all

% HW 12: spectral imaging
% MJMurdoch 202010, using Roy Berns' homework and images

% load all the target images into a "stack"
img(:,:,:,1) = single(imread('slot0.tif'))./(2.^16-1);
img(:,:,:,2) = single(imread('slot1.tif'))./(2.^16-1);
img(:,:,:,3) = single(imread('slot2.tif'))./(2.^16-1);
img(:,:,:,4) = single(imread('slot3.tif'))./(2.^16-1);
img(:,:,:,5) = single(imread('slot4.tif'))./(2.^16-1);
img(:,:,:,6) = single(imread('slot5.tif'))./(2.^16-1);
img(:,:,:,7) = single(imread('slot6.tif'))./(2.^16-1);

% load all the white images into a "stack"
imgw(:,:,:,1) = single(imread('white0.tif'))./(2.^16-1);
imgw(:,:,:,2) = single(imread('white1.tif'))./(2.^16-1);
imgw(:,:,:,3) = single(imread('white2.tif'))./(2.^16-1);
imgw(:,:,:,4) = single(imread('white3.tif'))./(2.^16-1);
imgw(:,:,:,5) = single(imread('white4.tif'))./(2.^16-1);
imgw(:,:,:,6) = single(imread('white5.tif'))./(2.^16-1);
imgw(:,:,:,7) = single(imread('white6.tif'))./(2.^16-1);

%% combine RGB to monochrome by adding RGB, then flatfield

im = squeeze(sum(img,3)); % sum across the 3rd dimension (R,G,B)
imw = squeeze(sum(imgw,3));
imageff = im ./ imw;

%% show images

% uncomment if you dare: will take time and memory to display
% figure
% for i = 1:7
%    subplot(4,7,i)
%    imshow(imresize(img(:,:,:,i),1/4));
%    if i==1, ylabel('Original'); end
%    subplot(4,7,i+7)
%    imshow(imresize(imgw(:,:,:,i),1/4));
%    if i==1, ylabel('White'); end
%    subplot(4,7,i+14)
%    imshow(imresize(im(:,:,i),1/4));
%    if i==1, ylabel('Raw single channel'); end
%    subplot(4,7,i+21)
%    imshow(imresize(imageff(:,:,i),1/4)); 
%    if i==1, ylabel('FF single channel'); end
% end


%% rescale

scalar = (max(max(imageff)));
scalar = scalar.*(.5);

% scale everything together to keep color balance
image = imageff ./ max(max(max(scalar)));

% change to a 2-d variable
[nr nc nb]=size(image);
pixels=reshape(image, nr*nc, nb);

%% use mask for sampling patches (find avg camera signals for each patch)

mask = imread('mask_ccsg.tif');
mask7 = repmat(mask,1,1,7);

% image of CCSG is rotated CCW 90 degrees
% patchC will be a vector N patches by 7 channels
patchCraw = zeros(max(max(mask)),7);
for i = 1:max(max(mask))
    patchCraw(i,:) = mean(reshape(image(mask7==i),[],7));
end
% reshape to 3D array and "unrotate"
patchC3 = rot90(permute(reshape(patchCraw,[10 14 7]),[2 1 3]),3);
% convert to 2D list in CCSG order
patchC = reshape(permute(patchC3,[2 1 3]),[],7);



%% load CCSG spectral data

[num,txt] = xlsread('ccsg.xlsx');
%reflectance factor
SGref = num(:,5:end)';
%wavelength
SGwl = (380:10:730)';

[num,txt] = xlsread('all_1nm_data.xlsx');
cmf2 = interp1(num(:,1),num(:,6:8),SGwl);
D65 = interp1(num(:,1),num(:,3),SGwl);
SGXYZ = ((cmf2' * diag(D65) * SGref) ./ (cmf2(:,2)' * D65) )';

%% plot CCSG spectral reflectance
figure (1);
plot(SGwl,SGref)
xlabel('wl'); ylabel('reflectance factor');
xlim([380 730]);

%% Compute CCSG D65 CIELAB 
% render sRGB colors: clipped to gamut-map
SGRGB = min(max(xyz2rgb(SGXYZ),0),1);

SGLab = xyz2lab(SGXYZ);
figure (2);
scatter3(SGLab(:,2),SGLab(:,3),SGLab(:,1),50,SGRGB,'filled');
axis equal
xlabel('a*'); ylabel('b*'); zlabel('L*');


%% find a matrix transform from 7-channel to XYZ
%Question 1: [In the m-file there is an "example" matrix to convert from 
%7-channel image to an approximation of XYZ. You should replace that matrix
%with a least-squares fit matrix (use the \ operator) that gives a better 
%estimate of XYZ. Include your matrix in your PDF, and report the min, max, and 
%mean ΔE00 between your estimate and the actual CCSG.]

% Compute XYZ matrix from CCSG camera data and known XYZ
% Use Matlab matrix-left-divide operator: \

%Calculating the matrix with a least-squares fit matrix (use the \ operator) 
%that gives a better estimate of XYZ. 
%Tristimulus values are already calculated which is in "SGXYZ"
% "patchC" contains data for 140 CCSG patches for 7-channel
M = (patchC\SGXYZ)';

% example M: this will provide a starting point, but is not accurate
% M = [-0.5   1.6  -2.4   1.8   0.5   0.5  -0.2
%      -0.4   1.0  -1.0   2.2   0.0   0.3  -0.2
%      -0.3   3.7  -3.1   1.9  -0.8   0.5  -0.2 ];

% replace M with your optimized matrix

% compute Delta E (D65 white already)
SGEstXYZ = (M*patchC')';
SGEstLab = xyz2lab(SGEstXYZ);
SGDE = deltaE00(SGEstLab',SGLab')';

%"SGDE" is already calculated the ΔE00 between the estimate and the actual CCSG.
%The min, max, and mean of "SGDE"
deltaE00_max = max(SGDE);
deltaE00_mean = mean(SGDE);
deltaE00_min = min(SGDE);

%Combining min, max, and mean of ΔE00 
deltaE00_max_mean_min = [deltaE00_max deltaE00_mean deltaE00_min];

%% Plot difference in CIELAB

figure (3);
% plot actual CCSG colors
scatter3(SGLab(:,2),SGLab(:,3),SGLab(:,1),50,SGRGB,'filled');
axis equal
xlabel('a*'); ylabel('b*'); zlabel('L*')

% plot estimated CCSG colors
hold on
scatter3(SGEstLab(:,2),SGEstLab(:,3),SGEstLab(:,1),50,SGRGB);
legend('Actual','Estimate')


%% convert the whole image to XYZ and then sRGB

% go from Cam Signals to XYZ to sRGB
pixXYZ = (pixels)*M';
pixRGB = xyz2rgb(pixXYZ);
imgRGB = reshape(pixRGB,[nr nc 3]);

%% Question 2: 
% Generate an sRGB image from the XYZ pixel data and save it as a JPG.
figure (4);
imshow(imgRGB)
imwrite(imgRGB,'ImageOutput_sRGB.jpg');


%% spectral estimation
%Question 3:
%Create a matrix to estimate spectra for each of the 140 CCSG patches. 
%Make a plot of your matrix, like Fig. 6.42 in the book. 
%And, make a plot of the actual and estimated spectra for the 6 patches near 
%the middle of the CCSG: blue, green, red, yellow, magenta, cyan.

%Equation 6.27
%This Matrix (M) contains the wavelengths are in a rows and
%the camera 7-channel are in a columns.
Matrix = SGref*pinv(patchC');

%Initializing for each channel to distinguish the plot
C1 = Matrix(:,1);
C2 = Matrix(:,2);
C3 = Matrix(:,3);
C4 = Matrix(:,4);
C5 = Matrix(:,5);
C6 = Matrix(:,6);
C7 = Matrix(:,7);

%Plotting the matrix (Matrix), like Fig. 6.42 in the book.
figure (5);
plot(SGwl,C1,'magenta','LineWidth', 1.5);
hold on
plot(SGwl,C2,'blue', 'LineWidth', 1.5);
hold on
plot(SGwl,C3,'cyan', 'LineWidth', 1.5);
hold on
plot(SGwl,C4,'green','LineWidth', 1.5);
hold on
plot(SGwl,C5,'yellow', 'LineWidth', 1.5);
hold on
plot(SGwl,C6,'red', 'LineWidth', 1.5);
hold on
plot(SGwl,C7,'black','LineWidth', 1.5);
xlim([380 730]);
xlabel("Wavelength (nm)");
ylabel("Matrix coefficient");
title("Camera Signals to Spectral Reflectance","FontSize",18);
set(gca,'Color',[.8 .8 .8]);
legend('Channel 1', 'Channel 2', 'Channel 3', 'Channel 4',...
    'Channel 5', 'Channel 6', 'Channel 7');

%Estimating spectra for the 6 patches near the middle of the CCSG: 
%blue, green, red, yellow, magenta, cyan.
%actual patches
actual_patches = SGref(:,47:52);

%Equation 6.26
%Estimating spectra for all patches of the CCSG
R_estimate = (Matrix *patchC');

%Extracting only the patches we need which are the estimated spectra for 
%the 6 patches near the middle of the CCSG: 
%blue, green, red, yellow, magenta, cyan.
estimated_patches = R_estimate(:,47:52);

%Initializing for each actual patch to distinguish the plot
blue_patch_A = actual_patches(:,1);
green_patch_A = actual_patches(:,2);
red_patch_A = actual_patches(:,3);
yellow_patch_A = actual_patches(:,4);
magenta_patch_A = actual_patches(:,5);
cyan_patch_A = actual_patches(:,6);

%Initializing for each estimated patch to distinguish the plot
blue_patch_E = estimated_patches(:,1);
green_patch_E = estimated_patches(:,2);
red_patch_E = estimated_patches(:,3);
yellow_patch_E = estimated_patches(:,4);
magenta_patch_E = estimated_patches(:,5);
cyan_patch_E = estimated_patches(:,6);

%Plotting of the actual and estimated spectra for the 6 patches near 
%the middle of the CCSG: blue, green, red, yellow, magenta, cyan.
figure (6);
plot(SGwl,blue_patch_A, "blue", 'LineWidth', 1);
hold on
plot(SGwl,green_patch_A , "green",'LineWidth', 1);
hold on
plot(SGwl,red_patch_A, "red", 'LineWidth', 1);
hold on
plot(SGwl,yellow_patch_A, "yellow",'LineWidth', 1);
hold on
plot(SGwl,magenta_patch_A, "magenta", 'LineWidth', 1);
hold on
plot(SGwl,cyan_patch_A, "cyan",'LineWidth', 1);
hold on
plot(SGwl,blue_patch_E, '--bs', 'LineWidth',1, 'MarkerSize',4,...
    'MarkerEdgeColor','blue', 'MarkerFaceColor','blue');
hold on
plot(SGwl,green_patch_E, '--gs', 'LineWidth',1, 'MarkerSize',4,...
    'MarkerEdgeColor','green', 'MarkerFaceColor','green');
hold on
plot(SGwl,red_patch_E, '--rs', 'LineWidth',1,'MarkerSize',4,...
    'MarkerEdgeColor','red', 'MarkerFaceColor','red');
hold on
plot(SGwl,yellow_patch_E, '--ys', 'LineWidth',1, 'MarkerSize',4,...
    'MarkerEdgeColor','yellow', 'MarkerFaceColor','yellow');
hold on
plot(SGwl,magenta_patch_E,'--ms', 'LineWidth',1, 'MarkerSize',4,...
    'MarkerEdgeColor','magenta', 'MarkerFaceColor','magenta');
hold on
plot(SGwl,cyan_patch_E,'--cs', 'LineWidth',1, 'MarkerSize',4,...
    'MarkerEdgeColor','cyan', 'MarkerFaceColor','cyan');
xlim([380 730]);
xlabel("Wavelength (nm)");
ylabel("Reflectance factor");
title("Reflectance Factor Measurements","FontSize",18);
set(gca,'Color',[.8 .8 .8]);
legend("Actual: blue patch", "Actual: green patch", "Actual: red patch",...
    "Actual: yellow patch", "Actual: magenta patch", "Actual: cyan patch",...
    "Estimated: blue patch","Estimated: green patch","Estimated: red patch",...
    "Estimated: yellow patch", "Estimated: magenta patch", ...
    "Estimated: cyan patch",'Location','northwest');
