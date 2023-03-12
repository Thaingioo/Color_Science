%Thain-Gi Oo
%Homework 5: Display RGB/CMF 
 
close all
clear all
 
%Set the variable of the three_conesLMS equal to LMS.
%LMS data has the information for the L, M, and S and they represent the three types of cones with their high sensitivities in the long, medium, and short wavelength regions.
three_conesLMS = xlsread('Opponency_Data.xlsx', 'LMS');
 
%Removing first column and grabbing the rest from the LMS cones data 
three_cones=three_conesLMS(:,2:4);
 
%Set the variable of the three_primariesRGB equal to DisplaySPD.
three_primariesRGB = xlsread('DisplaySPD_Data.xlsx', 'DisplaySPD');
 
%Removing first column and grabbing the rest from the DisplaySPD data 
three_primaries=three_primariesRGB(:,2:4);
 
%wavelength
wavelength = three_conesLMS(:,1);
%Color-Matching function
%Question #1: Compute the LMS values of the display's three primaries.
%We can compute by using Equation4.4
LMS_RGB = three_cones'*three_primaries;

%------------------------------
%Question #2:Compute the color matching functions (CMFs, also known as r-bar, g-bar, b-bar) of the display’s primaries, and plot them versus wavelength. Label your plot clearly. See Eq. 4.5.
%Equation4.5: the color-matching functions = inverse of LMS_RGB * lms sensitivities transpose
%RGB Color-Matching Functions for the CIE 1931 standard observer
%It is not required to average r_bar, g_bar, and b_bar.
%Since R has one value for each wavelength, 
%G has one value for each wavelength, and
%B has one value for each wavelength.  

RGB = inv(LMS_RGB) * three_cones';

%grabbing each r, g, and b values from each row.
r = RGB(1,:);
g = RGB(2,:);
b = RGB(3,:);
%plotting a figure before normalizing 
figure(1)
plot(wavelength,r, 'red');
hold on
plot(wavelength,g, 'green');
plot(wavelength,b, 'blue');
hold on
yline( 0 );
xlim([380 730]);
xlabel('Wavelength(nm)');
ylabel('Tristimulus values');
legend('r-bar','g-bar','b-bar');
title('RGB Color-Matching Functions');
hold off 

%The RGB equation will result in the maximum tristimulus values of
%50, but we want the maximum tristimulus values to be 1, so we can rescale by using
%the MATLAB code RGB./(max(max(RGB))).
RGB_rescale=RGB./(max(max(RGB)));

%defining each of the color matching functions so we can color each of the
%functions as we want.
r_bar = RGB_rescale(1,:);
g_bar = RGB_rescale(2,:);
b_bar = RGB_rescale(3,:);



%plot for RGB Color-Matching Functions for the CIE 1931 standard observer
figure(2)
plot(wavelength,r_bar, 'red');
hold on
plot(wavelength,g_bar, 'green');
plot(wavelength,b_bar, 'blue');
hold on
yline( 0 );
xlim([380 730]);
xlabel('Wavelength(nm)');
ylabel('Tristimulus values');
legend('r-bar','g-bar','b-bar');
title('RGB Color-Matching Functions');
hold off 
%------------------------------
r1 = three_primariesRGB(:,2);
g1 = three_primariesRGB(:,3);
b1 = three_primariesRGB(:,4);
figure(3)
plot(wavelength,r1, 'red');
hold on
plot(wavelength,g1, 'green');
plot(wavelength,b1, 'blue');
xlim([380 730]);
xlabel('Wavelength(nm)');
ylabel('Radiance(W/m^2Sr)');
legend('r','g','b');
title('Spectral radiance of each primary at unit amount');
hold off 


%------------------------------
%Question #4:4. Create a 3x3 matrix to linearly transform the CMFs you computed to approximate the CIE 1931 standard colorimetric observer. 
%Reading data from CIE 1931 of DisplaySPD_Data 
CIE_value = xlsread('DisplaySPD_Data.xlsx', 'CIE 1931');
%Set the variable of the CIE_values equal to CIE 1931.
%CIE_1931 will have a data of x_bar, y_bar, and z_bar for
%each wavelength (380-730).
CIE_1931 = CIE_value(:,2:4);

%Equation 4.10: Normalizing r-bar, g-bar, and b-bar. 
r_Normalize= normalize(RGB_rescale(1,:), 'norm', 1);
g_Normalize= normalize(RGB_rescale(2,:), 'norm', 1);
b_Normalize= normalize(RGB_rescale(3,:), 'norm', 1);
%sum(abs(R_Normalize(:)))
 
% Making a matrix by combining each normalized 1X36 matrix such as r_Normalize, g_Normalize, and b_Normalize.
%Set r_Normalize in first row, g_Normalize in second row, and b_Normalize in third row.
rgb_normalized = [r_Normalize;g_Normalize;b_Normalize];
 
%rgb_normalized(matrix) x M (matrix) = (x_bar, y_bar, z_bar)  to find M (matrix) we will
%have to use \ operation. 
%Equation 4.13. M (matrix)= rgb_normalized(matrix)\(x_bar, y_bar, z_bar).
Matrix = rgb_normalized'\CIE_1931;
 
 
%Equation (4.14): Calculating CIE XYZ 1931 standard colorimetric observer
%We are approximating each x_bar, y_bar, and z_bar 
%from CIE XYZ 1931 standard colorimetric observer.
% xyz_approximatation=rgb_normalized'*Matrix*normalized 
xyz_approximatation=rgb_normalized'*Matrix;
 
%defining each function of x_bar, y_bar, and z_bar from xyz_approximation.
x_bar = xyz_approximatation(:,1);
y_bar = xyz_approximatation(:,2);
z_bar = xyz_approximatation(:,3);

 
%plotting the approximated x_bar, y_bar, and z_bar values that used the 
%color-matching function.
figure(4)
plot(wavelength,x_bar,'red.');
hold on
plot(wavelength,y_bar,'green.');
plot(wavelength,z_bar,'blue.');
xlim([380 730]);
ylim([0 1.8]);
xlabel('Wavelength(nm)');
ylabel('Tristimulus values');
legend('x-bar: LCD','y-bar: LCD','z-bar: LCD');
title('Calculated Color Matching Functions for LCD Display');

%Question #5
%setting a variable of each x_bar, y_bar, and z_bar value from CIE 1931 data
X_bar = CIE_1931(:,1);
Y_bar = CIE_1931(:,2);
Z_bar = CIE_1931(:,3);
%plotting CIE XYZ 1931 standard colorimetric observer and the approximated
% x_bar, y_bar, and z_bar values that used the color-matching function.
figure(5)
plot(wavelength,x_bar,'r.');
hold on
plot(wavelength,y_bar,'g.');
hold on
plot(wavelength,z_bar,'b.');
hold on
plot(wavelength,X_bar,'r');
hold on
plot(wavelength,Y_bar,'g');
hold on
plot(wavelength,Z_bar,'b');
xlim([380 780]);
ylim([0 1.8]);
xlabel('Wavelength(nm)');
ylabel('Tristimulus values');
legend('x-bar: LCD','y-bar: LCD','z-bar: LCD','x-bar: CIE 1931','y-bar: CIE 1931','z-bar: CIE 1931');
title('Standard Observer vs. Calculated Color Matching Functions for LCD Display');

