
close all
clear all

%Importing data
%TwoDegChromaticity contains (x,y)
%lab1 contains selected 10 color patches
%lab2 contains selected 1 color patches with 10 measurements
%Illuminant Data contains (A, D50, D55, D65, D75)
%StdObsFuncs contains CIE (1931) 2-deg and 10-deg color matching functions
deg2_Chro = xlsread('TwoDegChromaticity.xlsx');
ten_patches = xlsread('lab1.xlsx');
green_patch1 = xlsread('lab2.xlsx');
illu = xlsread('Illuminant Data.xlsx');
stdObs = xlsread('StdObsFuncs.xlsx');


%Removing the the name (wavelengths) and choosing data between 380-780
green_patch = (green_patch1(:,3:42))';
deg2_Chroma = deg2_Chro(:,2:3);
illumin = illu(17:97,2:6);
deg2_CMF = stdObs(5:85,2:4);
wl=deg2_Chro(:,1);
colors = ten_patches(2:11,3:42);
G_patch  = green_patch1(2:11,3:42);
wl1 = ten_patches(1,3:42);

%Interpolating and Extrapolating of green patch 
GP = interp1(green_patch(:,1),green_patch(:,2:end),380:5:780,'linear','extrap');

%% Question 1
%plotting 10 patches of their reflectance factors
figure(1)
plot(wl1,colors);
legend('Red-27', 'Cyan-18', 'Magenta-29', 'Brown-40', 'Green-14',...
    'Light Blue-32', 'Black:42', 'Yellow-28', 'White-41', 'Light Skin-2');
xlabel('Wavelength (nm)');
ylabel('Reflectance factor');
title('Measurement of 10 Samples');
xlim([360 900]);
%% Question 2
%plotting a green patch with 10 measurements of their reflectance factors
figure(2)
plot(wl1,G_patch);
legend('trial1','trial2','trial3','trial4','trial5','trial6',...
    'trial7','trial8','trial9','trial10');
xlabel('Wavelength (nm)');
ylabel('Reflectance factor');
title('10 Measurements of 1 Sample');

%% Question 3
%Create a variable using the data from 2-deg color matching functions (CMF) 
x_bar2=deg2_CMF(:,1);
y_bar2=deg2_CMF(:,2);
z_bar2=deg2_CMF(:,3);

%Create a variable using the data from the CIE standard D65 illuminant.
D65=illumin(:,4);

%Since D65 are 81x1 vectors, we can transform to 81x81 by using
%diagonalization method.
Dia_D65=diag(D65);

%calculating delta
delta=mean(diff(wl));

%2-DEG standard observer under illuminant D65.
%Calculating k constant for 2-deg standard observers under illuminant D65. 
%Equation 4.21 (page 63).
k_2_D65=100./(D65.'*y_bar2.*delta);

%2-DEG standard observer and D65 illuminant
%Calculating CIE XYZ tristimulus values for 2-deg standard observer under
%D65 illuminant.
%Equation 4.18 to 4.21 (page 62-63).
X_deg2_D65 = k_2_D65.*((Dia_D65*GP)'*x_bar2).*delta;
Y_deg2_D65 = k_2_D65.*((Dia_D65*GP)'*y_bar2).*delta;
Z_deg2_D65 = k_2_D65.*((Dia_D65*GP)'*z_bar2).*delta;

%Creating a matrix for calculated CIE XYZ tristimulus values for 2-deg 
%standard observer under D65 illuminant.
XYZ_deg2_D65 = [X_deg2_D65, Y_deg2_D65, Z_deg2_D65];

%Calculating xy chromaticity coordinates for 2-deg standard observer 
%under D65 illuminant.
%Equation 4.22 to 4.24 (page 65).
x_deg2_D65 = X_deg2_D65./(X_deg2_D65+Y_deg2_D65+Z_deg2_D65);
y_deg2_D65 = Y_deg2_D65./(X_deg2_D65+Y_deg2_D65+Z_deg2_D65);

%Creating a matrix for the xy chromaticity coordinates for 
%2-deg standard observer under D65 illuminant.
xy_deg2_D65 = [x_deg2_D65 y_deg2_D65];


%Creating a line between the spectrum locus's curve initial and 
%end points. 
x = [0.174 0.7368];
y = [0.0048 0.2632];
%---------------------------------------------------------------------
%Plotting the 2-deg (x,y) results for Illum D65 on a 
%chromaticity diagram
figure(3)
labels_a = {'1','2','3','4','5','6','7','8','9','10'};
labels_b = {'Spectrum Locus'};
labels_c = {'480 (nm)'};
labels_d = {'520 (nm)'};
labels_e = {'560 (nm)'};
labels_f = {'780 (nm)'};

scatter(x_deg2_D65,y_deg2_D65,45,'MarkerEdgeColor','blue',...
              'MarkerFaceColor','cyan',...
              'LineWidth',1);
hold on 

%Spectrum locus
plot(deg2_Chro(:,2),deg2_Chro(:,3),'black');
hold on 
plot(x,y,'magenta');
labels_P = {'Purple Line'};
text(x(:,1),y(:,1),labels_P,'VerticalAlignment','bottom','HorizontalAlignment','left');
text(x_deg2_D65,y_deg2_D65,labels_a,'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',7.5);
text(deg2_Chro(35,2),deg2_Chro(35,3),labels_b,'VerticalAlignment','bottom','HorizontalAlignment','left')
text(deg2_Chro(21,2),deg2_Chro(21,3),labels_c,'VerticalAlignment','top','HorizontalAlignment','right',FontSize=6);
text(deg2_Chro(29,2),deg2_Chro(29,3),labels_d,'VerticalAlignment','top','HorizontalAlignment','left',FontSize=6);
text(deg2_Chro(37,2),deg2_Chro(37,3),labels_e,'VerticalAlignment','top','HorizontalAlignment','left',FontSize=6);
text(deg2_Chro(81,2),deg2_Chro(81,3),labels_f,'VerticalAlignment','top','HorizontalAlignment','left',FontSize=6);
title('Chromaticity Diagram (2-DEG)','FontSize',14);
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
axis equal
xlim([0 1]);
ylim([0 1]);

%% CIELAB
%Calculating CIELAB values (L*, a*, b*, and C*). 
%Using the CIE 1931 under D65.

%Calculating the white points, using the 2-deg XYZ values for 
%D65 illuminant. 
%Equation 4.18 to 4.21, without using reflectance factor. (page 62-63)
X_n_D65=k_2_D65.*(D65'*x_bar2).*delta;
Y_n_D65=k_2_D65.*(D65'*y_bar2).*delta;
Z_n_D65=k_2_D65.*(D65'*z_bar2).*delta;

%Creating matrix for white points (D65).
XYZ_n_D65 = [X_n_D65,Y_n_D65,Z_n_D65];

%Calculating X', Y', and Z' for D65 illuminant.
%Equation 4.64 to 4.66. (page 75)
X_prime_D65 = X_deg2_D65/X_n_D65;
Y_prime_D65 = Y_deg2_D65/Y_n_D65;
Z_prime_D65 = Z_deg2_D65/Z_n_D65;

%Creating matrix for X', Y', and Z' (D65). 
XYZ_prime_D65 = [X_prime_D65, Y_prime_D65, Z_prime_D65];

%Calculating the function from equation (4.73). (page 75)
minVal = (24/116)^3;

xyz_D65 = XYZ_prime_D65;
if (xyz_D65  > minVal)
    y_D65=xyz_D65.^(1/3);
elseif (xyz_D65  <= (minVal))
    y_D65=(841/108)*xyz_D65 +16/116;
end

%Calculating the L*, a*, b*, and C* from equations (4.70-4.72 and 4.80).
%(pages 75 and 76)
%D65 illuminant
L_star_D65 = (116.*y_D65(:,2)-16);
a_star_D65 = (500.*(y_D65(:,1)-y_D65(:,2)));
b_star_D65 = (200.*(y_D65(:,2)-y_D65(:,3)));
C_star_D65 = sqrt(a_star_D65.^2+b_star_D65.^2);
LabC_star_D65 = [L_star_D65, a_star_D65, b_star_D65, C_star_D65];

%Plotting the a*, b* values of each patch on the a*b* planes in
%CIELAB space for D65 illuminant. Note, it used 2-deg XYZ values.
figure(4)
labels3 = {'1','2','3','4','5','6','7','8','9','10'};
scatter(a_star_D65,b_star_D65,45,'MarkerEdgeColor','blue',...
              'MarkerFaceColor','cyan',...
              'LineWidth',1);
hold on
xline(0);
hold on
yline(0);
title('CIELAB','FontSize',14);
xlabel('a^*',fontsize=14);
ylabel('b^*',fontsize=14);
text(a_star_D65,b_star_D65,labels3,'VerticalAlignment','bottom','HorizontalAlignment','left');
axis equal
xlim([-200 200]);
ylim([-50 300]);


%Plotting the a*, b* values of each patch on the a*b* planes in
%CIELAB space for D65 illuminant. Note, it used 2-deg XYZ values.
%x-limit ([-120.2 -119.7]) and y-limit ([255 255.5])
figure(5)
labels3 = {'1','2','3','4','5','6','7','8','9','10'};
scatter(a_star_D65,b_star_D65,45,'MarkerEdgeColor','blue',...
              'MarkerFaceColor','cyan',...
              'LineWidth',1);
title('CIELAB','FontSize',14);
xlabel('a^*',fontsize=14);
ylabel('b^*',fontsize=14);
text(a_star_D65,b_star_D65,labels3,'VerticalAlignment','bottom','HorizontalAlignment','left');
axis equal
xlim([-120.2 -119.7]);
ylim([255 255.5]);


%Plotting the L*, C* values of each patch on the L*C* planes 
%in CIELAB space for D65 illuminant.
figure (6)
labels4 = {'1','2','3','4','5','6','7','8','9','10'};
scatter(C_star_D65,L_star_D65,45,'MarkerEdgeColor','black',...
              'MarkerFaceColor','yellow',...
              'LineWidth',1);
title('CIELAB','FontSize',14);
xlabel('C^*_{ab}',fontsize=14);
ylabel('L^*',fontsize=14);
text(C_star_D65,L_star_D65,labels4,'VerticalAlignment','bottom','HorizontalAlignment','left');
axis equal
xlim([0 400]);
ylim([0 400]);

%Plotting the L*, C* values of each patch on the L*C* planes 
%in CIELAB space for D65 illuminant.
figure (7)
labels4 = {'1','2','3','4','5','6','7','8','9','10'};
scatter(C_star_D65,L_star_D65,45,'MarkerEdgeColor','black',...
              'MarkerFaceColor','yellow',...
              'LineWidth',1);
title('CIELAB','FontSize',14);
xlabel('C^*_{ab}',fontsize=14);
ylabel('L^*',fontsize=14);
text(C_star_D65,L_star_D65,labels4,'VerticalAlignment','bottom','HorizontalAlignment','left');
axis equal

%Calculating the mean (average) of L*, a*, and b*
L_bar_star= mean(L_star_D65);
a_bar_star= mean(a_star_D65);
b_bar_star= mean(b_star_D65);

%Equation 6.5 (page 125)
n=10;
MCDM_value = value(L_star_D65,a_star_D65,b_star_D65,...
    L_bar_star,a_bar_star,b_bar_star,n);
function MCDM = value(L,a,b,L_bar,a_bar,b_bar,n)
    MCDM = (sum(((L-L_bar).^2+(a-a_bar).^2+(b-b_bar).^2).^(1/2)))./n;
end










