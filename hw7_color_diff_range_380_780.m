
close all
clear all

%Importing data
%TwoDegChromaticity data contains (x,y)
%Colorchecker data contains 24 color patches
%Illuminant data contains (A, D50, D55, D65, D75)
%StdObsFuncs data contains CIE (1931) 2-deg and 10-deg color matching functions
%InkjetColorChecker data contains 12 color patches
deg2_Chro = xlsread('TwoDegChromaticity.xlsx');
colorChecker = xlsread('MacbethColorChecker.xlsx');
illu = xlsread('Illuminant Data.xlsx');
stdObs = xlsread('StdObsFuncs.xlsx');
inkjetcolorChecker = xlsread('InkjetColorChecker.xlsx');

%Removing the name (wavelengths) and choosing data between 380-780
deg2_Chroma = deg2_Chro(:,2:3);
%Color checker from 13 to 24.
macbeth = colorChecker(2:82,14:end);

%Interpolating and Extrapolating inkjetcolorChecker
inkjetcolor1 = interp1(inkjetcolorChecker(:,1),inkjetcolorChecker(:,2:end),380:5:780,'linear','extrap');

%Normalize inkjetColorCheckers
inkjetcolor = inkjetcolor1/100;
ink = normalize(inkjetcolor1,'range');
ink12 = [inkjetcolor ink];

illumin = illu(17:97,2:6);
deg2_CMF = stdObs(5:85,2:4);
wl=deg2_Chro(1:81,1);


%Create a variable using the data from 2-deg color matching functions (CMF) 
xyz_bar2=[deg2_CMF(:,1) deg2_CMF(:,2) deg2_CMF(:,3)];

%Create a variable using the data from the CIE standard D50 illuminant, 
%A illuminant, and D65 illuminant.
A=illumin(:,1);
D65=illumin(:,4);

%Since D50, A, and D65 are 81x1 vectors, we can transform to 81x81 by using
%diagonalization method.
Dia_A=diag(A);
Dia_D65=diag(D65);

%calculating delta lambda
delta=mean(diff(wl));

%2-DEG standard observer under A illuminant and D65 illuminant.
%Calculating k constant for 2-deg standard observers under A illuminant 
%and D65 illuminant. 
%Equation 4.21 (page 63).
k_2_A=100./(A.'*xyz_bar2(:,2).*delta);
k_2_D65=100./(D65.'*xyz_bar2(:,2).*delta);


%% Macbeth ColorChecker 
%2-DEG standard observer and A illuminant
%Calculating CIE XYZ tristimulus values for 2-deg standard observer under
%A illuminant.
%Equations 4.18 to 4.21 (page 62-63).
XYZ_deg2_A_Macbeth = k_2_A.*((Dia_A*macbeth)'*xyz_bar2).*delta;

%Calculating xy chromaticity coordinates for 2-deg standard observer 
%under A illuminant.
%Equations 4.22 to 4.24 (page 65).
xy_deg2_A_Macbeth = [XYZ_deg2_A_Macbeth(:,1)./(XYZ_deg2_A_Macbeth(:,1)+XYZ_deg2_A_Macbeth(:,2)+XYZ_deg2_A_Macbeth(:,3)) XYZ_deg2_A_Macbeth(:,2)./(XYZ_deg2_A_Macbeth(:,1)+XYZ_deg2_A_Macbeth(:,2)+XYZ_deg2_A_Macbeth(:,3))];

%2-DEG standard observer and D65 illuminant
%Calculating CIE XYZ tristimulus values for 2-deg standard observer under
%D65 illuminant.

%Macbeth ColorChecker (D65 illuminant)
%Equations 4.18 to 4.21 (page 62-63).
XYZ_deg2_D65_Macbeth = k_2_D65.*((Dia_D65*macbeth)'*xyz_bar2).*delta;

%Calculating xy chromaticity coordinates for 2-deg standard observer 
%under D65 illuminant.
%Equations 4.22 to 4.24 (page 65).
xy_deg2_D65_Macbeth = [XYZ_deg2_D65_Macbeth(:,1)./(XYZ_deg2_D65_Macbeth(:,1)+XYZ_deg2_D65_Macbeth(:,2)+XYZ_deg2_D65_Macbeth(:,3)) XYZ_deg2_D65_Macbeth(:,2)./(XYZ_deg2_D65_Macbeth(:,1)+XYZ_deg2_D65_Macbeth(:,2)+XYZ_deg2_D65_Macbeth(:,3))];


%Calculating the white points using the 2-deg XYZ values for A illuminant.
%Equations 4.18 to 4.21, without using reflectance factor. (page 62-63)
XYZ_n_A_Macbeth=(k_2_A.*(A'*xyz_bar2).*delta);
 
%Calculating X', Y', and Z' for A illuminant.
%Equations 4.64 to 4.66. (page 75)
XYZ_prime_A_Macbeth = [XYZ_deg2_A_Macbeth(:,1)/XYZ_n_A_Macbeth(:,1) XYZ_deg2_A_Macbeth(:,2)/XYZ_n_A_Macbeth(:,2) XYZ_deg2_A_Macbeth(:,3)/XYZ_n_A_Macbeth(:,3)];

%Calculating the white points, using the 2-deg XYZ values for 
%D65 illuminant. 
%Equations 4.18 to 4.21, without using reflectance factor. (page 62-63)
XYZ_n_D65_Macbeth=k_2_D65.*(D65'*xyz_bar2).*delta;

%Calculating X', Y', and Z' for D65 illuminant.
%Equations 4.64 to 4.66. (page 75)
XYZ_prime_D65_Macbeth = [XYZ_deg2_D65_Macbeth(:,1)/XYZ_n_D65_Macbeth(:,1) XYZ_deg2_D65_Macbeth(:,2)/XYZ_n_D65_Macbeth(:,2) XYZ_deg2_D65_Macbeth(:,3)/XYZ_n_D65_Macbeth(:,3)];

%Calculating the function from equation (4.73). (page 75)
minVal = (24/116)^3;

xyz_A = XYZ_prime_A_Macbeth;
if (xyz_A > minVal)
    y_A=xyz_A.^(1/3);
elseif (xyz_A <= (minVal))
    y_A=(841/108)*xyz_A+16/116;
end

xyz_D65 = XYZ_prime_D65_Macbeth;
if (xyz_D65  > minVal)
    y_D65=xyz_D65.^(1/3);
elseif (xyz_D65  <= (minVal))
    y_D65=(841/108)*xyz_D65 +16/116;
end


%A illuminant
%Calculating the L*, a*, b*, and C* from equations (4.70-4.72 and 4.80).
%(pages 75 and 76)
L_Std_A = (116.*y_A(:,2)-16);
a_Std_A = (500.*(y_A(:,1)-y_A(:,2)));
b_Std_A = (200.*(y_A(:,2)-y_A(:,3)));
C_ab_Std_A = sqrt(a_Std_A.^2+b_Std_A.^2);
CIELAB_A_macbeth = [L_Std_A a_Std_A b_Std_A C_ab_Std_A]';

%D65 illuminant
%Calculating the L*, a*, b*, and C* from equations (4.70-4.72 and 4.80).
%(pages 75 and 76)
L_Std_D65 = (116.*y_D65(:,2)-16);
a_Std_D65 = (500.*(y_D65(:,1)-y_D65(:,2)));
b_Std_D65 = (200.*(y_D65(:,2)-y_D65(:,3)));
C_ab_Std_D65 = sqrt(a_Std_D65.^2+b_Std_D65.^2);
CIELAB_D65_macbeth = [L_Std_D65 a_Std_D65 b_Std_D65 C_ab_Std_D65]';

%% Inkjet ColorChecker 
% Question 1: Compute and make tables of CIELAB values of the inkjet-printed patches from the spectral reflectance data for Illuminants A and D65. 
%Equations 4.18 to 4.21 (page 62-63).
XYZ_deg2_A_inkjet = k_2_A.*((Dia_A*inkjetcolor)'*xyz_bar2).*delta;
%Calculating xy chromaticity coordinates for 2-deg standard observer 
%under A illuminant.
%Equations 4.22 to 4.24 (page 65).
xy_deg2_A_inkjet= [XYZ_deg2_A_inkjet(:,1)./(XYZ_deg2_A_inkjet(:,1)+XYZ_deg2_A_inkjet(:,2)+XYZ_deg2_A_inkjet(:,3)) XYZ_deg2_A_inkjet(:,2)./(XYZ_deg2_A_inkjet(:,1)+XYZ_deg2_A_inkjet(:,2)+XYZ_deg2_A_inkjet(:,3))];
%2-DEG standard observer and D65 illuminant
%Calculating CIE XYZ tristimulus values for 2-deg standard observer under
%D65 illuminant.

%Equations 4.18 to 4.21 (page 62-63).
XYZ_deg2_D65_inkjet = k_2_D65.*((Dia_D65*inkjetcolor)'*xyz_bar2).*delta;

%Calculating xy chromaticity coordinates for 2-deg standard observer 
%under D65 illuminant.
%Equations 4.22 to 4.24 (page 65).
xy_deg2_D65_inkjet = [XYZ_deg2_D65_inkjet(:,1)./(XYZ_deg2_D65_inkjet(:,1)+XYZ_deg2_D65_inkjet(:,2)+XYZ_deg2_D65_inkjet(:,3)) XYZ_deg2_D65_inkjet(:,2)./(XYZ_deg2_D65_inkjet(:,1)+XYZ_deg2_D65_inkjet(:,2)+XYZ_deg2_D65_inkjet(:,3))];


%Calculating the white points using the 2-deg XYZ values for A illuminant.
%Equations 4.18 to 4.21, without using reflectance factor. (page 62-63)
XYZ_n_A_inkjet=(k_2_A.*(A'*xyz_bar2).*delta);
 
%Calculating X', Y', and Z' for A illuminant.
%Equations 4.64 to 4.66. (page 75)
XYZ_prime_A_inkjet = [XYZ_deg2_A_inkjet(:,1)/XYZ_n_A_inkjet(:,1) XYZ_deg2_A_inkjet(:,2)/XYZ_n_A_inkjet(:,2) XYZ_deg2_A_inkjet(:,3)/XYZ_n_A_inkjet(:,3)];

%Calculating the white points using the 2-deg XYZ values for 
%D65 illuminant. 
%Equations 4.18 to 4.21, without using reflectance factor. (page 62-63)
XYZ_n_D65_inkjet=k_2_D65.*(D65'*xyz_bar2).*delta;

%Calculating X', Y', and Z' for D65 illuminant.
%Equations 4.64 to 4.66. (page 75)
XYZ_prime_D65_inkjet = [XYZ_deg2_D65_inkjet(:,1)/XYZ_n_D65_inkjet(:,1) XYZ_deg2_D65_inkjet(:,2)/XYZ_n_D65_inkjet(:,2) XYZ_deg2_D65_inkjet(:,3)/XYZ_n_D65_inkjet(:,3)];

%Calculating the function from equation (4.73). (page 75)
minVal = (24/116)^3;

xyz_A_inkjet = XYZ_prime_A_inkjet;
if (xyz_A_inkjet > minVal)
    y_A_inkjet=xyz_A_inkjet.^(1/3);
elseif (xyz_A_inkjet <= (minVal))
    y_A_inkjet=(841/108)*xyz_A_inkjet+16/116;
end

xyz_D65_inkjet = XYZ_prime_D65_inkjet;
if (xyz_D65_inkjet  > minVal)
    y_D65_inkjet=xyz_D65_inkjet.^(1/3);
elseif (xyz_D65_inkjet  <= (minVal))
    y_D65_inkjet=(841/108)*xyz_D65_inkjet +16/116;
end


%A illuminant
%Calculating the L*, a*, b*, and C* from equations (4.70-4.72 and 4.80).
L_Bat_A = (116.*y_A_inkjet(:,2)-16);
a_Bat_A = (500.*(y_A_inkjet(:,1)-y_A_inkjet(:,2)));
b_Bat_A = (200.*(y_A_inkjet(:,2)-y_A_inkjet(:,3)));
C_ab_Bat_A = sqrt(a_Bat_A.^2+b_Bat_A.^2);

CIELAB_A_inkjet = [L_Bat_A a_Bat_A b_Bat_A C_ab_Bat_A]';

%D65 illuminant
%Calculating the L*, a*, b*, and C* from equations (4.70-4.72 and 4.80).
L_Bat_D65 = (116.*y_D65_inkjet(:,2)-16);
a_Bat_D65 = (500.*(y_D65_inkjet(:,1)-y_D65_inkjet(:,2)));
b_Bat_D65 = (200.*(y_D65_inkjet(:,2)-y_D65_inkjet(:,3)));
C_ab_Bat_D65 = sqrt(a_Bat_D65.^2+b_Bat_D65.^2);

CIELAB_D65_inkjet = [L_Bat_D65 a_Bat_D65 b_Bat_D65 C_ab_Bat_D65]';

%% Question 2
%Calculate and make tables of ΔL*, Δa*, Δb*, ΔC*, ΔH*, ΔE*ab and ΔE*00 color
%differences between the 12 corresponding patches on the inkjet-printed ColorChecker 
% (sample set) and the real ColorChecker (reference set). 

%Illuminant A
%Calculating the differences between two coordinates of lightness for batch
%and lightness for standard by using
%Equation 5.1
delta_L_star_A = L_Bat_A - L_Std_A;

%Calculating the differences between two coordinates of a^*_bat and a^*_Std
%by using
%Equation 5.2
delta_a_star_A = a_Bat_A - a_Std_A;

%Calculating the differences between two coordinates of b^*_bat and b^*_Std
%by using
%Equation 5.3
delta_b_star_A = b_Bat_A - b_Std_A;

%Calculating the differences between two coordinates of C^*_ab,bat and 
%C^*_ab,Std by using
%Equation 5.5
delta_C_ab_star_A = C_ab_Bat_A - C_ab_Std_A;

%Using the Seve formula to calculate delta_H_ab_star_A.
%Equation 5.8
delta_H_ab_star_A = ((a_Bat_A.*b_Std_A)- (a_Std_A.*b_Bat_A))./((0.5.*((C_ab_Bat_A.*C_ab_Std_A)+(a_Bat_A.*a_Std_A)+(b_Bat_A.*b_Std_A))).^(1/2));

%"The CIELAB Euclidean distance between two colors is defined as delta
%E_ab."
%Equation 5.9
delta_E_ab_star_A = (delta_L_star_A.^2+delta_a_star_A.^2+delta_b_star_A.^2).^(1/2);


%______________________________________________
%Illuminant D65

%Calculating the differences between two coordinates of lightness for batch
%and lightness for standard by using
%Equation 5.1
delta_L_star_D65 = L_Bat_D65 - L_Std_D65;

%Calculating the differences between two coordinates of a^*_bat and a^*_Std
%by using
%Equation 5.2
delta_a_star_D65 = a_Bat_D65 - a_Std_D65;

%Calculating the differences between two coordinates of b^*_bat and b^*_Std
%by using
%Equation 5.3
delta_b_star_D65 = b_Bat_D65 - b_Std_D65;


%Calculating the differences between two coordinates of C^*_ab,bat and 
%C^*_ab,Std by using
%Equation 5.3
delta_C_ab_star_D65 = C_ab_Bat_D65 - C_ab_Std_D65;

%Using the Seve formula to calculate delta_H_ab_star_A.
%Equation 5.8
delta_H_ab_star_D65 = (a_Bat_D65.*b_Std_D65- a_Std_D65.*b_Bat_D65)./((0.5.*(C_ab_Bat_D65.*C_ab_Std_D65+a_Bat_D65.*a_Std_D65+b_Bat_D65.*b_Std_D65)).^(1/2));


%"The CIELAB Euclidean distance between two colors is defined as delta
%E_ab."
%Equation 5.9
delta_E_ab_star_D65 = (delta_L_star_D65.^2+delta_a_star_D65.^2+delta_b_star_D65.^2).^(1/2);

%_____________________________________________________
%Calculating Delta E^*_00 for illuminant A and D65 using the given funtion.
delta_E00_A = deltaE00(CIELAB_A_inkjet,CIELAB_A_macbeth);
delta_E00_D65 = deltaE00(CIELAB_D65_inkjet,CIELAB_D65_macbeth);

deltas_A = [delta_L_star_A delta_a_star_A delta_b_star_A ...
    delta_C_ab_star_A delta_H_ab_star_A delta_E_ab_star_A delta_E00_A'];

deltas_D65 = [delta_L_star_D65 delta_a_star_D65  delta_b_star_D65 ...
    delta_C_ab_star_D65 delta_H_ab_star_D65 delta_E_ab_star_D65 delta_E00_D65'];

%% plot standard (a*-b* plane)
figure;
labels1 = {'13','14','15','16','17','18','19','20','21','22','23','24'};
scatter(a_Std_A ,b_Std_A,45,'MarkerEdgeColor','blue',...
              'MarkerFaceColor','cyan',...
              'LineWidth',1);
hold on
xline(0);
hold on
yline(0);
title('CIELAB (Standard)','FontSize',14);
xlabel('a^*',fontsize=14);
ylabel('b^*',fontsize=14);
text(a_Std_A,b_Std_A,labels1,'VerticalAlignment','bottom','HorizontalAlignment','left')
legend('A');
axis equal
xlim([-90 90]);
ylim([-90 90]);

figure;
labels3 = {'13','14','15','16','17','18','19','20','21','22','23','24'};
scatter(a_Std_D65,b_Std_D65,45,'MarkerEdgeColor','blue',...
              'MarkerFaceColor','cyan',...
              'LineWidth',1);
hold on
xline(0);
hold on
yline(0);
title('CIELAB (Standard)','FontSize',14);
xlabel('a^*',fontsize=14);
ylabel('b^*',fontsize=14);
text(a_Std_D65,b_Std_D65,labels3,'VerticalAlignment','bottom','HorizontalAlignment','left');
legend('D65');
axis equal
xlim([-90 90]);
ylim([-90 90]);

%% plot batch (a*-b* plane)
figure;
labels1 = {'25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36'};
scatter(a_Bat_A ,b_Bat_A,45,'MarkerEdgeColor','blue',...
              'MarkerFaceColor','cyan',...
              'LineWidth',1);
hold on
xline(0);
hold on
yline(0);
title('CIELAB (Batch)','FontSize',14);
xlabel('a^*',fontsize=14);
ylabel('b^*',fontsize=14);
text(a_Bat_A,b_Bat_A,labels1,'VerticalAlignment','bottom','HorizontalAlignment','left')
legend('A');
axis equal
xlim([-90 90]);
ylim([-90 90]);

figure;
labels3 = {'25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36'};

scatter(a_Bat_D65,b_Bat_D65,45,'MarkerEdgeColor','blue',...
              'MarkerFaceColor','cyan',...
              'LineWidth',1);
hold on
xline(0);
hold on
yline(0);
title('CIELAB (Batch)','FontSize',14);
xlabel('a^*',fontsize=14);
ylabel('b^*',fontsize=14);
text(a_Bat_D65,b_Bat_D65,labels3,'VerticalAlignment','bottom','HorizontalAlignment','left');
legend('D65');
axis equal
xlim([-90 90]);
ylim([-90 90]);
