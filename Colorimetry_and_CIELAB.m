
close all
clear all

%Importing data
%TwoDegChromaticity contains (x,y)
%Colorchecker contains 24 color patches
%Illuminant Data contains (A, D50, D55, D65, D75)
%StdObsFuncs contains CIE (1931) 2-deg and 10-deg color matching functions
deg2_Chro = xlsread('TwoDegChromaticity.xlsx');
colorChecker = xlsread('MacbethColorChecker.xlsx');
illu = xlsread('Illuminant Data.xlsx');
stdObs = xlsread('StdObsFuncs.xlsx');

%Removing the the name (wavelengths) and choosing data between 380-780
deg2_Chroma = deg2_Chro(:,2:3);
macbeth = colorChecker(2:82,2:end);
illumin = illu(17:97,2:6);
deg2_CMF = stdObs(5:85,2:4);
deg10_CMF = stdObs(5:85,7:9);
wl=deg2_Chro(:,1);

%Create a variable using the data from 2-deg color matching functions (CMF) 
x_bar2=deg2_CMF(:,1);
y_bar2=deg2_CMF(:,2);
z_bar2=deg2_CMF(:,3);

%Create a variable using the data from 10-deg color matching functions (CMF) 
x_bar10=deg10_CMF(:,1);
y_bar10=deg10_CMF(:,2);
z_bar10=deg10_CMF(:,3);

%Create a variable using the data from the CIE standard D50 illuminant, 
%A illuminant, and D65 illuminant.
D50 = illumin(:,2);
A=illumin(:,1);
D65=illumin(:,4);

%Since D50, A, and D65 are 81x1 vectors, we can transform to 81x81 by using
%diagonalization method.
Dia_D50=diag(D50);
Dia_A=diag(A);
Dia_D65=diag(D65);

%calculating delta
delta=mean(diff(wl));

%% Question #1
%Calculate CIE XYZ tristimulus values and chromaticity coordinates for 
%each ColorChecker patch using the CIE standard D50 
%illuminant and each of the 2-deg and 10-deg standard observers. 


%Calculating k constant for 2-deg and 10-deg standard observers under D50 
%illuminant. 
%Equation 4.21 (page 63).
k_2_D50=100./(D50.'*y_bar2.*delta);
k_10_D50=100./(D50.'*y_bar10.*delta);
%---------------------------------------------------------------------

%2-DEG standard observer and D50 illuminant
%Calculating CIE XYZ tristimulus values for 2-deg standard observer under 
%D50 illuminant.
%Equation 4.18 to 4.21 (page 62-63).
X_deg2_D50 = k_2_D50.*((Dia_D50*macbeth)'*x_bar2).*delta;
Y_deg2_D50 = k_2_D50.*((Dia_D50*macbeth)'*y_bar2).*delta;
Z_deg2_D50 = k_2_D50.*((Dia_D50*macbeth)'*z_bar2).*delta;

%Creating a matrix for calculated CIE XYZ tristimulus values for 2-deg 
%standard observer under D50 illuminant.
XYZ_deg2_D50 = [X_deg2_D50, Y_deg2_D50, Z_deg2_D50];

%Calculating xy chromaticity coordinates for 2-deg standard observer 
%under D50 illuminant.
%Equation 4.22 to 4.24 (page 65).
x_deg2_D50 = X_deg2_D50./(X_deg2_D50+Y_deg2_D50+Z_deg2_D50);
y_deg2_D50 = Y_deg2_D50./(X_deg2_D50+Y_deg2_D50+Z_deg2_D50);

%Creating a matrix for the xy chromaticity coordinates for 
%2-deg standard observer under D50 illuminant.
xy_deg2_D50 = [x_deg2_D50 y_deg2_D50];
%---------------------------------------------------------------------

%10-DEG standard observer and D50 illuminant
%Calculating CIE XYZ tristimulus values for 10-deg standard observer under
%D50 illuminant.
%Equation 4.18 to 4.21 (page 62-63).
X_deg10_D50 = k_10_D50.*((Dia_D50*macbeth)'*x_bar10).*delta;
Y_deg10_D50 = k_10_D50.*((Dia_D50*macbeth)'*y_bar10).*delta;
Z_deg10_D50 = k_10_D50.*((Dia_D50*macbeth)'*z_bar10).*delta;

%Creating a matrix for calculated CIE XYZ tristimulus values for 10-deg 
%standard observer under D50 illuminant.
XYZ_deg10_D50 = [X_deg10_D50, Y_deg10_D50, Z_deg10_D50];

%Calculating xy chromaticity coordinates for 10-deg standard observer
%under D50 illuminant.
%Equation 4.22 to 4.24 (page 65).
x_deg10_D50 = X_deg10_D50./(X_deg10_D50+Y_deg10_D50+Z_deg10_D50);
y_deg10_D50 = Y_deg10_D50./(X_deg10_D50+Y_deg10_D50+Z_deg10_D50);

%Creating a matrix for the xy chromaticity coordinates for 
%10-deg standard observer under D50 illuminant.
xy_deg10_D50 = [x_deg10_D50 y_deg10_D50];

%% Question 2
%Calculate CIE XYZ tristimulus values and chromaticity coordinates for each
%ColorChecker patch using the CIE 2-deg standard observer and the A and D65
%illuminants. Please put these in a table, one row per patch.

%2-DEG standard observer under A illuminant and D65 illuminant.
%Calculating k constant for 2-deg standard observers under A illuminant 
%and D65 illuminant. 
%Equation 4.21 (page 63).
k_2_A=100./(A.'*y_bar2.*delta);
k_2_D65=100./(D65.'*y_bar2.*delta);
%---------------------------------------------------------------------

%2-DEG standard observer and A illuminant
%Calculating CIE XYZ tristimulus values for 2-deg standard observer under
%A illuminant.
%Equation 4.18 to 4.21 (page 62-63).
X_deg2_A = k_2_A.*((Dia_A*macbeth)'*x_bar2).*delta;
Y_deg2_A = k_2_A.*((Dia_A*macbeth)'*y_bar2).*delta;
Z_deg2_A = k_2_A.*((Dia_A*macbeth)'*z_bar2).*delta;

%Creating a matrix for calculated CIE XYZ tristimulus values for 2-deg 
%standard observer under A illuminant.
XYZ_deg2_A = [X_deg2_A, Y_deg2_A, Z_deg2_A];

%Calculating xy chromaticity coordinates for 2-deg standard observer 
%under A illuminant.
%Equation 4.22 to 4.24 (page 65).
x_deg2_A = X_deg2_A./(X_deg2_A+Y_deg2_A+Z_deg2_A);
y_deg2_A = Y_deg2_A./(X_deg2_A+Y_deg2_A+Z_deg2_A);

%Creating a matrix for the xy chromaticity coordinates for 
%2-deg standard observer under A illuminant.
xy_deg2_A = [x_deg2_A y_deg2_A];
%---------------------------------------------------------------------

%2-DEG standard observer and D65 illuminant
%Calculating CIE XYZ tristimulus values for 2-deg standard observer under
%D65 illuminant.
%Equation 4.18 to 4.21 (page 62-63).
X_deg2_D65 = k_2_D65.*((Dia_D65*macbeth)'*x_bar2).*delta;
Y_deg2_D65 = k_2_D65.*((Dia_D65*macbeth)'*y_bar2).*delta;
Z_deg2_D65 = k_2_D65.*((Dia_D65*macbeth)'*z_bar2).*delta;

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

%% Question 3
%[Plot all of your 2-deg (x,y) results (for Illum A, D50, and D65) on a 
%chromaticity diagram, with a legend to make it clear which points go with
%which illuminant. You can label the points (1-24) using the text function 
%in Matlab, and you can plot the spectrum locus using the (x,y) coordinates
%provided. In chromaticity diagrams, it is important to scale the x-axis 
%and y-axis equally (think about why). In Matlab please use the command 
%axis equal.]

%Creating a line between the spectrum locus's curve initial and 
%end points. 
x = [0.174 0.7368];
y = [0.0048 0.2632];
%---------------------------------------------------------------------
% %Plotting the 2-deg (x,y) results (for Illum A, D50, and D65) on a 
% %chromaticity diagram
figure(1)
labels_a = {'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24'};
labels_b = {'Spectrum Locus'};
%labels_c = {'380','385','390','395','400','405','410','415','420','425','430','435','440','445','450','455','460','465','470','475','480','485','490','495','500','505','510','515','520','525','530','535','540','545','550','555','560','565','570','575','580','585','590','595','600','605','610','615','620','625','630','635','640','645','650','655','660','665','670','675','680','685','690','695','600','705','710','715','720','725','730','735','740','745','750','755','760','765','770','775','780'};
labels_c = {'480 (nm)'};
labels_d = {'520 (nm)'};
labels_e = {'560 (nm)'};
labels_f = {'780 (nm)'};
scatter(x_deg2_D50,y_deg2_D50,45,'MarkerEdgeColor','black',...
              'MarkerFaceColor','magenta',...
              'LineWidth',1);
hold on
scatter(x_deg2_A,y_deg2_A,45,'MarkerEdgeColor','black',...
              'MarkerFaceColor','yellow',...
              'LineWidth',1);
hold on
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
text(x_deg2_D50,y_deg2_D50,labels_a,'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',7.5);
text(x_deg2_A,y_deg2_A,labels_a,'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',7.5);
text(x_deg2_D65,y_deg2_D65,labels_a,'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',7.5);
text(deg2_Chro(35,2),deg2_Chro(35,3),labels_b,'VerticalAlignment','bottom','HorizontalAlignment','left')
text(deg2_Chro(21,2),deg2_Chro(21,3),labels_c,'VerticalAlignment','top','HorizontalAlignment','right',FontSize=6);
text(deg2_Chro(29,2),deg2_Chro(29,3),labels_d,'VerticalAlignment','top','HorizontalAlignment','left',FontSize=6);
text(deg2_Chro(37,2),deg2_Chro(37,3),labels_e,'VerticalAlignment','top','HorizontalAlignment','left',FontSize=6);
text(deg2_Chro(81,2),deg2_Chro(81,3),labels_f,'VerticalAlignment','top','HorizontalAlignment','left',FontSize=6);
title('Chromaticity Diagram (2-DEG)','FontSize',14);
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
legend('D50','A','D65');
axis equal
xlim([0 1]);
ylim([0 1]);

%% Question 4
%What do you observe about the relative chromaticity coordinates for 
%the different light sources?

%Plotting A, D50, D55, D65, and D75 illuminants
figure(2)
plot(wl,illumin);
xlabel('Wavelength (nm)');
ylabel('Illuminant');
legend('A','D50','D55','D65','D75');



%% Question 5
%Why did I not ask you to plot the 10-deg (x,y) values with the rest?

%Calculating k constant for 10-deg standard observers under A illuminant 
%and D65 illuminant.
%Equation 4.21 (page 63).
k_10_A=100./(A.'*y_bar10.*delta);
k_10_D65=100./(D65.'*y_bar10.*delta);

%10-DEG standard observer and A illuminant
%Calculating CIE XYZ tristimulus values for 10-deg standard observer under
%A illuminant.
%Equation 4.18 to 4.21 (page 62-63).
X_deg10_A = k_10_A.*((Dia_A*macbeth)'*x_bar10).*delta;
Y_deg10_A = k_10_A.*((Dia_A*macbeth)'*y_bar10).*delta;
Z_deg10_A = k_10_A.*((Dia_A*macbeth)'*z_bar10).*delta;

%Calculating xy chromaticity coordinates for 10-deg standard observer
%under A illuminant.
%Equation 4.22 to 4.24 (page 65).
x_deg10_A = X_deg10_A./(X_deg10_A+Y_deg10_A+Z_deg10_A);
y_deg10_A = Y_deg10_A./(X_deg10_A+Y_deg10_A+Z_deg10_A);

%10-DEG standard observer and D65 illuminant.
%Calculating CIE XYZ tristimulus values for 2-deg standard observer under
%D65 illuminant.
%Equation 4.18 to 4.21 (page 62-63).
X_deg10_D65 = k_10_D65.*((Dia_D65*macbeth)'*x_bar10).*delta;
Y_deg10_D65 = k_10_D65.*((Dia_D65*macbeth)'*y_bar10).*delta;
Z_deg10_D65 = k_10_D65.*((Dia_D65*macbeth)'*z_bar10).*delta;

%Calculating xy chromaticity coordinates for 10-deg standard observer 
%under D65 illuminant.
%Equation 4.22 to 4.24 (page 65).
x_deg10_D65 = X_deg10_D65./(X_deg10_D65+Y_deg10_D65+Z_deg10_D65);
y_deg10_D65 = Y_deg10_D65./(X_deg10_D65+Y_deg10_D65+Z_deg10_D65);


%plotting degree 2 standard observer and degree 10 standard observer
%simultaneously under 3 different light conditions such as A, D50, and D65 
%illuminant.
figure(3)
labels_g = {'Spectrum Locus'};
labels_h = {'480 (nm)'};
labels_i = {'520 (nm)'};
labels_j = {'560 (nm)'};
labels_k = {'780 (nm)'};
scatter(x_deg2_D50,y_deg2_D50,45,'MarkerEdgeColor','black',...
              'MarkerFaceColor','red',...
              'LineWidth',1);
hold on
scatter(x_deg2_A,y_deg2_A,45,'MarkerEdgeColor','black',...
              'MarkerFaceColor','red',...
              'LineWidth',1);
hold on
scatter(x_deg2_D65,y_deg2_D65,45,'MarkerEdgeColor','black',...
              'MarkerFaceColor','red',...
              'LineWidth',1);
hold on
scatter(x_deg10_D50,y_deg10_D50,45,'MarkerEdgeColor','black',...
              'MarkerFaceColor','cyan',...
              'LineWidth',1);
hold on
scatter(x_deg10_A,y_deg10_A,45,'MarkerEdgeColor','black',...
              'MarkerFaceColor','cyan',...
              'LineWidth',1);
hold on
scatter(x_deg10_D65,y_deg10_D65,45,'MarkerEdgeColor','black',...
              'MarkerFaceColor','cyan',...
              'LineWidth',1);
plot(deg2_Chro(:,2),deg2_Chro(:,3),'black');
hold on 
plot(x,y,'magenta');
title('Chromaticity Diagram (2-DEG & 10-DEG)','FontSize',14);
labels_L = {'Purple Line'};
text(x(:,1),y(:,1),labels_L,'VerticalAlignment','bottom','HorizontalAlignment','left');
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
text(deg2_Chro(35,2),deg2_Chro(35,3),labels_g,'VerticalAlignment','bottom','HorizontalAlignment','left');
text(deg2_Chro(21,2),deg2_Chro(21,3),labels_h,'VerticalAlignment','top','HorizontalAlignment','right',FontSize=6);
text(deg2_Chro(29,2),deg2_Chro(29,3),labels_i,'VerticalAlignment','top','HorizontalAlignment','left',FontSize=6);
text(deg2_Chro(37,2),deg2_Chro(37,3),labels_j,'VerticalAlignment','top','HorizontalAlignment','left',FontSize=6);
text(deg2_Chro(81,2),deg2_Chro(81,3),labels_k,'VerticalAlignment','top','HorizontalAlignment','left',FontSize=6);
legend('2D-D50','2D-A','2D-D65','10D-D50','10D-A','10D-D65');
axis equal
xlim([0 1]);
ylim([0 1]);


%% Question 7
%[Using the 2-deg XYZ values for both D65 and A illuminants, calculate
%CIELAB values (L*, a*, b*, and C*). Please put these in a table, one row
%per patch. Note that you will need to determine what you will use as your
%white XnYnZn for each calculation, and please list the values in your 
%answer.]


%Calculating the white points, using the 2-deg XYZ values for A illuminant.
%Equation 4.18 to 4.21, without using reflectance factor. (page 62-63)
X_n_A=(k_2_A.*(A'*x_bar2).*delta);
Y_n_A=(k_2_A.*(A'*y_bar2).*delta);
Z_n_A=k_2_A.*(A'*z_bar2).*delta;

%Creating matrix for white points (A).
XYZ_n_A = [X_n_A,Y_n_A,Z_n_A];

%Calculating X', Y', and Z' for A illuminant.
%Equation 4.64 to 4.66. (page 75)
X_prime_A = X_deg2_A/X_n_A;
Y_prime_A = Y_deg2_A/Y_n_A;
Z_prime_A = Z_deg2_A/Z_n_A;

%Creating matrix for X', Y', and Z' (D65). 
XYZ_prime_A = [X_prime_A, Y_prime_A, Z_prime_A];

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

xyz_A = XYZ_prime_A;
if (xyz_A > minVal)
    y_A=xyz_A.^(1/3);
elseif (xyz_A <= (minVal))
    y_A=(841/108)*xyz_A+16/116;
end

xyz_D65 = XYZ_prime_D65;
if (xyz_D65  > minVal)
    y_D65=xyz_D65.^(1/3);
elseif (xyz_D65  <= (minVal))
    y_D65=(841/108)*xyz_D65 +16/116;
end

%Calculating the L*, a*, b*, and C* from equations (4.70-4.72 and 4.80).
%(pages 75 and 76)

%A illuminant
L_star_A = (116.*y_A(:,2)-16);
a_star_A = (500.*(y_A(:,1)-y_A(:,2)));
b_star_A = (200.*(y_A(:,2)-y_A(:,3)));
C_star_A = sqrt(a_star_A.^2+b_star_A.^2);
LabC_star_A = [L_star_A, a_star_A, b_star_A, C_star_A];

%D65 illuminant
L_star_D65 = (116.*y_D65(:,2)-16);
a_star_D65 = (500.*(y_D65(:,1)-y_D65(:,2)));
b_star_D65 = (200.*(y_D65(:,2)-y_D65(:,3)));
C_star_D65 = sqrt(a_star_D65.^2+b_star_D65.^2);
LabC_star_D65 = [L_star_D65, a_star_D65, b_star_D65, C_star_D65];

%% Question 8
%[Plot the L*, C* and a*, b* values of each patch on the a*b* and L*C* 
%planes in CIELAB space. Make two sets of plots: one for D65 and one for A. 
%Please label the points with patch number and/or color them appropriately.
%In these plots of a uniform color space, it is important to scale the 
%x-axis and y-axis equally (think about why). In Matlab please use the 
%command axis equal.]

%Plotting the a*, b* values of each patch on the a*b* planes in
%CIELAB space for A illuminant. Note, it used 2-deg XYZ values.
figure(4)
labels1 = {'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24'};
scatter(a_star_A,b_star_A,45,'MarkerEdgeColor','blue',...
              'MarkerFaceColor','cyan',...
              'LineWidth',1);
hold on
xline(0);
hold on
yline(0);
title('CIELAB','FontSize',14);
xlabel('a^*',fontsize=14);
ylabel('b^*',fontsize=14);
text(a_star_A,b_star_A,labels1,'VerticalAlignment','bottom','HorizontalAlignment','left')
legend('A');
axis equal
xlim([-90 90]);
ylim([-90 90]);

%Plotting the L*, C* values of each patch on the L*C* 
%planes in CIELAB space for A illuminant.
figure (5)
labels2 = {'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24'};
scatter(C_star_A,L_star_A,45,'MarkerEdgeColor','black',...
              'MarkerFaceColor','yellow',...
              'LineWidth',1);
title('CIELAB','FontSize',14);
xlabel('C^*_{ab}',fontsize=14);
ylabel('L^*',fontsize=14);
legend('A');
text(C_star_A,L_star_A,labels2,'VerticalAlignment','bottom','HorizontalAlignment','left')
axis equal
xlim([0 100]);
ylim([0 100]);

%Plotting the a*, b* values of each patch on the a*b* planes in
%CIELAB space for D65 illuminant. Note, it used 2-deg XYZ values.
figure(6)
labels3 = {'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24'};
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
legend('D65');
axis equal
xlim([-90 90]);
ylim([-90 90]);

%Plotting the L*, C* values of each patch on the L*C* planes 
%in CIELAB space for D65 illuminant.
figure (7)
labels4 = {'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24'};
scatter(C_star_D65,L_star_D65,45,'MarkerEdgeColor','black',...
              'MarkerFaceColor','yellow',...
              'LineWidth',1);
title('CIELAB','FontSize',14);
xlabel('C^*_{ab}',fontsize=14);
ylabel('L^*',fontsize=14);
text(C_star_D65,L_star_D65,labels4,'VerticalAlignment','bottom','HorizontalAlignment','left')
legend('D65');
axis equal
xlim([0 100]);
ylim([0 100]);

