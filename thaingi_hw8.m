
%Thain-Gi Oo
%Homework 8

%Import data
DuPont = xlsread('RIT_DuPont_Data_2019.xlsx');

%Extract the delta V values of 243 colors from DuPont Data
delta_V = DuPont(:,4);

%Extract the L*, a*, and b* values of 243 colors from DuPont Data
%Calculating Standard C*
L_star_Std = DuPont(:,5);
a_star_Std = DuPont(:,6);
b_star_Std = DuPont(:,7);
C_ab_star_Std = (a_star_Std.^2+b_star_Std.^2).^(1/2);
CIELAB_Std = [L_star_Std a_star_Std b_star_Std C_ab_star_Std]';

%Extract the L*, a*, and b* values of 243 colors from DuPont Data
%Calculating TR C*
L_star_TR = DuPont(:,8);
a_star_TR = DuPont(:,9);
b_star_TR = DuPont(:,10);
C_ab_star_TR = (a_star_TR.^2+b_star_TR.^2).^(1/2);
CIELAB_TR = [L_star_TR a_star_TR b_star_TR C_ab_star_TR]';

%Calculating the differences between two coordinates of lightness for TR
%and lightness for standard by using
%Equation 5.1
delta_L_star = L_star_TR - L_star_Std;

%Calculating the differences between two coordinates of a^*_TR and a^*_Std
%by using
%Equation 5.2
delta_a_star = a_star_TR - a_star_Std;

%Calculating the differences between two coordinates of b^*_TR and b^*_Std
%by using
%Equation 5.3
delta_b_star = b_star_TR - b_star_Std;


%Calculating the differences between two coordinates of C^*_ab,TR and 
%C^*_ab,Std by using
%Equation 5.3
delta_C_ab_star = C_ab_star_TR - C_ab_star_Std;

%Using the Seve formula to calculate delta_H_ab_star_A.
%Equation 5.8
delta_H_ab_star = (a_star_TR.*b_star_Std- a_star_Std.*b_star_TR)./...
    ((0.5.*(C_ab_star_TR.*C_ab_star_Std+a_star_TR.*a_star_Std+b_star_TR.*b_star_Std))...
    .^(1/2));

%% Question 3

%Calculating color differences between two colors by using the method of Euclidean distance
%Equation 5.9
delta_E_ab = (delta_L_star.^2+delta_a_star.^2+delta_b_star.^2).^(1/2);

%Calculate the color differences by using CIEDE2000
%Equation 5.35
delta_E_00 = deltaE00(CIELAB_TR,CIELAB_Std)';

%Calculate the color differences by using CIE94
%Equation 5.30
delta_E_94 = deltaE94(delta_L_star,C_ab_star_Std,C_ab_star_TR,delta_C_ab_star,delta_H_ab_star);

%Calculating maximum, mininmum, and mean of delta Eab.
delta_E_ab_max = max(delta_E_ab);
delta_E_ab_min = min(delta_E_ab);
delta_E_ab_avg = mean(delta_E_ab);

delta_E_ab_mma = [delta_E_ab_max delta_E_ab_min delta_E_ab_avg];

%Calculating maximum, mininmum, and mean of delta E00.
delta_E_00_max = max(delta_E_00);
delta_E_00_min = min(delta_E_00);
delta_E_00_avg = mean(delta_E_00);

delta_E_00_mma = [delta_E_00_max delta_E_00_min delta_E_00_avg];

%Calculating maximum, mininmum, and mean of delta E94.
delta_E_94_max = max(delta_E_94);
delta_E_94_min = min(delta_E_94);
delta_E_94_avg = mean(delta_E_94);

delta_E_94_mma = [delta_E_94_max delta_E_94_min delta_E_94_avg];

%% Question 4

%Calculating STRESS for Eab, E00, and E94.
stress_delta_E_ab = stress(delta_E_ab,delta_V);
stress_delta_E_00 = stress(delta_E_00,delta_V);
stress_delta_E_94 = stress(delta_E_94,delta_V);
stress3 = [stress_delta_E_ab stress_delta_E_94 stress_delta_E_00];

%% –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

%Extract the delta V values of 156 colors from DuPont Data
delta_V_156 = DuPont(1:156,4);

dV = delta_V_156;

%Extract the L*, a*, and b* values of 156 colors from DuPont Data
%Calculating Standard C*
L_star_Std_156 = DuPont(1:156,5);
a_star_Std_156 = DuPont(1:156,6);
b_star_Std_156 = DuPont(1:156,7);
C_ab_star_Std_156 = (a_star_Std_156.^2+b_star_Std_156.^2).^(1/2);
CIELAB_Std_156 = [L_star_Std_156 a_star_Std_156 b_star_Std_156 C_ab_star_Std_156]';

%Extract the L*, a*, and b* values of 156 colors from DuPont Data
%Calculating TR C*
L_star_TR_156 = DuPont(1:156,8);
a_star_TR_156 = DuPont(1:156,9);
b_star_TR_156 = DuPont(1:156,10);
C_ab_star_TR_156 = (a_star_TR_156.^2+b_star_TR_156.^2).^(1/2);
CIELAB_TR_156 = [L_star_TR_156 a_star_TR_156 b_star_TR_156 C_ab_star_TR_156]';

%Calculating the differences between two coordinates of lightness for TR
%and lightness for standard by using
%Equation 5.1
delta_L_star_156 = L_star_TR_156 - L_star_Std_156;

%Calculating the differences between two coordinates of a^*_TR and a^*_Std
%by using
%Equation 5.2
delta_a_star_156 = a_star_TR_156 - a_star_Std_156;

%Calculating the differences between two coordinates of b^*_TR and b^*_Std
%by using
%Equation 5.3
delta_b_star_156 = b_star_TR_156 - b_star_Std_156;

%Calculating the differences between two coordinates of C^*_ab,bat and 
%C^*_ab,Std by using
%Equation 5.3
delta_C_ab_star_156 = C_ab_star_TR_156 - C_ab_star_Std_156;

%Using the Seve formula to calculate delta_H_ab_star_A.
%Equation 5.8
delta_H_ab_star_156 = (a_star_TR_156.*b_star_Std_156- a_star_Std_156.*b_star_TR_156)./...
    ((0.5.*(C_ab_star_TR_156.*C_ab_star_Std_156+a_star_TR_156.*a_star_Std_156+b_star_TR_156.*b_star_Std_156))...
    .^(1/2));

%Calculating color differences between two colors by using the method of Euclidean distance
%Equation 5.9
delta_E_ab_156 = (delta_L_star_156.^2+delta_a_star_156.^2+delta_b_star_156.^2).^(1/2);

%Calculate the color differences by using CIEDE2000
%Equation 5.35
delta_E_00_156 = deltaE00(CIELAB_TR_156,CIELAB_Std_156)';

%Calculate the color differences by using CIE94
%Equation 5.30
delta_E_94_156 = deltaE94(delta_L_star_156,C_ab_star_Std_156,C_ab_star_TR_156,delta_C_ab_star_156,delta_H_ab_star_156);

%Calculating STRESS for Eab, E00, and E94.
stress_delta_E_ab_156 = stress(delta_E_ab_156,delta_V_156);
stress_delta_E_00_156 = stress(delta_E_00_156,delta_V_156);
stress_delta_E_94_156 = stress(delta_E_94_156,delta_V_156);
stress3_156 = [stress_delta_E_ab_156 stress_delta_E_94_156 stress_delta_E_00_156];

%% Question 1

%Function to implement STRESS by using equations 5.20 and 5.21
function s = stress(deltaE,deltaV)
    F = sum(deltaE.^2)./sum(deltaE.*deltaV);
    s = 100*(sum((deltaE-(F.*deltaV)).^2)./sum((F.^2).*(deltaV.^2))).^(1/2);
end

%% Question 2

%Function to implement delta E94 by using equations 5.30-5.34
function E94 = deltaE94(dL,CS,CT,dC,dH)
    C_ab = sqrt(CS.*CT);
    S_L = 1;
    S_C = 1 + (0.045.*C_ab);
    S_H = 1 + (0.015.*C_ab);
    k_L = 1;
    k_C = 1;
    k_H = 1;
    E94 = (((dL)./(k_L.*S_L)).^2+((dC)./(k_C.*S_C)).^2+...
        ((dH)./(k_H.*S_H)).^2).^(1/2);
end




