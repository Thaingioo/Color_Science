%Thain-Gi Oo
%Homework 11
%PCA Assignment

close all
clear all

%% load CCSG spectral data
[num,txt] = xlsread('ccsg.xlsx');
%reflectance factor
SGref = num(:,5:end);
%wavelength
SGwl = (380:10:730)';

[num,txt] = xlsread('all_1nm_data.xlsx');
cmf2 = interp1(num(:,1),num(:,6:8),SGwl);
D65 = interp1(num(:,1),num(:,3),SGwl);
SGXYZ = ((cmf2' * diag(D65) * SGref') ./ (cmf2(:,2)' * D65) )';
SGLab = xyz2lab(SGXYZ);


%% Question 1
%%%Use the Matlab function pca to compute 36 principal components 
%%%of the CCSG. Make a plot of the first 5 principal components versus
%%%wavelength, and label them with a legend. Make a plot of the cumulative 
%%%percentage of variance explained as a function of number of components.

%latent is the variance, coeff is the eigenvectors, mu is the mean, score
%is principal component scores 
[coeff,score,latent,tsquared,explained,mu] = pca(SGref);

%coeff is a 36 principal components of the CCSG.
coeff;

%Plotting the pca of CCSG (the first 5 principal components vs. wavelengths).
figure;
plot(SGwl,coeff(:,1:5));
hold on
yline(0);
xlabel('Wavelength (nm)');
ylabel('Reflectance factor');
title('First Five Eigenvectors', FontSize=16);
xlim([380 730]);
ylim([-0.5 0.5]);
legend('e_1','e_2','e_3','e_4','e_5');


%Finding the cumulative percentage of variance
CV = cumsum(latent);
CVP = (CV /CV(end,:))*100;

num_eig_vec = 1:36;
%Plotting the cumulative percentage of variance
figure;
plot(num_eig_vec,CVP, '.');
plot(num_eig_vec,CVP, '-b.', 'LineWidth',1, 'MarkerSize',15,...
    'MarkerEdgeColor','blue', 'MarkerFaceColor','blue');
xlabel('Number of Eigenvectors');
ylabel('Cummulative variance (%)');

%% Question 2
%%%Use the pca coefficients, scores, and mu (mean) to reconstruct the 140 spectra 
%%%of CCSG patches. Confirm that they match. 

%coeff is a pca coefficients
coeff;

%scores is a pca scores 
score;

% mu is mean 
mu;

%Reconstructing the data.
Recon_SGref = (score*coeff')+mu;

%% Question 3
%%%Reconstructing the 140 CCSG patches by using the first pca coefficient (vector)
%%%with the mu. 
Recon_SGref_1 = score(:,1)*coeff(:,1)'+mu;

%%%Compute the DE00 color difference (for 2-deg observer 
%%%under D65 illuminant) between the original and reconstructed; report the mean 
%%%and max DE00.

%Computing the reconstructed 140 CCSG patches of XYZ values
Recon_SGXYZ_1 = ((cmf2' * diag(D65) * Recon_SGref_1') ./ (cmf2(:,2)' * D65) )';

%Computing the  reconstructed 140 CCSG patches of Lab values
Recon_SGLab_1 = xyz2lab(Recon_SGXYZ_1);

%Computing deltaE00 values between the original and reconstructed 
%140 patches of their reflectance factors 
deltaE00_1 = deltaE00(Recon_SGLab_1',SGLab');
deltaE00_mean_1 = mean(deltaE00_1);
deltaE00_max_1 = max(deltaE00_1);

%% Question 4
%%%Reconstruct the 140 CCSG patches with the first two pca coefficients, and 
%report mean and max DE00.

%Reconstruct the 140 CCSG patches with the first two pca coefficients
Recon_SGref_2 = score(:,1:2)*coeff(:,1:2)'+mu;

%Computing the reconstructed 140 CCSG patches of XYZ values
Recon_SGXYZ_2 = ((cmf2' * diag(D65) * Recon_SGref_2') ./ (cmf2(:,2)' * D65) )';

%Computing the reconstructed 140 CCSG patches of Lab values
Recon_SGLab_2 = xyz2lab(Recon_SGXYZ_2 );

%Computing deltaE00 values between the original and reconstructed 
%140 patches of their reflectance factors 
deltaE00_2 = deltaE00(Recon_SGLab_2',SGLab');
deltaE00_mean_2 = mean(deltaE00_2);
deltaE00_max_2 = max(deltaE00_2);

%% Bonus
%Keep going with this series, and determine how many pca coefficients are 
%necessary to reduce the max DE00 below 0.5.

%Reconstruct the 140 CCSG patches with the first three pca coefficients

Recon_SGref_3 = score(:,1:3)*coeff(:,1:3)'+mu;

Recon_SGXYZ_3 = ((cmf2' * diag(D65) * Recon_SGref_3') ./ (cmf2(:,2)' * D65) )';
Recon_SGLab_3 = xyz2lab(Recon_SGXYZ_3 );

deltaE00_3 = deltaE00(Recon_SGLab_3',SGLab');
deltaE00_mean_3 = mean(deltaE00_3);
deltaE00_max_3 = max(deltaE00_3);

%Reconstruct the 140 CCSG patches with the first four pca coefficients

Recon_SGref_4 = score(:,1:4)*coeff(:,1:4)'+mu;

Recon_SGXYZ_4 = ((cmf2' * diag(D65) * Recon_SGref_4') ./ (cmf2(:,2)' * D65) )';
Recon_SGLab_4 = xyz2lab(Recon_SGXYZ_4 );

deltaE00_4 = deltaE00(Recon_SGLab_4',SGLab');
deltaE00_mean_4 = mean(deltaE00_4);
deltaE00_max_4 = max(deltaE00_4);

%Reconstruct the 140 CCSG patches with the first five pca coefficients

Recon_SGref_5 = score(:,1:5)*coeff(:,1:5)'+mu;

Recon_SGXYZ_5 = ((cmf2' * diag(D65) * Recon_SGref_5') ./ (cmf2(:,2)' * D65) )';
Recon_SGLab_5 = xyz2lab(Recon_SGXYZ_5 );

deltaE00_5 = deltaE00(Recon_SGLab_5',SGLab');
deltaE00_mean_5 = mean(deltaE00_5);
deltaE00_max_5 = max(deltaE00_5);

%Reconstruct the 140 CCSG patches with the first six pca coefficients

Recon_SGref_6 = score(:,1:6)*coeff(:,1:6)'+mu;

Recon_SGXYZ_6 = ((cmf2' * diag(D65) * Recon_SGref_6') ./ (cmf2(:,2)' * D65) )';
Recon_SGLab_6 = xyz2lab(Recon_SGXYZ_6 );

deltaE00_6 = deltaE00(Recon_SGLab_6',SGLab');
deltaE00_mean_6 = mean(deltaE00_6);
deltaE00_max_6 = max(deltaE00_6);

%Reconstruct the 140 CCSG patches with the first seven pca coefficients

Recon_SGref_7 = score(:,1:7)*coeff(:,1:7)'+mu;

Recon_SGXYZ_7 = ((cmf2' * diag(D65) * Recon_SGref_7') ./ (cmf2(:,2)' * D65) )';
Recon_SGLab_7 = xyz2lab(Recon_SGXYZ_7 );

deltaE00_7 = deltaE00(Recon_SGLab_7',SGLab');
deltaE00_mean_7 = mean(deltaE00_7);
deltaE00_max_7 = max(deltaE00_7);

%Reconstruct the 140 CCSG patches with the first eight pca coefficients

Recon_SGref_8 = score(:,1:8)*coeff(:,1:8)'+mu;

Recon_SGXYZ_8 = ((cmf2' * diag(D65) * Recon_SGref_8') ./ (cmf2(:,2)' * D65) )';
Recon_SGLab_8 = xyz2lab(Recon_SGXYZ_8);

deltaE00_8 = deltaE00(Recon_SGLab_8',SGLab');
deltaE00_mean_8 = mean(deltaE00_8);
deltaE00_max_8 = max(deltaE00_8);

%Reconstruct the 140 CCSG patches with the first nine pca coefficients

Recon_SGref_9 = score(:,1:9)*coeff(:,1:9)'+mu;

Recon_SGXYZ_9 = ((cmf2' * diag(D65) * Recon_SGref_9') ./ (cmf2(:,2)' * D65) )';
Recon_SGLab_9 = xyz2lab(Recon_SGXYZ_9);

deltaE00_9 = deltaE00(Recon_SGLab_9',SGLab');
deltaE00_mean_9 = mean(deltaE00_9);
deltaE00_max_9 = max(deltaE00_9);

%Reconstruct the 140 CCSG patches with the first ten pca coefficients

Recon_SGref_10 = score(:,1:10)*coeff(:,1:10)'+mu;

Recon_SGXYZ_10 = ((cmf2' * diag(D65) * Recon_SGref_10') ./ (cmf2(:,2)' * D65) )';
Recon_SGLab_10 = xyz2lab(Recon_SGXYZ_10);

deltaE00_10 = deltaE00(Recon_SGLab_10',SGLab');
deltaE00_mean_10 = mean(deltaE00_10);
deltaE00_max_10 = max(deltaE00_10);

%Reconstruct the 140 CCSG patches with the first eleven pca coefficients

Recon_SGref_11 = score(:,1:11)*coeff(:,1:11)'+mu;

Recon_SGXYZ_11 = ((cmf2' * diag(D65) * Recon_SGref_11') ./ (cmf2(:,2)' * D65) )';
Recon_SGLab_11 = xyz2lab(Recon_SGXYZ_11);

deltaE00_11 = deltaE00(Recon_SGLab_11',SGLab');
deltaE00_mean_11 = mean(deltaE00_11);
deltaE00_max_11 = max(deltaE00_11);

%The first eleven pca coefficients of the max DE00 is 0.3455.
%Therefore, there are 11 pca coefficients that are necessary to reduce 
%the max DE00 below 0.5.

%% OR
%We can use the for loop to find the number pca coefficients are 
%necessary to reduce the max DE00 below 0.5.
for i = 1:36
    Re_SGref = score(:,1:i)*coeff(:,1:i)'+mu;
    Re_SGXYZ = ((cmf2' * diag(D65) * Re_SGref') ./ (cmf2(:,2)' * D65) )';
    Re_SGLab = xyz2lab(Re_SGXYZ);
    deltaE00_iter = deltaE00(Re_SGLab',SGLab');
    deltaE00_mean_iter= mean(deltaE00_iter);
    deltaE00_max_iter = max(deltaE00_iter);
    if deltaE00_max_iter < 0.5
       fprintf('Iteration %d: Coefficients=1 to %d and deltaE00_max=%f\n', i, i, deltaE00_max_iter);
       break;
    end
end
