%% Quantifying Fixational Eye Stability with BCEA(bivariate contour ellipse area) and ISOA(isoline area) methods
% Refecences:
%
% 1) Castet, E., & Crossland, M. (2012). Quantifying eye stability during a
% fixation task: a review of definitions and methods. Seeing and
% Perceiving, 25(5), 449-469.
%
% 2) Timberlake, G. T., Sharma, M. K., Grose, S. A., Gobert, D. V., Gauch,
% J. M., & Maino, J. H. (2005). Retinal location of the preferred retinal
% locus relative to the fovea in scanning laser ophthalmoscope images.
% Optometry and vision science, 82(3), E177-E187.
%
% 3)Zdravko Botev (2021). kernel density estimation
% (https://www.mathworks.com/matlabcentral/fileexchange/17204-kernel-density-estimation),
% MATLAB Central File Exchange. Retrieved November 20, 2021.
%
%
% 2021.11.20 by CH Du
% https://www.researchgate.net/profile/Chenghang-Du
%% A simple example
clc;clear;
load('S1_data.mat');
% delete missing data (blink)
delx = find(ismissing(xy_deg(1,:))==1);
dely = find(ismissing(xy_deg(2,:))==1);
del = union(delx,dely);
xy_deg(:,del)=[];
% default probability to calculate is 0.682.  
% default probability to plot is 0.682(red line) and 0.95(blue line).
% BCEA
figure;
[BCEA_formula,BCEA] = calBCEA(xy_deg,3500,10^-6,1); % biger n, closer between BCEA and BCEA_formula
disp((BCEA-BCEA_formula)/BCEA_formula); %nearly 0.05%
% ISOA
figure;
ISOA = calIsoA(xy_deg,2^11,10^-5,1);
% ISOA_fast 
% Program running time is saved by over 90%, the error is about 0.05% compared with ISOA method
figure;
ISOA_fast = calIsoA_fast(xy_deg,2^11,10^-5,1); % n has to be a power of 2
disp((ISOA_fast-ISOA)/ISOA);%nearly 0.05%