
datapath=[pwd,filesep,'results',filesep];%pwd is the current work directory
addpath(datapath); 
% load('ArabidopsisGOAinternotH201610311_NFA_bp.mat');
% load('ArabidopsisGOAinternotH201610311_NFA_cc.mat');
% load('ArabidopsisGOAinternotH201610311_NFA_mf.mat');
% load('YeastGOAinternotH201610311_NFA_bp.mat');
% load('YeastGOAinternotH201610311_NFA_cc.mat');
load('YeastGOAinternotH201610311_NFA_mf.mat');

precision1=precision;
% load('ArabidopsisGOAinternotH201610311_NoisyGOA_EC_hr_bp.mat');
% load('ArabidopsisGOAinternotH201610311_NoisyGOA_EC_hr_cc.mat');
% load('ArabidopsisGOAinternotH201610311_NoisyGOA_EC_hr_mf.mat');
% load('YeastGOAinternotH201610311_NoisyGOA_EC_hr_bp.mat');
% load('YeastGOAinternotH201610311_NoisyGOA_EC_hr_cc.mat');
load('YeastGOAinternotH201610311_NoisyGOA_EC_hr_mf.mat');



for ii=17:19
 a=precision1{ii};
 b=precision{ii};
[h,sig,ci]=ttest2(a,b);
fprintf('%d \n',h);
end
