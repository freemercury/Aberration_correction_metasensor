%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  High-resolution reconstruction with wave-optics model.
%  The Code is created based on the method described in the following paper 
%   [1]  JIAMIN WU, YUDUO GUO and CHAO DENG etc,
%        An integrated imaging sensor for aberration-corrected 3D photography
%        Nature, 2022. 
% 
%   Contact: YUDUO GUO (gyd20@mails.tsinghua.edu.cn)
%   Date  : 6/10/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;close all;
%% Parameters
for num = [2 3 1]
    ang_R = 6; %% The diameter of the angle numbers is 2*ang_R+1
    para.initpsfpath = '../Data/PSF/psf_init.mat';  %% Input the initial PSF w/o aberration
    para.Metaimage_name = ['../Data/Metaimage/demo',num2str(num),'.tif'];
    para.savefolder = ['../Data/Reconstruction/Demo',num2str(num),'/'];
    para.Nshift = 5;
    para.Nnum = 15;
    para.anglenum = ang_R.^2;
    para.iter_ratio = 2.1;
    para.iter_times = 5;
    para.start_iteration = 1;
    para.ab_size = 205;
    para.x3objspace = [0]+1e-8;

    if ~exist(para.savefolder,'file')
        mkdir(para.savefolder);
    end
    META_subfunc(para);
end

