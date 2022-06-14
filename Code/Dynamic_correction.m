%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Optical-flow-based Motion correction.
%  The Code is created based on the method described in the following paper 
%   [1]  JIAMIN WU, YUDUO GUO and CHAO DENG etc,
%        An integrated imaging sensor for aberration-corrected 3D photography
%        Nature, 2022. 
% 
%    Contact: YUDUO GUO (gyd20@mails.tsinghua.edu.cn)
%    Date  : 6/10/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath(pwd));
gpuDevice(1);
Nshift=5;
Nnum=15;
ang_R = 6;
scanning_pos = 17;
Metaimage_name   = '../Data/Metaimage/demo_dynamic.tif';
save_path        = '../Data/Reconstruction/Dynamic/';
save_path_single = [save_path,'singleWigner/'];
mkdir(save_path_single);
%% Preparameters
if Nshift==3
    a1= [1,1,1,2,3,3,3,2,2];
    a2= [3,2,1,1,1,2,3,3,2];
elseif Nshift==5
    a1=[1,1,1,1,1,2,3,4,5,5,5,5,5,4,4,4,4,3,2,2,3,3,3,2,2];
    a2=[5,4,3,2,1,1,1,1,1,2,3,4,5,5,4,3,2,2,2,3,3,4,5,5,4];
end
%% Input loading
Meta_image = double(imread(Metaimage_name,1));
%% Dynamic correction
index=zeros(size(Meta_image));
now_flag=mod(scanning_pos - 1,Nshift*Nshift);
index1=[a1(now_flag+1:end),a1(1:now_flag)];
index2=[a2(now_flag+1:end),a2(1:now_flag)];
for time_aver=1:1:Nshift*Nshift
    index(index1(time_aver):Nshift:end,index2(time_aver):Nshift:end) =0.99^(abs(time_aver-round(Nshift^2/2)));
end
znew = Dynamic_separate(Meta_image,index,index1,index2,Nshift);
imwriteTFmeta(single(znew),[save_path_single,'test_No0.tif']);
%% Optical flow based correction
alpha = 0.012;
ratio = 0.85;
minWidth = 10;
nOuterFPIterations = 7;
nInnerFPIterations = 1;
nSORIterations = 60;
para = [alpha,ratio,minWidth,nOuterFPIterations,nInnerFPIterations,nSORIterations];
[Meta_image_DC,WarpStack] = Meta_OF(Meta_image,index,index1,index2,Nshift,para);