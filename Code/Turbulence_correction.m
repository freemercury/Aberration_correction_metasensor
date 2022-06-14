%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Turbulence correction with multi-site DAO.
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
gpuDevice(1);
Nshift= 5;
Nnum= 15;
scanning_pos = 1;
contrl_point1 = 60;
contrl_point2 = 80;
epcho = 80;
lr = 5e-4;
Metaimage_name  = '../Data/Metaimage/demo_turbulence.tif';
save_path       = '../Data/Reconstruction/Turbulence/';
save_path_single= [save_path,'singleWigner/'];
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
Meta_image = double(imread(Metaimage_name));
Meta_image = Meta_image./max(Meta_image(:));
%% Dynamic correction
index=zeros(size(Meta_image));
now_flag=mod(scanning_pos - 1,Nshift*Nshift);
index1=[a1(now_flag+1:end),a1(1:now_flag)];
index2=[a2(now_flag+1:end),a2(1:now_flag)];
for time_aver=1:1:Nshift*Nshift
    index(index1(time_aver):Nshift:end,index2(time_aver):Nshift:end) =0.99^(abs(time_aver-round(Nshift^2/2)));
end
%% Separate Meta image by time
znew = Dynamic_separate(Meta_image,index,index1,index2,Nshift);
imwriteTFmeta(single(znew),[save_path_single,'test_No_1_1.tif']);
%% Machine Learning based correction
command = sprintf('python Turbulence_register.py %d %d %d %d %d %d %f %s %s ', ...
    contrl_point1,contrl_point2,Nshift,0,0,epcho,lr,save_path_single,save_path);
system(command);
%% Apply distortion Map
Meta_image_TC = Meta_TurbCorrect(Meta_image,index,a1,a2,Nshift,1,1,save_path,contrl_point1,contrl_point2,epcho,lr);

