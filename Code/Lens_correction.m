%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Lens aberration correction.
%  The Code is created based on the method described in the following paper 
%   [1]  JIAMIN WU, YUDUO GUO and CHAO DENG etc,
%        An integrated imaging sensor for aberration-corrected 3D photography
%        Nature, 2022. 
% 
%    Contact: YUDUO GUO (gyd20@mails.tsinghua.edu.cn)
%    Date  : 6/10/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;
addpath(genpath(pwd));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Shiftmap_estimation_from_calibrated_meta_iamge
ang_R = 6;
para.Nnum = 15;
para.anglenum = 100;
para.anglenum2 = ang_R^2;
para.contrl_pointx = 20;
para.contrl_pointy = 15;
para.epoch = 80;
para.lr = 4e-4;
para.Metaimage_name = '../Data/Calibration/demo.tif';
para.Single_meta_imgname = '../Data/Calibration/TMP/Single/';
para.Savename = '../Data/Calibration/TMP/Calibrate_result/';
para.Phase_savename = '../Data/Calibration/Aberration/';
para.Shiftmap_namex = '../Data/Calibration/Shift_map/Mapx/';
para.Shiftmap_namey = '../Data/Calibration/Shift_map/Mapy/';
para.BFCorrect_Metaimage_name = '../Data/Metaimage/ISO_Whole_noART.tif';
para.AFCorrect_Metaimage_name = '../Data/Metaimage/Correct_ISO_Whole_noART.tif';
mkdir(para.Savename);
mkdir(para.Single_meta_imgname);
mkdir(para.Shiftmap_namex);
mkdir(para.Shiftmap_namey);
mkdir(para.Phase_savename);
Registration_in_pytorch(para);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Post_process_of_the_shift_map
weight = zeros(para.Nnum);
Shift_map = zeros(para.contrl_pointy,para.contrl_pointx,para.Nnum,para.Nnum,2);
count = 1;
for u = 1:para.Nnum
    for v = 1:para.Nnum
        if (u-(para.Nnum+1)/2)^2+(v-(para.Nnum+1)/2)^2 <= para.anglenum 
            weight(u,v) = 1;
            if u ~= (para.Nnum+1)/2 || v ~= (para.Nnum+1)/2
                name_img = [para.Savename,'Meta_',num2str(u),'_',num2str(v),'/',...
                    num2str(para.contrl_pointy),'_',num2str(para.contrl_pointx),...
                    '_',num2str(para.lr),'/Warped_Img/warped_Phantom_',num2str(para.epoch-1),'.tif'];
                stack_warp(:,:,count) = double(imread(name_img));
                count = count+1;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                name_map = [para.Savename,'Meta_',num2str(u),'_',num2str(v),'/',...
                    num2str(para.contrl_pointy),'_',num2str(para.contrl_pointx),...
                    '_',num2str(para.lr),'/Shift_map/ShiftMap_',num2str(para.epoch-1),'.mat'];
                Shift_map(:,:,u,v,:) = load(name_map).ShiftMap;
            end
        end
    end
end
imwriteTFmeta(single(stack_warp),'../Data/Calibration/Warp_result.tif');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dimension change of shift map
scale = 4;
Shift_map = imresize(Shift_map, scale, 'cubic');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
weight2 = zeros(size(weight));
for u = 1:para.Nnum
    for v = 1:para.Nnum
        if (u-(para.Nnum+1)/2)^2+(v-(para.Nnum+1)/2)^2 <= para.anglenum2 
            weight2(u,v) = 1;
        end
    end
end
Ratiox = size(stack_warp,2)/2;
Ratioy = size(stack_warp,1)/2;
Phase_map = zeros(para.Nnum,para.Nnum,2,para.contrl_pointy,para.contrl_pointx);
Phase_map_clean = zeros(size(Phase_map));
for uu = 1:para.contrl_pointy * scale
    for vv = 1:para.contrl_pointx * scale
        TMP = squeeze(Shift_map(uu,vv,:,:,:));
        Phase_map(:,:,:,uu,vv) = TMP;
%         Phase_map_clean(:,:,1,uu,vv) = weight2 .* TMP(:,:,1)*Ratiox;
%         Phase_map_clean(:,:,2,uu,vv) = weight2 .* TMP(:,:,2)*Ratioy;
        Phase_map_clean(:,:,1,uu,vv) = weight2 .* filloutliers(TMP(:,:,1)*Ratiox,'spline');
        Phase_map_clean(:,:,2,uu,vv) = weight2 .* filloutliers(TMP(:,:,2)*Ratioy,'spline');
        imwriteTFmeta(single(Phase_map_clean(:,:,1,uu,vv)),[para.Shiftmap_namex,'SMx_',num2str(uu),'_',num2str(vv),'.tif']);
        imwriteTFmeta(single(Phase_map_clean(:,:,2,uu,vv)),[para.Shiftmap_namey,'SMy_',num2str(uu),'_',num2str(vv),'.tif']);
        save('../Data/Calibration/Shift_map/Shift_map.mat','Phase_map_clean','-v7.3');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimate_Phase_map_for_Reconstruction
Nnum = para.Nnum;
phase_ratio = 4.6;
ps_smallsize = 99;
ps_phasesize = 205;
aberration = zeros(ps_phasesize,ps_phasesize,para.contrl_pointy* scale,para.contrl_pointx* scale);
indices = [];
for n = 0:8
    for m = -n:2:n
        indices = [indices; n m];
    end
end
desired_indices = 1:45;
zernikeMats = zernike_mats(zeros(ps_phasesize,ps_phasesize),indices);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:para.contrl_pointy * scale
    for jj = 1:para.contrl_pointx * scale
        map_wavshape(:,:,1) = Phase_map_clean(:,:,2,ii,jj);
        map_wavshape(:,:,2) = Phase_map_clean(:,:,1,ii,jj);
        expand = 5;
        waveShape = -squeeze(map_wavshape) .* weight2;
        waveShape_expand = zeros(expand*Nnum,expand*Nnum,2);
        waveShape_expand(:,:,1)=imresize(waveShape(:,:,1),[expand*Nnum,expand*Nnum],'nearest');
        waveShape_expand(:,:,2)=imresize(waveShape(:,:,2),[expand*Nnum,expand*Nnum],'nearest');
        [xx,yy] = meshgrid(-(expand*Nnum-1)/2:(expand*Nnum-1)/2,-(expand*Nnum-1)/2:(expand*Nnum-1)/2);
        r_actual = 2*sqrt(para.anglenum2)+1;
        mask_shiftmap = xx.^2+yy.^2 <= ((expand*r_actual/2).^2);
        waveShape_expand = waveShape_expand.*mask_shiftmap;
        waveShape_expand = waveShape_expand((end+1)/2-round(expand*r_actual/2)+1:(end+1)/2+round(expand*r_actual/2)-1,...
            (end+1)/2-round(expand*r_actual/2)+1:(end+1)/2+round(expand*r_actual/2)-1,:);
        [x1,y1] = meshgrid(1:size(waveShape_expand,1),1:size(waveShape_expand,2));
        [x2,y2] = meshgrid(linspace(1,size(waveShape_expand,1),ps_smallsize),linspace(1,size(waveShape_expand,1),ps_smallsize));
        calcu_dephase = zeros(ps_smallsize,ps_smallsize,2);
        calcu_dephase(:,:,1)  = interp2(x1,y1,waveShape_expand(:,:,1),x2,y2,'cubic');
        calcu_dephase(:,:,2)  = interp2(x1,y1,waveShape_expand(:,:,2),x2,y2,'cubic');
        ra = (ps_smallsize-1)/2;
        [xx,yy]=meshgrid(-ra:ra,-ra:ra);
        mask_phase = xx.^2+yy.^2<=(ra^2);
        [xx, yy] = meshgrid(linspace(-1,1,ps_smallsize),linspace(-1,1,ps_smallsize));
        calcu_phase = mask_phase.*sli2q(calcu_dephase(:,:,2).*mask_phase, calcu_dephase(:,:,1).*mask_phase, xx, yy);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        a1 = zernike_moments(phase_ratio*calcu_phase,indices,double(mask_phase));
        for i = 1:length(desired_indices)
            aberration(:,:,ii,jj) = aberration(:,:,ii,jj) + a1(desired_indices(i)).*zernikeMats(:,:,desired_indices(i));
        end
        aberration_tmp = aberration(:,:,ii,jj);

        mask_aber = aberration_tmp;
        mask_aber(~isnan(mask_aber)) = 1;
         
        aberration_tmp(isnan(aberration_tmp)) = 0;
        aberration_tmp = aberration_tmp - median(aberration_tmp,'all')-1e-4;
        aberration_tmp = aberration_tmp .* mask_aber;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%SAVE Image%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        filename = [para.Phase_savename,num2str(ii),'_',num2str(jj),'.png'];
        imshow(aberration_tmp,[]);colormap jet;caxis([-6*2*pi,6*2*pi]);
        set(gca,'color','none');set(gcf,'color','none');export_fig(filename,'-transparent');
    end
end