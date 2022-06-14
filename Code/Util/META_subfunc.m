function META_subfunc(para)
%ps_phasesize = 2NA/lambda*d/Nnum/M/OSR*Phase_size[2*(HALF_ML_NUM+1)*OSR*Nnum+1]%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%Parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psf = load(para.initpsfpath).psf_z;
anglenum = para.anglenum;
Nshift = para.Nshift;
Nnum = para.Nnum;
iter_times = para.iter_times;
data_name = para.Metaimage_name;
Genpath = para.savefolder;
iter_ratio = para.iter_ratio;
ps_phasesize = para.ab_size;
ps_smallsize = 99;
para.list = 0;
para.weight = [];
maxIter_est_ab = 2;
maxIter_RL = 5;
weight=ones(Nnum,Nnum,2);
for u=1:Nnum
    for v=1:Nnum
        if ((u-round(Nnum/2))^2+(v-round(Nnum/2))^2) > anglenum
            weight(u,v,:)=0;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%Read in%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:Nnum^2
    u = ceil(i/Nnum);v = mod(i,Nnum);
    if v == 0
        v = Nnum;
    end
    if weight(u,v,1) == 1
        disp(['Demo || u = ', num2str(u), ', v = ', num2str(v), 'has been read!']);
        tmp = double(imread(data_name,i));
        if i == 1
            x_wdf = size(tmp,2);
            y_wdf = size(tmp,1);
            Meta_data = zeros(y_wdf,x_wdf,Nnum,Nnum);
        end
        if Nshift ~= Nnum
            tmp = imresize(tmp,Nnum/Nshift,'cubic');
        end
        Meta_data(:,:,u,v) = tmp;
    end
end
Meta_data = Meta_data./max(Meta_data(:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%Init Para%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aberration_iter = zeros(ps_phasesize);
AOstar = 1;
indices = [];
for n = 0:8
    for m = -n:2:n
        indices = [indices; n m];
    end
end
desired_indices = 4:45;
zernikeMats = zernike_mats(zeros(ps_phasesize,ps_phasesize),indices);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%Phase iteration%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for Iter_count = para.start_iteration:1:iter_times
    disp(['1Step -- Iter',num2str(Iter_count),'/',num2str(iter_times),':Reconstruction.']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    savefolder_recon = [Genpath,'/Iterations_',num2str(Iter_count),'/Recon/'];
    mkdir(savefolder_recon);

    if Iter_count ~= 1
        clear psf;
        PSFsavepathgen = [Genpath,'/Iterations_',num2str(Iter_count-1),'/PSF/'];
        psf = load([PSFsavepathgen,'psf_pupil',num2str(Iter_count-1),'_layer_1.mat']).psf_z;
    end

    if Iter_count == iter_times
        tmp_kernel = psf(:,:,round(Nnum/2),round(Nnum/2));
        [centerX,centerY] = find(tmp_kernel == max(tmp_kernel(:)));
        kernel = tmp_kernel(centerX-round(Nnum/2):centerX+round(Nnum/2),...
            centerY-round(Nnum/2):centerY+round(Nnum/2));
        for uu = 1:Nnum
            for vv = 1:Nnum
                if weight(uu,vv,1) == 1
                    for tape_count = 1:2
                        Meta_data(:,:,uu,vv) = edgetaper(Meta_data(:,:,uu,vv), kernel);
                    end
                end
            end
        end
        RLDeconv_metaimaging_fast(maxIter_RL,Meta_data,psf,savefolder_recon,0,anglenum,para);
        return;
    else
        lambda = 5e4;
        FastDeconv_metaimaging_fast(maxIter_est_ab,Meta_data,psf,lambda,savefolder_recon,AOstar,anglenum,para);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(['2Step -- Iter',num2str(Iter_count),'/',num2str(iter_times),':Phase Generation.']);
    expand = 5;
    savefolder_Phase = [Genpath,'/Iterations_',num2str(Iter_count),'/Phase/'];
    mkdir(savefolder_Phase);
    Mapsavename = [savefolder_recon,'shift_map_ite',num2str(maxIter_est_ab),'.mat'];
    map_wavshape = load(Mapsavename).map_wavshape;
    waveShape = -squeeze(map_wavshape) .* weight;
    waveShape_expand = zeros(expand*Nnum,expand*Nnum,2);
    waveShape_expand(:,:,1)=imresize(waveShape(:,:,1),[expand*Nnum,expand*Nnum],'nearest');
    waveShape_expand(:,:,2)=imresize(waveShape(:,:,2),[expand*Nnum,expand*Nnum],'nearest');
    [xx,yy] = meshgrid(-(expand*Nnum-1)/2:(expand*Nnum-1)/2,-(expand*Nnum-1)/2:(expand*Nnum-1)/2);
    r_actual = 2*sqrt(anglenum)+1;
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a1 = zernike_moments(iter_ratio*calcu_phase,indices,double(mask_phase));
    aber_tmp = zeros(ps_phasesize,ps_phasesize);
    for i = 1:length(desired_indices)
        aber_tmp = aber_tmp + a1(desired_indices(i)).*zernikeMats(:,:,desired_indices(i));
    end
    if Iter_count > 1
        clear aberration_iter;
        aberration_iter = load([Genpath,'/Iterations_',num2str(Iter_count-1),...
            '/Phase/pupil_ite',num2str(Iter_count-1),'_ratio_',num2str(iter_ratio),'.mat']).aberration_iter;
    end
    mask_aber = aber_tmp;
    mask_aber(~isnan(aber_tmp)) = 1;
    mask_aber(isnan(aber_tmp)) = 0;
    aber_tmp(isnan(aber_tmp)) = 0;
    aberration_iter = aberration_iter + aber_tmp;
    save([savefolder_Phase,'pupil_ite',num2str(Iter_count),'_ratio_',num2str(iter_ratio),'.mat'],...
        'a1','ps_smallsize','waveShape','aberration_iter');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(['3Step -- Iter',num2str(Iter_count),'/',num2str(iter_times),':PSF Generation.']);
    PSFsavepathgen = [Genpath,'/Iterations_',num2str(Iter_count),'/PSF/'];
    mkdir(PSFsavepathgen);
    PSF_inner_applyphase(aberration_iter, para.x3objspace, PSFsavepathgen, Iter_count, mask_aber)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end