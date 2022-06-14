function PSF_inner_applyphase(aberration, x3objspace, savepathgen, iteration,mask_aber)
M =         1;
MLPitch =   69*1e-6;
Nnum =      15;
OSR =       3;
n =         1;
NA =        1/20;
fml =       MLPitch/(2*NA);
lambda =    525*1e-9;

eqtol = 1e-10;% threshold
tol = 0.0005; % threshold

k = 2*pi*n/lambda; % k
k0 = 2*pi*1/lambda; % k
d = fml;   % optical distance between the microlens and the sensor
ftl = 50e-3;        % focal length of tube lens

if mod(Nnum,2)==0
    error(['Nnum should be an odd number']);
end
pixelPitch = MLPitch/Nnum; %% pitch of virtual pixels

% define object space
% test PSF on line on the NATIVE image plane

x1objspace = [0];
x2objspace = [0];

objspace = ones(length(x1objspace),length(x2objspace),length(x3objspace));% discrete object space

IMGsize=Nnum*19;
IMGsize_L = Nnum*19;
HALF_ML_NUM = 8;

validpts = find(objspace>eqtol);% find non-zero points
numpts = length(validpts);%
[p1indALL p2indALL p3indALL] = ind2sub( size(objspace), validpts);% index to subcripts
p1ALL = x1objspace(p1indALL)';% effective obj points x location
p2ALL = x2objspace(p2indALL)';% effective obj points y location
p3ALL = x3objspace(p3indALL)';% effective obj points z location

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% projection from points on z axis %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Start Calculating PSF...']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Compute Light Field PSFs (light field) %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pixelPitch_OSR = MLPitch/OSR/Nnum; %simulated pixel size after OSR
fov=length(-(HALF_ML_NUM+1)*OSR*Nnum:(HALF_ML_NUM+1)*OSR*Nnum)*pixelPitch_OSR;   %the size of field of view for the PSF
pixelSize_OSR=length(-(HALF_ML_NUM+1)*OSR*Nnum:(HALF_ML_NUM+1)*OSR*Nnum); %the number of the pixels for each PSF
k2=2*pi/lambda;
% fx_sinalpha = lambda / 2 / pixelPitch2 / lambda;
fx_sinalpha = 1/(2*pixelPitch_OSR);
fx_step = 1/fov ;
fx_max = fx_sinalpha ;
fx= -fx_max+fx_step/2 : fx_step : fx_max;
[fxcoor fycoor] = meshgrid( fx , fx );
fx2coor=fxcoor.*fxcoor;
fy2coor=fycoor.*fycoor;

aperture_mask=((fx2coor+fy2coor)<=((NA/(lambda*M)).^2));
psfWAVE2=ones(pixelSize_OSR,pixelSize_OSR).*aperture_mask;
psfWAVE2 = gpuArray(single(psfWAVE2));

x1MLspace = (pixelPitch/OSR)* [-(Nnum*OSR-1)/2 : 1 : (Nnum*OSR-1)/2]; % total x space per ML
x2MLspace = (pixelPitch/OSR)* [-(Nnum*OSR-1)/2 : 1 : (Nnum*OSR-1)/2]; % total y space per ML
x1MLdist = length(x1MLspace);%
x2MLdist = length(x2MLspace);%

x1space = (pixelPitch/OSR)*[-(HALF_ML_NUM+1)*Nnum*OSR:1:(HALF_ML_NUM+1)*Nnum*OSR]; % x space
x2space = (pixelPitch/OSR)*[-(HALF_ML_NUM+1)*Nnum*OSR:1:(HALF_ML_NUM+1)*Nnum*OSR]; % y space
x1length = length(x1space);%x
x2length = length(x2space);%y

[MLARRAY,MLARRAYab] = calcML(fml, k0, x1MLspace, x2MLspace, x1space, x2space); % micro array phase mask
MLARRAY = gpuArray(single(MLARRAY));

x1objspace = (pixelPitch/M)*[-floor(Nnum/2):1:floor(Nnum/2)];% corresponding object space x1
x2objspace = (pixelPitch/M)*[-floor(Nnum/2):1:floor(Nnum/2)];% corresponding object space x2
XREF = ceil(length(x1objspace)/2);
YREF = ceil(length(x1objspace)/2);
centerPT = ceil(length(x1space)/2);
halfWidth = HALF_ML_NUM*Nnum*OSR;

CP = ( (centerPT-1)/OSR+1 - halfWidth/OSR :1: (centerPT-1)/OSR+1 + halfWidth/OSR  );%
H_z = zeros(length(CP),length(CP),Nnum,Nnum);
psf_z = zeros(IMGsize,IMGsize,Nnum,Nnum,'single');

filepath = savepathgen;
mkdir(filepath);

for eachpt = 1
    aa = tic;
    if(eachpt<0)
        continue;
    else
        disp(['calcu point #',num2str(eachpt),' ...............']);
        time_s = tic;
        p1 = p1ALL(eachpt); % object point #eachpt x
        p2 = p2ALL(eachpt);
        p3 = p3ALL(eachpt);
        
        timeWAVE = tic;
        tempP=k2*n*p3*realsqrt((1-(fxcoor.*lambda./n.*M).^2-(fycoor.*lambda./n.*M).^2).*aperture_mask);
        tempP = gpuArray(single(tempP));
        psfWAVE_fAFTERNO=psfWAVE2.*exp(1j*tempP);
        psfWAVE_AFTERNO=fftshift(ifft2(ifftshift(squeeze(psfWAVE_fAFTERNO))));
        tempaperture = 2*NA/lambda*MLPitch/Nnum/M/OSR*size(psfWAVE_fAFTERNO,1);
        
        padr = round((size(psfWAVE_fAFTERNO,1) - size(aberration,1))/2);
        aberration = padarray(aberration,[padr,padr],0,'both');
        mask_aber = padarray(mask_aber,[padr,padr],0,'both');
        psfWAVE_AFTERNO = ifft2(ifftshift(fftshift(fft2(psfWAVE_AFTERNO)).*exp(1i* aberration) .* mask_aber));
        
        tt = toc(timeWAVE);
        disp(['calcu psfWAVE take ',num2str(tt),' sec....']);
        timeFre = tic;
        for b1 = 1:length(x2objspace)
            for a1 = 1:length(x1objspace)
                timein = tic;
                psfSHIFT0= im_shift2_GPU(psfWAVE_AFTERNO, OSR*(a1-XREF), OSR*(b1-YREF) );%
                f1=fresnel2D_GPU(psfSHIFT0.*MLARRAY.*MLARRAYab, pixelPitch/OSR, fml,lambda);%
                f1= im_shift2_GPU(f1, -OSR*(a1-XREF), -OSR*(b1-YREF) );%
                [f1_AP_resize, x1shift, x2shift] = pixelBinning_GPU(abs(f1).^2, OSR);
                f1_CP = f1_AP_resize( CP - x1shift, CP-x2shift );
                H_z(:,:,a1,b1) = gather(f1_CP);%
                tt = toc(timein);
%                 disp(['calcu one point H take ',num2str(tt),' sec....']);
            end
        end
        tt = toc(timeFre);
        disp(['calcu H take ',num2str(tt),' sec....']);
        
        H4Dslice = H_z;
        H4Dslice(H4Dslice< (tol*max(H4Dslice(:))) ) = 0;% remove noise
        H_z = H4Dslice;
        
        disp(['normalize...NA_',num2str(NA)]);
        
        sss = H_z(:,:,(Nnum+1)/2,(Nnum+1)/2);
        for b2 = 1:length(x2objspace)
            for a2 = 1:length(x1objspace)
                H_z(:,:,a2,b2) = H_z(:,:,a2,b2)./sum(sss(:));
            end
        end
        %% OSR = 5 
        if OSR == 5
            temp = zeros(size(H_z,1)+1,size(H_z,2)+1,size(H_z,3),size(H_z,4));
            temp(2:end,2:end,:,:) = H_z;
            H_z = temp(1:end-1,1:end-1,:,:,:);
        end
        %%
        disp('split');
        border=fix(IMGsize_L/2)-fix(size(H_z,1)/2);
        blur_image=zeros(IMGsize_L,IMGsize_L,size(H_z,3),size(H_z,4));
        for i=1:size(H_z,3)
            for j=1:size(H_z,4)
                temp=zeros(IMGsize_L,IMGsize_L);
                temp(border+1:end-border,border+1:end-border)=squeeze(H_z(:,:,i,j));
                blur_image(:,:,i,j)=(im_shift3d(temp,i-((Nnum+1)/2),j-((Nnum+1)/2)));
            end
        end
        blur_image(isnan(blur_image)) = 0;
        %%
        IMGsize=size(H_z,1)-mod((size(H_z,1)-Nnum),2*Nnum);   %%%%% dertermine the psf size            
        sLF=zeros(IMGsize,IMGsize,Nnum,Nnum);
        index1=round(size(H_z,1)/2)-fix(size(sLF,1)/2);
        index2=round(size(H_z,1)/2)+fix(size(sLF,1)/2);
        for ii=1:size(H_z,3)
            for jj=1:size(H_z,4)
                sLF(:,:,ii,jj)=im_shift3(squeeze(H_z(index1:index2,index1:index2,ii,jj)),ii-((Nnum+1)/2), jj-(Nnum+1)/2);
            end
        end     
        bb=zeros(Nnum,Nnum,size(sLF,1)/size(H_z,3),size(sLF,2)/size(H_z,4),Nnum,Nnum);
        for i=1:size(H_z,3)
            for j=1:size(H_z,4)
                for a=1:size(sLF,1)/size(H_z,3)
                    for b=1:size(sLF,2)/size(H_z,4)
                        bb(i,j,a,b,:,:)=squeeze(sLF((a-1)*Nnum+i,(b-1)*Nnum+j,:,:));
                    end
                end
            end
        end
        WDF=zeros(  size(sLF,1),size(sLF,2),Nnum,Nnum  );
        for a=1:size(sLF,1)/size(H_z,3)
            for c=1:Nnum
                x=Nnum*a+1-c;
                for b=1:size(sLF,2)/size(H_z,4)
                    for d=1:Nnum
                        y=Nnum*b+1-d;
                        WDF(x,y,:,:)=squeeze(bb(:,:,a,b,c,d));
                    end
                end
            end
        end
        psf_z=WDF;
        save([savepathgen,'psf_pupil',num2str(iteration),'_layer_1.mat'],'psf_z','-v7.3');
        onetime = toc(aa);
        disp(['idz = ',num2str(eachpt),'taketime ',num2str(onetime),' sec......']);
    end
end