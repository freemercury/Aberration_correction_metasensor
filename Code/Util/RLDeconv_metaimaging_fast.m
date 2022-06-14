function Xguess = RLDeconv_metaimaging_fast(maxIter,meta_image,psf,savepath,AOstar,anglenum,para)
Nnum = para.Nnum;
border = 3;
shift_threshhold = 50;
cor_sidelobe=shift_threshhold+10;
Xguess=ones(size(meta_image,1),size(meta_image,2));
Xguess=Xguess./sum(Xguess(:)).*sum(Xguess(:))./(size(Xguess,3)*size(Xguess,4));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[img_r,img_c,~] = size(meta_image);
[psf_r,~,~,~,allz] = size(psf);
if isempty(para.weight)
    weight = ones(Nnum,Nnum);
    for i = 1:Nnum
        for j = 1:Nnum
            if ((i-round(Nnum/2))^2+(j-round(Nnum/2))^2 > anglenum)
                weight(i,j) = 0;
            end
        end
    end
else
    weight = para.weight;
end
weight=squeeze(weight./sum(weight(:)));
iter_weight = 0.9*1/max(weight(:)); %% max iter weight
seq = load(['./Util/seq',num2str(Nnum),'.mat']).seq;
psf = single(psf);
largepsf = gpuArray(single(zeros(img_r,img_c,allz,'single')));
[coordinate1,coordinate2]=meshgrid(1:size(meta_image,2),1:size(meta_image,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:maxIter
    if AOstar == 1
        map_wavshape=zeros(Nnum,Nnum,2);
        if i>1
            for u_2=1:Nnum
                for v_2=1:Nnum
                    [v,u] = find((u_2-1)*Nnum+v_2 == seq);
                    if weight(u,v)==0 || sum(sum(psf(:,:,u,v)))==0
                        continue;
                    else
                        tmp = gpuArray(psf(:,:,u,v,:));
                        largepsf((img_r+1)/2-(psf_r-1)/2:(img_r+1)/2+(psf_r-1)/2,(img_c+1)/2-(psf_r-1)/2:(img_c+1)/2+(psf_r-1)/2,:) ...
                            = squeeze(tmp);
                        clear tmp;
                        FFTPSF = fft2(flip((largepsf),3));
                        tmptmp = fftshift(fftshift(Xguess,1),2);
                        HXguess =sum((fft2(tmptmp).*FFTPSF),3)./allz;
                        HXguess =abs(ifft2(HXguess));
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        sub_HXguess1=HXguess(1+cor_sidelobe:end-cor_sidelobe,1+cor_sidelobe:end-cor_sidelobe);
                        sub_blur_image=gpuArray(squeeze(meta_image(:,:,u,v)));
                        corr_map=gather(normxcorr2(sub_HXguess1,sub_blur_image));
                        [testa,testb]=find(corr_map==max(corr_map(:)));
                        map_wavshape(u,v,1)=-(size(sub_blur_image,1)-testa-cor_sidelobe);%%y
                        map_wavshape(u,v,2)=-(size(sub_blur_image,2)-testb-cor_sidelobe);%%x
                        disp(['u = ',num2str(u),',v = ',num2str(v),' has been AOed!']);
                    end
                end
            end
            cx=map_wavshape(round(Nnum/2),round(Nnum/2),1);
            cy=map_wavshape(round(Nnum/2),round(Nnum/2),2);
            map_wavshape(:,:,1)=(squeeze(map_wavshape(:,:,1))-cx);
            map_wavshape(:,:,2)=(squeeze(map_wavshape(:,:,2))-cy);
            for u=1:Nnum
                for v=1:Nnum
                    map_wavshape(u,v,1)=min(max(map_wavshape(u,v,1),-shift_threshhold),shift_threshhold);
                    map_wavshape(u,v,2)=min(max(map_wavshape(u,v,2),-shift_threshhold),shift_threshhold);
                end
            end
            save(strcat(savepath,'/shift_map_ite',num2str(i),'.mat'),'map_wavshape');
        end
    end

    for u_2=1:Nnum
        for v_2=1:Nnum
            [v,u] = find((u_2-1)*Nnum+v_2 == seq);
            if weight(u,v)==0 || ismember((u-1)*Nnum+v,para.list)
                continue;
            else
                if i>1 && AOstar == 1
                    map_wavshape_x=squeeze(map_wavshape(u,v,1));
                    map_wavshape_y=squeeze(map_wavshape(u,v,2));
                    blur_image_uv=interp2(coordinate1,coordinate2,gpuArray(single(meta_image(:,:,u,v))),coordinate1+map_wavshape_y,coordinate2+map_wavshape_x,'cubic',0);
                else
                    blur_image_uv=gpuArray(single(meta_image(:,:,u,v)));
                end
                tmp = gpuArray(single(psf(:,:,u,v,:)));
                largepsf((img_r+1)/2-(psf_r-1)/2:(img_r+1)/2+(psf_r-1)/2,(img_c+1)/2-(psf_r-1)/2:(img_c+1)/2+(psf_r-1)/2,:) ...
                    = squeeze(tmp);
                %%
                FFTPSF = fft2_new(largepsf);
                HXguess =fft2_new(Xguess).*FFTPSF;
                FFTPSF2 = fft2_new(rot90(largepsf,2));
                HXguessBack = ifft2_new(HXguess.*FFTPSF2);
                tmp = single(squeeze(blur_image_uv));
                errorBack = ifft2_new(fft2_new(tmp).*FFTPSF2);
                errorBack = real(errorBack./HXguessBack);
                Xguess = Xguess.*errorBack.*weight(u,v)*iter_weight+(1-weight(u,v)*iter_weight).*Xguess;
                Xguess(isnan(Xguess)) = 1e-8;
                Xguess(Xguess<1e-8) = 0;
                Xguess = real(Xguess);
            end
            disp(['  iter ' num2str(i) ' | ' num2str(maxIter),' (u=',num2str(u), ', v=',num2str(v), '),  Energy=' num2str(sum(Xguess(:)))]);
        end
    end
    tt = gather(Xguess(border+1:end-border,border+1:end-border));
    output = tt./max(tt(:));
    imwriteTFmeta(single(output), strcat(savepath,'/RL_AO_',num2str(AOstar),'_iter_',num2str(i),'.tif'));
end
end
% 
function output = fft2_new(input)
output = fftshift(fft2(ifftshift(input)));
end

function output = ifft2_new(input)
output = fftshift(ifft2(ifftshift(input)));
end


