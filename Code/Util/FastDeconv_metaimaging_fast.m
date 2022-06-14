function Xguess = FastDeconv_metaimaging_fast(maxIter,meta_image,psf,lambda,savepath,AOstar,anglenum,para)
if AOstar ~= 1
    maxIter = 1;
end
Nnum = para.Nnum;
border = 3;
shift_threshhold = 50;
cor_sidelobe=shift_threshhold+10;
alpha = 2/3;
% R = round(Nnum/2);
[img_r,img_c,~] = size(meta_image);
[psf_r,~,~,~,allz] = size(psf);
weight = ones(Nnum,Nnum);
for i = 1:Nnum
    for j = 1:Nnum
        if ((i-round(Nnum/2))^2+(j-round(Nnum/2))^2 > anglenum)
            weight(i,j) = 0;
        end
    end
end
weight=squeeze(weight./sum(weight(:)));
seq = load(['./Util/seq',num2str(Nnum),'.mat']).seq;
psf = single(psf);
largepsf = gpuArray(single(zeros(img_r,img_c,allz,'single')));
[coordinate1,coordinate2]=meshgrid(1:size(meta_image,2),1:size(meta_image,1));
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
                        FFTPSF = fftn(flip((largepsf),3));
                        tmptmp = fftshift(fftshift(Xguess,1),2);
                        HXguess =sum((fftn(tmptmp).*FFTPSF),3)./allz;
                        HXguess =abs(ifftn(HXguess));
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

    Xguess_iter = gpuArray.zeros(size(meta_image,1),size(meta_image,2),Nnum*Nnum);
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
                blur_image_uv(isnan(blur_image_uv)) = 1e-7;
                blur_image_uv(blur_image_uv<0) = 1e-7;
                tmp = gpuArray(single(psf(:,:,u,v,:)));
                y = blur_image_uv;
%                 y = padarray(blur_image_uv, [1 1]*R, 'replicate', 'both');
                Xguess = fast_deconv(gpuArray(single(y)), tmp, lambda, alpha, gpuArray(single(y)));
                Xguess_tmp = Xguess.* weight(u,v);
                Xguess_tmp(isnan(Xguess_tmp)) = 1e-8;
                Xguess_tmp(Xguess_tmp<1e-8) = 0;
                Xguess_iter(:,:,(u-1)*Nnum+v) = real(Xguess_tmp);
%                 Xguess_iter(:,:,(u-1)*Nnum+v) = real(Xguess_tmp(R+1:end-R,R+1:end-R));
            end
            disp(['  iter ' num2str(i) ' | ' num2str(maxIter),' (u=',num2str(u), ', v=',num2str(v), '),  Energy=' num2str(sum(Xguess(:)))]);
        end
    end
    Xguess = sum(Xguess_iter,[3,4]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tt = gather(Xguess(border+1:end-border,border+1:end-border));
    output = tt./max(tt(:));
    imwriteTFmeta(single(output), strcat(savepath,'/Fast_AO_',num2str(AOstar),'_iter_',num2str(i),'_lambda_',num2str(lambda),'.tif'));
end
end

function [yout] = fast_deconv(yin, k, lambda, alpha, Init)
beta = 1;
beta_rate = 2*sqrt(2);
beta_max = 2^8;

% number of inner iterations per outer iteration
mit_inn = 1;

[m,n] = size(yin);
% initialize with input
yout = Init;

% padding zeros if H is not odd sized
% & make sure H is a 3 x 3 cell array
if ((mod(size(k, 1), 2) ~= 1) || (mod(size(k, 2), 2) ~= 1))
    fprintf('Error - blur kernel k must be odd-sized.\n');
end
ks = floor((size(k, 1)-1)/2);

% compute constant quantities
% see Eqn. (3) of paper
[Nomin1, Denom1, Denom2] = computeDenominator(gather(yin), gather(k));

% x and y gradients of yout (with circular boundary conditions)
% other gradient filters may be used here and their transpose will then need to
% be used within the inner loop (see comment below)

% store some of the statistics
costfun = [];
Outiter = 0;

youtx = [diff(yout, 1, 2), yout(:,1) - yout(:,n)];
youty = [diff(yout, 1, 1); yout(1,:) - yout(m,:)];
%% Main loop
while beta < beta_max


    Outiter = Outiter + 1;
    %     fprintf('Outer iteration %d; beta %.3g\n',Outiter, beta);

    gamma = beta/lambda;
    Denom = Denom1 + gamma*Denom2;
    Inniter = 0;
    stopci = 0;

    for Inniter = 1:mit_inn
        if (0)
            %%% Compute cost function
            youtk = conv2(yout, k, 'same');
            % likelihood term
            lh = sum(sum((youtk - yin).^2 ));

            if (alpha == 1)
                cost = (lambda/2)*lh +  sum(abs(youtx(:))) + sum(abs(youty(:)));
            else
                cost = (lambda/2)*lh +  sum(abs(youtx(:)).^alpha) + sum(abs(youty(:)).^alpha);
            end
            fprintf('Inniter iteration %d; cost %.3g\n', Inniter, cost);

            costfun = [costfun, cost];
        end
        %
        % w-subproblem: eqn (5) of paper
        %
        Wx = solve_image_fast_gyd(youtx, beta, alpha);
        Wy = solve_image_fast_gyd(youty, beta, alpha);

        %
        %   x-subproblem: eqn (3) of paper
        %
        % The transpose of x and y gradients; if other gradient filters
        % (such as higher-order filters) are to be used, then add them
        % below the comment above as well

        Wxx = [Wx(:,n) - Wx(:, 1), -diff(Wx,1,2)];
        Wxx = Wxx + [Wy(m,:) - Wy(1, :); -diff(Wy,1,1)];

        Fyout = (Nomin1 + gamma*fft2(Wxx))./Denom;
        yout = real(ifft2(Fyout));
    end %inner
    beta = beta*beta_rate;
end %Outer

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Nomin1, Denom1, Denom2] = computeDenominator(y, k)
%
% computes denominator and part of the numerator for Equation (3) of the
% paper
%
% Inputs:
%  y: blurry and noisy input
%  k: convolution kernel
%
% Outputs:
%      Nomin1  -- F(K)'*F(y)
%      Denom1  -- |F(K)|.^2
%      Denom2  -- |F(D^1)|.^2 + |F(D^2)|.^2
%

sizey = size(y);
otfk  = psf2otf(k, sizey);
Nomin1 = conj(otfk).*fft2(y);
Denom1 = abs(otfk).^2;
Denom2 = abs(psf2otf([1,-1],sizey)).^2 + abs(psf2otf([1;-1],sizey)).^2;
end
