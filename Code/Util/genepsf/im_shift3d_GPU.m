function new_im = im_shift3d_GPU( img, SHIFTX, SHIFTY)
%UNTITLED2 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

eqtol = 1e-10;

xlength = size(img,1);
ylength = size(img,2);

if abs(mod(SHIFTX,1))>eqtol | abs(mod(SHIFTY,1))>eqtol
   error('SHIFTX and SHIFTY should be integer numbers');
end
% 
% SHIFTX = round(SHIFTX);
% SHIFTY = round(SHIFTY);

new_im = gpuArray(single(zeros(xlength, ylength,size(img,3))));

if SHIFTX >=0 & SHIFTY >= 0,
    new_im( (1+SHIFTX:end), (1+SHIFTY:end) ,:) = img( (1:end-SHIFTX), (1:end-SHIFTY) ,:);
elseif SHIFTX >=0 & SHIFTY < 0,
    new_im( (1+SHIFTX:end), (1:end+SHIFTY) ,:) = img( (1:end-SHIFTX), (-SHIFTY+1:end) ,:);
elseif SHIFTX <0 & SHIFTY >= 0,
    new_im( (1:end+SHIFTX), (1+SHIFTY:end) ,:) = img( (-SHIFTX+1:end), (1:end-SHIFTY) ,:);
else
    new_im( (1:end+SHIFTX), (1:end+SHIFTY) ,:) = img( (-SHIFTX+1:end), (-SHIFTY+1:end) ,:);
end


end

