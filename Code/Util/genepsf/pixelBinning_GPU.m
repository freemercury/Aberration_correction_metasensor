function [OIMG, x1shift, x2shift] = pixelBinning_GPU(SIMG, OSR)

x1length = size(SIMG,1);
x2length = size(SIMG,2);

x1center = (x1length-1)/2 + 1;
x2center = (x2length-1)/2 + 1;

x1centerinit = x1center - (OSR-1)/2;
x2centerinit = x2center - (OSR-1)/2;
x1init = x1centerinit -  floor(x1centerinit/OSR)*OSR ;
x2init = x2centerinit -  floor(x2centerinit/OSR)*OSR ;


x1shift = 0;
x2shift = 0;
if x1init<1,
    x1init = x1init + OSR;
    x1shift = 1;
end
if x2init<1,
    x2init = x2init + OSR;
    x2shift = 1;
end


% SIMG_crop = SIMG( (x1init:1:end-OSR+1), (x2init:1:end-OSR+1) );
% SIMG_crop = SIMG_crop( (1:1: floor(size(SIMG_crop,1)/OSR)*OSR) ,  (1:1: floor(size(SIMG_crop,2)/OSR)*OSR) );
halfWidth = length( (x1init:x1center-1) );
% SIMG_crop = SIMG( [ (x1init:x1center-1) x1center x1center+1:x1center+halfWidth ],  [ (x2init:x2center-1) x2center x2center+1:x2center+halfWidth ] );
SIMG_crop = SIMG( x1init:x1center+halfWidth,  x2init:x2center+halfWidth  );


%%%%%%%%%%%%%%%%%% PIXEL BINNING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,n]=size(SIMG_crop); %M is the original matrix
SIMG_crop = sum( reshape(SIMG_crop,OSR,[]) ,1 );
SIMG_crop=reshape(SIMG_crop,m/OSR,[]).'; %Note transpose
SIMG_crop=sum( reshape(SIMG_crop,OSR,[]) ,1);
OIMG =reshape(SIMG_crop,n/OSR,[]).'; %Note transpose
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%