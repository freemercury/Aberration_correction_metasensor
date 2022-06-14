function Backprojection = backwardProjectGPU(psf_t,projection)

Backprojection=gpuArray.zeros(size(projection,1),size(projection,2),size(psf_t,3),'single');
[ra ca]=size(projection);
[rb cb]=size(psf_t(:,:,1));
r = ra+rb-1;c=ca+cb-1; p1 = (r-ra)/2;
b1 = gpuArray.zeros(r,r,'single');
a1 = gpuArray.zeros(r,r,'single');
z = 1;
a1(1:ra,1:ca) = projection(:,:) ;
b1(1:rb,1:cb) = psf_t(:,:,z) ;
clear con1;
con1 = ifft2(fft2(a1) .* fft2(b1));
Backprojection(:,:,z) = real(con1(p1+1:r-p1,p1+1:r-p1));

end