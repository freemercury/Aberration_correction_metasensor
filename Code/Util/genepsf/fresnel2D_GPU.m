function [f1]=fresnel2D_GPU(f0,dx0,z,lambda)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx = size(f0,1);
Ny = size(f0,2);
k = 2*pi/lambda;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
du = 1./(Nx*dx0);
u = gpuArray(single([0:ceil(Nx/2)-1 ceil(-Nx/2):-1]*du)); 
dv = 1./(Ny*dx0);
v = gpuArray(single([0:ceil(Ny/2)-1 ceil(-Ny/2):-1]*dv)); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% timeh1 = tic;
H = exp(-i*2*pi^2*(repmat(u',1,length(v)).^2+repmat(v,length(u),1).^2)*z/k); 
% tt = toc(timeh1);
% disp(['___calcu one point H1 take ',num2str(tt),' sec....']);
% H = gpuArray(single(H));
% timef1 = tic;
f1 = exp(i*k*z)*ifft2( fft2(f0) .* H );
% tt = toc(timef1);
% disp(['____calcu one point f1 take ',num2str(tt),' sec....']);
% f2 = ifft2( fft2(f0) .* H );
% dx1 = dx0;
% x1 = [-Nx/2:Nx/2-1]*dx1; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%












