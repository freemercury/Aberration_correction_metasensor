function [psfLine] = calcPSFFT(p3, fobj, NA, x1space, scale, lambda, fml, M, n)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
k = 2*pi*n/lambda;
alpha = asin(NA/n);
p1 = 0; 
p2 = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
psfLine = zeros(1,length(x1space));
parfor a=1:length(x1space)
    x1 = x1space(a);
    x2 = 0;        
    xL2normsq = (((x1+M*p1)^2+(x2+M*p2)^2)^0.5)/M;               
    v = k*xL2normsq*sin(alpha);    
    u = 4*k*p3*(sin(alpha/2)^2);
    
    Koi = M/((fobj*lambda)^2)*exp(-1i*u/(4*(sin(alpha/2)^2)));
    intgrand = @(theta) (sqrt(cos(theta))) .* (1+cos(theta))  .*  (exp(-(1i*u/2)* (sin(theta/2).^2) / (sin(alpha/2)^2)))  .*  (besselj(0, sin(theta)/sin(alpha)*v))  .*  (sin(theta));
    I0 = integral(@(theta)intgrand (theta),0,alpha);  
                
    psfLine(a) =  Koi*I0;
end
psfLine = abs(psfLine.^2)/max(abs(psfLine.^2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
