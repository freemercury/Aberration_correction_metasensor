function f = calc_f(map_wavshape)
[rr,cc,~] = size(map_wavshape);
dfx=squeeze(map_wavshape(:,:,1));
dfy=squeeze(map_wavshape(:,:,2));
ra = (rr-1)/2;

[X,Y]=meshgrid([-ra:ra],[-ra:ra]);
%[X1,Y1]=meshgrid([-ra:0.1:ra],[-ra:0.1:ra]);
dfx( (X.^2+Y.^2)>ra^2 )=0;
dfy( (X.^2+Y.^2)>ra^2 )=0;
f=zeros(size(map_wavshape,1),size(map_wavshape,2));
Nnum=rr;
%dx = 1.4*2/(Nnum*525*1e-9);
dx = 1;
%ddx = 1.3082e-7*5;
ddx = 1;
ite = 10000;
for u=1:Nnum
    for v=1:Nnum
        temp1=0;
        temp2=0;
        
        record = 0;
        for i = 1:ite
            step = zeros(1,u+v-2);
            sr = v-1;
            sd = u-1;
            cs=randperm(length(step));
            step(cs(1:sr)) = 1;
            r2 = 0;
            rr = 1;
            cc = 1;
            for ids = 1:length(step)
                if(step(ids)==1)
                    r2 = r2+dfy(rr,cc)*dx*ddx;
                    cc = cc+1;
                else
                    r2 = r2+dfx(rr,cc)*dx*ddx;
                    rr = rr+1;
                end
            end
            record = record+r2;
        end
        
        f(u,v) = record/ite;
        disp(['record = ',num2str(record/ite)]);
        disp(['u = ',num2str(u),'; v = ',num2str(v)]);
        
        %{
        for kk=1:u
            temp1=temp1+dfx(kk,1)*dx*ddx;
        end
        for kk=2:v
            temp1=temp1+dfy(u,kk)*dx*ddx;
        end
        for kk=1:v
            temp2=temp2+dfy(1,kk)*dx*ddx;
        end
        for kk=2:u
            temp2=temp2+dfx(kk,v)*dx*ddx;
        end
        f(u,v)=(temp1+temp2)/2;
        f1(u,v)=temp1;
        f2(u,v)=temp2;
        %}
    end
end
%f=f-f(7,7);
f( (X.^2+Y.^2)>ra^2 )=0;
end

