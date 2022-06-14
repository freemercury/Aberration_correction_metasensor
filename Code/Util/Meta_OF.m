function [znew,map] = Meta_OF(ref,index,a1,a2,Nshift,para)
ref = ref./max(ref(:));

[x,y] = find(index(1:Nshift,1:Nshift) == max(max(index(1:Nshift,1:Nshift))));
[X,Y] = meshgrid(linspace(1,size(ref,2)/Nshift,size(ref,2)/Nshift),linspace(1,size(ref,1)/Nshift,size(ref,1)/Nshift));

XX = zeros(size(ref));
YY = zeros(size(ref));
[X_out,Y_out] = meshgrid(linspace(1,size(ref,2)/Nshift,size(ref,2)),linspace(1,size(ref,1)/Nshift,size(ref,1)));

for i = 1:Nshift^2
    if a1(i) ~= x || a2(i) ~= y
        [vx,vy,~] = Coarse2FineTwoFrames(ref(a1(i):Nshift:end,a2(i):Nshift:end),ref(x:Nshift:end,y:Nshift:end),para);
        map(:,:,i,1) = vx;
        map(:,:,i,2) = vy;
        XX(a1(i):Nshift:end,a2(i):Nshift:end) = X+vx;
        YY(a1(i):Nshift:end,a2(i):Nshift:end) = Y+vy;
    else
        XX(a1(i):Nshift:end,a2(i):Nshift:end) = X;
        YY(a1(i):Nshift:end,a2(i):Nshift:end) = Y;
    end
end
znew = griddata(XX,YY,ref,X_out,Y_out,'cubic');

end