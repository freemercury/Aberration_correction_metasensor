function znew= Meta_TurbCorrect(Meta_image,index,a1,a2,Nshift,u,v,save_name,contrl_point1,contrl_point2,epcho,lr)
[x,y] = find(index(1:Nshift,1:Nshift) == max(max(index(1:Nshift,1:Nshift))));
[X,Y] = meshgrid(linspace(-1, 1, size(Meta_image,2)/Nshift),linspace(-1, 1, size(Meta_image,1)/Nshift));

XX = zeros(size(Meta_image));
YY = zeros(size(Meta_image));
[X_out,Y_out] = meshgrid(linspace(-1, 1, size(Meta_image,2)),linspace(-1, 1, size(Meta_image,1)));

for i = 1:Nshift^2
    if a1(i) ~= x || a2(i) ~= y
        load([save_name,'Wigner_',num2str(u),'_',num2str(v),'/',num2str(i),'/',num2str(contrl_point1),...
            '_',num2str(contrl_point2),'_',num2str(lr),'/Shift_map/ShiftMap_',num2str(epcho-1),'.mat']);
        ShiftMap = imresize(ShiftMap,[size(Meta_image,1)/Nshift,size(Meta_image,2)/Nshift],'bicubic');
        XX(a1(i):Nshift:end,a2(i):Nshift:end) = X-ShiftMap(:,:,1);
        YY(a1(i):Nshift:end,a2(i):Nshift:end) = Y-ShiftMap(:,:,2);
    else
        XX(a1(i):Nshift:end,a2(i):Nshift:end) = X;
        YY(a1(i):Nshift:end,a2(i):Nshift:end) = Y;
    end
end
znew = griddata(XX,YY,Meta_image,X_out,Y_out,'cubic');

end
