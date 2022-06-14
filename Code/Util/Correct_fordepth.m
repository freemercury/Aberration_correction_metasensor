function Correct_fordepth(Shift_map,para)
Nnum = para.Nnum;
Meta_image_name = para.BFCorrect_Metaimage_name;
for i = 1:Nnum^2
    u = ceil(i/Nnum); v = mod(i,Nnum);
    if v == 0
        v = Nnum;
    end
    if (u-8)^2+(v-8)^2 <= para.anglenum2
        disp(['u = ', num2str(u), ', v = ', num2str(v), 'has been read!']);
        tmp = double(imread(Meta_image_name,i));
        tmp = tmp./max(tmp(:));
        [X,Y] = meshgrid(linspace(-1, 1, size(tmp,2)),linspace(-1, 1, size(tmp,1)));
        Shift_map_uv = imresize(Shift_map(:,:,:,u,v),[size(tmp,1),size(tmp,2)],'bicubic');
        tmp = interp2(X,Y,tmp,X+Shift_map_uv(:,:,1),Y+Shift_map_uv(:,:,2),'bicubic');
        Meta_image(:,:,u,v) = tmp;
    end
end
imwriteTFmeta(single(reshape(permute(Meta_image,[1,2,4,3]),[size(tmp,1),size(tmp,2),Nnum^2])),para.AFCorrect_Metaimage_name);

end