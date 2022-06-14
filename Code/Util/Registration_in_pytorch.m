function Registration_in_pytorch(para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nnum = para.Nnum;
anglenum = para.anglenum;
contrl_pointx = para.contrl_pointx;
contrl_pointy = para.contrl_pointy;
epoch = para.epoch;
lr = para.lr;
Metaimage_name = para.Metaimage_name;
Single_meta_imgname = para.Single_meta_imgname;
Savename = para.Savename;
Sidelobe = 15;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp = double(imread(Metaimage_name,((Nnum+1)/2-1)*Nnum+(Nnum+1)/2));
imwriteTFmeta(single(tmp),[Single_meta_imgname,'Meta_No_',num2str((Nnum+1)/2),'_',num2str((Nnum+1)/2),'.tif']);
for u = 1:Nnum
    parfor v = 1:Nnum
        if (u-(Nnum+1)/2)^2+(v-(Nnum+1)/2)^2 <= anglenum && (u ~= (Nnum+1)/2 || v ~= (Nnum+1)/2)
            tic;
            tmp = double(imread(Metaimage_name,(u-1)*Nnum+v));
            imwriteTFmeta(single(tmp),[Single_meta_imgname,'Meta_No_',num2str(u),'_',num2str(v),'.tif']);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            command = sprintf('python ./Lens_register.py %d %d %d %d %d %d %f %s %s %s', ...
                contrl_pointy,contrl_pointx,epoch,u,v,Sidelobe,lr,Single_meta_imgname,Savename);
            system(command);
            tt = toc;
            disp(['Registration with pytorch||u = ', num2str(u), ' & v = ', num2str(v),' has taken ',num2str(tt),' seconds!']);
        end
    end
end

end
