function Output_name = Meta_Calibrate(path,path_order,resize_ratio,rotate_ratio,save_path)
Output_name = [];
path2 = [path,'\*.tiff'];
namelist = dir(path2);
for idz = path_order
    count_num = 1;
    file = [path,'\',namelist(idz).name];
    if isempty(Output_name)
        Output_name = [save_path,namelist(idz).name];
    end
    info = imfinfo(file);
    Slice=size(info,1);
    Height = info.Height;
    Width = info.Width;
    if count_num == 1
        output = uint16(zeros(round(Height.*resize_ratio),round(Width.*resize_ratio),Slice));
    end
    for i = 1:Slice
        output(:,:,count_num) = imresize(imread(file,i),[round(Height.*resize_ratio),round(Width.*resize_ratio)],'bicubic');
        if rotate_ratio ~= 0
            output(:,:,count_num) = imrotate(output(:,:,count_num),rotate_ratio,'bicubic');
        end
        count_num = count_num+1;
    end
    disp([namelist(idz).name, ' has been Calibrated!']);
    imwriteTFmeta(uint16(output),[save_path,namelist(idz).name])
end

end

