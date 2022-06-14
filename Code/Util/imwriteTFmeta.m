function [] = imwriteTF_new(Volume,name)
    Volume = squeeze(gather(Volume));
    
    if(numel(size(Volume)) == 2 || numel(size(Volume)) == 3)
        imwriteTFSKNormal(Volume,name);
    elseif(numel(size(Volume)) == 4)
        [a1,a2,a3,a4] = size(Volume);
        vol_t = zeros(a1,a2,a3*a4, class(Volume));
        for i3 = 1:a3
            for i4 = 1:a4
                vol_t(:,:,(i3-1)*a4 + i4) = Volume(:,:,i3,i4);
            end
        end
        imwriteTFSKNormal(vol_t,name);
    else
        warning('dimension of volume is not supported yet!');
    end
end




function [] = imwriteTFSKNormal(Volume,name)
    sz = whos('Volume');
    b = sz.bytes;
    if(b > 4 * 1000 * 1000 * 1000)
        lastdot_pos = find(name == '.', 1, 'last');
        parts = ceil(b / (4 * 1000 * 1000 * 1000));
        splitP = ceil(size(Volume, 3) / parts);
        for i = 1:parts
            nameSub = sprintf('%s.%d.%s', name(1:lastdot_pos-1), i, name(lastdot_pos+1:end));
            if(i < parts)
                imwriteTFSKSub(Volume(:,:,splitP * (i - 1) + 1:splitP * i), nameSub);
            else
                imwriteTFSKSub(Volume(:,:,splitP * (i - 1) + 1:end), nameSub);
            end
        end
    else
        imwriteTFSKSub(Volume,name);
    end
end



function [] = imwriteTFSKSub(Volume,name)
%imwriteTFSK 写入3D-tiff文件（支持single/double)
%    imwriteTFSK(volume,name)
t = Tiff(name,'w'); % Filename by variable name
tagstruct.ImageLength       = size(Volume,1);
tagstruct.ImageWidth        = size(Volume,2);
tagstruct.Photometric       = Tiff.Photometric.MinIsBlack;
% if(strcmp(class(Volume), 'single'))
if(isa(Volume, 'single'))
    tagstruct.BitsPerSample	= 32;
    tagstruct.SampleFormat	= Tiff.SampleFormat.IEEEFP;
elseif(isa(Volume, 'double'))
    warning('ImageJ may not support double/64-bit tiff!');
    tagstruct.BitsPerSample	= 64;
    tagstruct.SampleFormat	= Tiff.SampleFormat.IEEEFP;
elseif(isa(Volume, 'uint8'))
    tagstruct.BitsPerSample	= 8;
    tagstruct.SampleFormat	= Tiff.SampleFormat.UInt;
elseif(isa(Volume, 'uint16'))
    tagstruct.BitsPerSample	= 16;
     tagstruct.SampleFormat	= Tiff.SampleFormat.UInt;
end
tagstruct.SamplesPerPixel	= 1;
tagstruct.RowsPerStrip      = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software          = 'MATLAB';
setTag(t,tagstruct)
write(t,squeeze(Volume(:,:,1)));
for i=2:size(Volume,3) % Write image data to the file
    writeDirectory(t);
    setTag(t,tagstruct)
    write(t,squeeze(Volume(:,:,i))); % Append
end
close(t);
end

