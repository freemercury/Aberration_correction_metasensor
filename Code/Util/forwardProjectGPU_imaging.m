function TOTALprojection = forwardProjectGPU_imaging(H, realspace, OSR)

Nshift = size(H,3);
zerospace = gpuArray.zeros(size(realspace,1),size(realspace,2),'single');
TOTALprojection = zerospace;

for aa=1:Nshift
    for bb=1:Nshift
        for cc=1:size(realspace,3)
            Hs = squeeze(H(:,:,aa,bb,cc));          
            tempspace = zerospace;
            tempspace((aa:Nshift:end),(bb:Nshift:end)) = realspace((aa:Nshift:end),(bb:Nshift:end),cc);
            projection = conv2(tempspace, Hs, 'same');
            TOTALprojection = TOTALprojection + projection;            
        end
    end
end
% TOTALprojection = gydbinning(TOTALprojection,OSR);
end

