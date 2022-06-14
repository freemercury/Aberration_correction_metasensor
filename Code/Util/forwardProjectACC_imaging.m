function TOTALprojection = forwardProjectACC_imaging( H, realspace, CAindex)

Nnum = size(H,3);
zerospace = gpuArray.zeros(size(realspace,1),size(realspace,2),'single');
TOTALprojection = zerospace;

for aa=1:Nnum
    for bb=1:Nnum
        for cc=1:size(realspace,3)
            Hs = squeeze(H(CAindex(cc,1):CAindex(cc,2),...
                CAindex(cc,1):CAindex(cc,2),...
                aa,bb,cc));          
            tempspace = zerospace;
            tempspace((aa:Nnum:end),(bb:Nnum:end)) = realspace((aa:Nnum:end),(bb:Nnum:end),cc);
            projection = conv2(tempspace, Hs, 'same');
            TOTALprojection = TOTALprojection + projection;            
        end
    end
end

