function Ht = calcHt(H)

Hsize = size(H);
Hsize1 = Hsize(1);
Nnum = Hsize(3);
x3length = Hsize(5);

tmpsize = ceil(size(H,1)/Nnum);
if mod(tmpsize,2) == 1
    imgsize = (tmpsize+2)*Nnum;
else
    imgsize = (tmpsize+3)*Nnum;
end

zeroprojection = zeros(imgsize, imgsize);
imcenter = ceil(imgsize/2);
imcenterinit = imcenter - ceil(Nnum/2);

Ht = zeros(Hsize);
for aa=1:1:Nnum
    for bb=1:1:Nnum
        temp = zeroprojection;
        temp( imcenterinit+aa, imcenterinit+bb ) = 1;
        
        tempback = backwardProject(H, temp, Nnum) ;
        tempback_cut = tempback( ( imcenter - (Hsize1-1)/2 - 0*Nnum :imcenter + (Hsize1-1)/2 + 0*Nnum), ( imcenter - (Hsize1-1)/2 - 0*Nnum :imcenter + (Hsize1-1)/2 + 0*Nnum) , :);
        
        for cc=1:x3length,
            tempback_shift(:,:,cc) = im_shift2(tempback_cut(:,:,cc), (ceil(Nnum/2)-aa) , (ceil(Nnum/2)-bb) );
        end
        Ht(:,:,aa,bb,:) = tempback_shift;
    end       
end
