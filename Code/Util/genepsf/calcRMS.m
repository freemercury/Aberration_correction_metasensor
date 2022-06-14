function rms = calcRMS(phase)
w=[];
for i=1:size(phase,1)
    for j=1:size(phase,2)
        if ((i-round(size(phase,1)/2))^2+(j-round(size(phase,2)/2))^2)<= fix(size(phase,1)/2)^2
            w=[w,phase(i,j)];
        end
    end
end
rms=std(w);
end

