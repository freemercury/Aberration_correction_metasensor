function znew = Dynamic_separate(Meta_image,index,a1,a2,Nshift)

[x,y] = find(index(1:Nshift,1:Nshift) == max(max(index(1:Nshift,1:Nshift))));
znew(:,:,1) = Meta_image(x:Nshift:end,y:Nshift:end);
for i = 1:Nshift^2
    tmp = Meta_image(a1(i):Nshift:end,a2(i):Nshift:end);   
    znew(:,:,i+1) = tmp;
end

end