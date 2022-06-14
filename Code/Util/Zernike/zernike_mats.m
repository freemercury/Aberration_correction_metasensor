function zernikeMatrices = zernike_mats(im,zernike_indices)
    
    % Returns an array of matrices corresponding to the Zernike polynomials
    % of degrees n and m as defined by the two columns of zernike_indices,
    % where the first column contains values of n and the second column
    % contains values of m. The dimensions of these matrices will be the
    % same as that of im.
    % This code is intended to be used as a helper function for
    % zernike_moments and zernike_recreation.
    %
    % Example:
    %     
    %     % Create some function z(x,y) over a 125x75 grid.
    %
    %     x=linspace(-1,1,75);
    %     y=linspace(1,-1,125);
    %     [x,y] = meshgrid(x,y);
    %     z = x.^2+y;
    %     figure,imagesc(elliptical_crop(z,1));
    %     colormap jet;
    %     title('z(x,y)');
    %
    %     % Specify the indices of Zernike polynomials of interest. In this
    %     % case, the indices correspond to the "piston", "x-tilt", "y-tilt",
    %     % and "power" Zernike polynomials.
    %
    %     indices = [0 0; 1 -1; 1 1; 2 0];
    %     zern_mats = zernike_mats(z,indices);
    %     for i = 1:size(zern_mats,3)
    %         figure, imagesc(elliptical_crop(zern_mats(:,:,i),1));
    %         colormap jet;
    %         title(indices(i,:));
    %     end
    %       
    %     % The images displayed show these four Zernike polynomials over
    %     % 125x75 grids.
    %
    % Functions required for use: zernike, zernike_radial, elliptical_crop 
    % See also: zernike_moments, zernike_recreation
    % 
    % Evan Czako, 8.14.2019
    % -------------------------------------------
    
    dim=size(im);
    x=linspace(-1,1,dim(2));
    y=linspace(1,-1,dim(1));
    [x,y] = meshgrid(x,y);
    [t,r] = cart2pol(x,y);
    
    zernikeMatrices = [];
    
    for i = 1:size(zernike_indices,1)
        zernikeMatrices(:,:,i) = zernike_gyd(r,t,zernike_indices(i,1),zernike_indices(i,2));
        zernikeMatrices(:,:,i) = elliptical_crop(zernikeMatrices(:,:,i),1);
    end
    
end