function [zernMoments] = zernike_moments(im,indices,mask)

    % Performs a least squares fit/decomposition of matrix im using Zernike 
    % polynomials of radial and azimuthal degrees described by indices, 
    % where the first column of indices contains values of n and the second 
    % column contains values of m.
    % Return a column vector containing the Zernike coefficients (or 
    % moments) that best describe matrix im.
    %
    % Example 1:
    % 
    %     % Create some function z(x,y) over a 100x100 grid.
    %     x=linspace(-1,1,100);
    %     y=linspace(1,-1,100);
    %     [x,y] = meshgrid(x,y);
    %     z = x.^2+y.^3;
    %     figure,imagesc(elliptical_crop(z,1));
    %     colormap jet;
    %     title('z(x,y)');
    % 
    %     % Specify desired indices of Zernike polynomials to be used in the
    %     % decomposition. In this case, include all polynomials up to and including 
    %     % radial degree 3.
    % 
    %     indices = [];
    %     for n = 0:3
    %         for m = -n:2:n
    %             indices = [indices; n m];
    %         end
    %     end
    % 
    %     % Display the coefficients associated with each Zernike polynomial used in
    %     % the fit and their corresponding (n,m) indices.
    %     disp(['    Coeff     ', 'n        ', 'm']);
    %     moments = zernike_moments(z,indices);
    %     disp([moments indices]);
    %
    %
    % Example 2:
    %
    %     % Some evidence that this function works.
    %     % Create list of all (n,m) indices up to n = 3:
    %     indices = [];
    %     for n = 0:3
    %         for m = -n:2:n
    %             indices = [indices; n m];
    %         end
    %     end
    % 
    %     % Create a 100x100 empty matrix for demonstration.
    %     zernike_sum_1 = zeros(100,100);
    % 
    %     % 3-dimensional matrix zern_mats contains the first 10 Zernike polynomials
    %     % as 100x100 grids.
    %     zern_mats = zernike_mats(zernike_sum_1,indices);
    % 
    %     % Add an equal amount of each of the first 10 Zernike polynomials to
    %     % zernike_sum_1. This way, the least-squares fit should return the same
    %     % coefficient for each of the polynomials used in the fit.
    %     for i = 1:size(zern_mats,3)
    %         zernike_sum_1 = zernike_sum_1 + zern_mats(:,:,i);
    %     end
    % 
    %     % The coefficients corresponding to each (n,m) pair used in the fit are the
    %     % same. In this case, they each have a value of 1.
    %     disp(['    Coeff     ', 'n        ', 'm']);
    %     disp([zernike_moments(zernike_sum_1,indices) indices]);
    %
    %     % In case you were wondering what an equal amount of Z1-10 looks like:
    %     figure,imagesc(zernike_sum_1);
    %     colormap jet;
    %
    %
    % Example 3:
    %
    %     % Some more evidence that this function works.
    %     % Create list of all (n,m) indices up to n = 3:
    %     indices = [];
    %     for n = 0:3
    %         for m = -n:2:n
    %             indices = [indices; n m];
    %         end
    %     end
    % 
    %     % Create a 100x100 empty matrix for demonstration.
    %     zernike_sum_2 = zeros(100,100);
    % 
    %     % 3-dimensional matrix zern_mats contains the first 10 Zernike polynomials
    %     % as 100x100 grids.
    %     zern_mats = zernike_mats(zernike_sum_2,indices);
    % 
    %     % Add an increasing amount of each of the first 10 Zernike polynomials to
    %     % zernike_sum_2. This way, the least-squares fit should return increasing 
    %     % values for the coefficients of each polynomial used in the fit.
    %     for i = 1:size(zern_mats,3)
    %         zernike_sum_2 = zernike_sum_2 + i*zern_mats(:,:,i);
    %     end
    % 
    % 
    %     % The coefficients corresponding to each (n,m) pair used in the fit increase 
    %     % incrementally from 1 to 10.
    %     disp(['    Coeff     ', 'n        ', 'm']);
    %     disp([zernike_moments(zernike_sum_2,indices) indices]);
    % 
    %     % In case you were wondering what this thing looks like:
    %     figure,imagesc(zernike_sum_2);
    %     colormap jet;
    %
    %
    % Functions required for use: zernike_mats, zernike, zernike_radial,
    % elliptical_crop
    % See also: zernike_recreation
    %
    % Evan Czako, 8.14.2019
    % -------------------------------------------
    
    dim=size(im);
    NumCols = dim(2);
    NumRows = dim(1);
    mask(mask == 0) = nan;
    
    zernikeMats = zernike_mats(im,indices).*mask;
 
    z_mats_reshaped = [];
    
    for i = 1:size(indices,1)
        z_mats_reshaped(:,i) = reshape(zernikeMats(:,:,i),NumCols*NumRows,1);
    end
    
    im(isnan(im))=0;   
    z_mats_reshaped(isnan(z_mats_reshaped))=0;

    image_reshaped = reshape(im,NumCols*NumRows,1);
    a = (z_mats_reshaped.'*z_mats_reshaped)^-1*(z_mats_reshaped.'*image_reshaped);
    zernMoments = a;
    
end