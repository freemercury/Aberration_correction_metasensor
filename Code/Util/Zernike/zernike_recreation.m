function [recreation,a_coeffs] = zernike_recreation(im,indices,desired_indices,mask2)
    
    % Returns a matrix that corresponds to the recreation of matrix im 
    % using Zernike polynomials of (n,m) values described by the columns of
    % indices, where the first column lists n values and the second column
    % lists m values.
    % desired_indices is a list of the rows of indices that are to be used
    % in the fit. For example, if only the first four (n,m) pairs listed in
    % indices are to be used in the recreation, desired_indices would be 
    % [1 2 3 4], or 1:4.
    %
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
    %     % Create list of all (n,m) pairs up to radial degree 3.
    % 
    %     indices = [];
    %     for n = 0:3
    %         for m = -n:2:n
    %             indices = [indices; n m];
    %         end
    %     end
    % 
    %     % Display recreations of z(x,y) using incrementally more Zernike
    %     % polynomials in each fit. As the polynomials become higher-order, the fits 
    %     % becomes more precise and are able to capture artifacts of higher spatial 
    %     % frequency.
    % 
    %     for i = 1:size(indices,1)
    %         figure,imagesc(zernike_recreation(z,indices,1:i));
    %         colormap jet;
    %         title(strcat('Recreation of z(x,y) using',{' '}, 'Z1-',string(i)));
    %     end
    % 
    %
    % Example 2:
    % 
    %     % Create some *other* function z(x,y) over a 100x100 grid.
    %     x=linspace(-1,1,100);
    %     y=linspace(1,-1,100);
    %     [x,y] = meshgrid(x,y);
    %     z = y.^2 + sin(6*x)+x;
    %     figure,imagesc(elliptical_crop(z,1));
    %     colormap jet;
    %     title('z(x,y)');
    % 
    %     % Create list of all (n,m) pairs up to radial degree 3.
    % 
    %     indices = [];
    %     for n = 0:3
    %         for m = -n:2:n
    %             indices = [indices; n m];
    %         end
    %     end
    % 
    %     % Display recreations of z(x,y) using incrementally more Zernike
    %     % polynomials in each fit. As the polynomials become higher-order, the fits 
    %     % becomes more precise and are able to capture artifacts of higher spatial 
    %     % frequency.
    % 
    %     for i = 1:size(indices,1)
    %         figure,imagesc(zernike_recreation(z,indices,1:i));
    %         colormap jet;
    %         title(strcat('Recreation of z(x,y) using',{' '}, 'Z1-',string(i)));
    %     end
    %
    %
    % Functions required for use: zernike_moments, zernike_mats, zernike,
    % zernike_radial, elliptical_crop
    %
    % Evan Czako, 8.14.2019
    % -------------------------------------------
    
    dim=size(im);
    NumCols = dim(2);
    NumRows = dim(1);
    
    zernikeMats = zernike_mats(im,indices);
    a_coeffs = zernike_moments(im,indices,mask2);

    recreation = zeros(NumRows,NumCols);   

    for i = 1:length(desired_indices)
        
        recreation = recreation + a_coeffs(desired_indices(i)).*zernikeMats(:,:,desired_indices(i)).* mask2;
        
    end
    
end