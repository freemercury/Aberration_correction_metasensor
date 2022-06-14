function zern = zernike_gyd(r,t,n,m)

    % Returns a matrix corresponding to the Zernike polynomial of radial
    % degree n and azimuthal degree m.
    % Matrices r and t correspond to the radii and angles of points on an 
    % x-y grid satisfying 1>x>-1 and 1>y>-1; they are obtained via 
    % cart2pol. See example for clarification.
    % Must satisfy n - m >= 0, n - m even, and n >= 0. n and m are 
    % both integers.
    % Function elliptical_crop should be used in this context to ignore all 
    % data outside of the unit circle.
    %
    % Example 1:
    %
    %     % Display the 'power' Zernike polynomial (n = 2, m = 0) over a 
    %     % 100x100 grid.
    %
    %     x=linspace(-1,1,100);
    %     y=linspace(1,-1,100);
    %     [x,y] = meshgrid(x,y);
    %     [t,r] = cart2pol(x,y);
    %     z_power = zernike(r,t,2,0);
    %     figure, imagesc(elliptical_crop(z_power,1));
    %     colormap jet;
    %
    % Example 2:
    %
    %     % Display all Zernike polynomials up to radial degree 4 over a
    %     % 100x100 grid.
    % 
    %     for n = 0:4
    %         for m = -n:2:n
    %             x=linspace(-1,1,100);
    %             y=linspace(1,-1,100);
    %             [x,y] = meshgrid(x,y);
    %             [t,r] = cart2pol(x,y);
    %             zern_mat = zernike(r,t,n,m);
    %             figure, imagesc(elliptical_crop(zern_mat,1));
    %             title(strcat('n = ',{' '},string(n),', m = ',{' '},string(m)));
    %             colormap jet;
    %         end
    %     end
    %
    %
    % Functions required for use: elliptical_crop, zernike_radial
    % See also: zernike_moments, zernike_recreation
    % 
    % Evan Czako, 8.14.2019
    % -------------------------------------------
    
    if mod(n-m,2) == 1
        error('n-m must be even');
    end
    if n < 0
        error('n must both be positive')
    end
    if floor(n) ~= n || floor(m) ~= m
        error('n and m must both be integers')
    end

    if m < 0
        zern = -zernike_radial(r,n,-m).*sin(m*t);
    else
        zern = zernike_radial(r,n,m).*cos(m*t);
    end
end