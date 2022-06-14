function radial = zernike_radial(r,n,m)

    % Returns an array corresponding to the radial component of the Zernike
    % polynomial of radial degree n and azimuthal degree m evaluated at all
    % points specified by matrix r.
    % 'r' is a matrix containing distances from the center of itself,
    % corresponding to an x-y grid satisfying 1>x>-1 and 1>y>-1. See
    % example for clarification.
    % Must satisfy n - m >= 0, n - m even, n >= 0, and m >= 0. n and m are
    % both integers.
    % Function elliptical_crop should be used in this context to ignore all 
    % data outside of the unit circle.
    % This function is utilized in function zernike to generate full 
    % Zernike polynomials.
    % Utilizes recursive algorithm described by Chong et al.
    % Pattern Recognition, 36;731-742 (2003).
    %
    % 
    % Example:
    %
    %     % Display the radial component of 'power' Zernike (n = 2, m = 0)  
    %     % over a 100x100 grid.
    %     x=linspace(-1,1,100);
    %     y=linspace(1,-1,100);
    %     [x,y] = meshgrid(x,y);
    %     [t,r] = cart2pol(x,y);
    %     zRad = zernike_radial(r,2,0);
    %     figure, imagesc(elliptical_crop(zRad,1));
    %     colormap jet;
    %
    %
    % Functions required for use: elliptical_crop
    % See also: zernike, zernike_moments, zernike_recreation
    % 
    % Evan Czako, 8.14.2019
    % -------------------------------------------

    if mod(n-m,2) == 1
        error('n-m must be even');
    end
    if n < 0 || m < 0
        error('n and m must both be positive in radial function')
    end
    if floor(n) ~= n || floor(m) ~= m
        error('n and m must both be integers')
    end
    if n == m
        radial = r.^n;
    elseif n - m == 2
        radial = n*zernike_radial(r,n,n)-(n-1)*zernike_radial(r,n-2,n-2);
    else
        H3 = (-4*((m+4)-2)*((m+4)-3)) / ((n+(m+4)-2)*(n-(m+4)+4));
        H2 = (H3*(n+(m+4))*(n-(m+4)+2)) / (4*((m+4)-1))  +  ((m+4)-2);
        H1 = ((m+4)*((m+4)-1) / 2)  -  (m+4)*H2  +  (H3*(n+(m+4)+2)*(n-(m+4))) / (8);
        radial = H1*zernike_radial(r,n,m+4) + (H2+H3 ./ r.^2).*zernike_radial(r,n,m+2);
        
        % Fill in NaN values that may have resulted from DIV/0 in prior
        % line. Evaluate these points directly (non-recursively) as they
        % are scarce if present.
        
        if sum(sum(isnan(radial))) > 0
            [row, col] = find(isnan(radial));
            c=1;
            while c<=length(row)
                x = 0;
                for k = 0:(n-m)/2
                    ((-1)^k*factorial(n-k))/(factorial(k)*factorial((n+m)/2-k)*factorial((n-m)/2-k))*0^(n-2*k);
                    x = x + ((-1)^k*factorial(n-k))/(factorial(k)*factorial((n+m)/2-k)*factorial((n-m)/2-k))*0^(n-2*k);
                end
                radial(row(c),col(c)) = x;
                c=c+1;
            end
        end

    end

end