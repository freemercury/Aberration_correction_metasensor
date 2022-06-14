function z_sli2q = sli2q(sx, sy, x, y)
%SLI2Q Spline-based Least-squares integration in quadrilateral geometry.
%   D * Z = G (G is mainly composed by spline estimated values).
%
%   Reference: 
%   L. Huang, J. Xue, B. Gao, C. Zuo, and M. Idir, "Spline based least 
%   squares integration for two-dimensional shape or wavefront 
%   reconstruction," Optics and Lasers in Engineering 91, 221-226 (2017)

%   Copyright since 2016 by Lei Huang. All Rights Reserved.
%   E-mail: huanglei0114@gmail.com
%   2016-10-01 Original Version
%   2016-11-01 Revised for x and y increasing directions.

% Check the number of arguments............................................
% Validate number of input arguments.
narginchk(4,4);
% Validate number of output arguments.
nargoutchk(1,1);

% Generate Matrix D and G..................................................
% Calculate size and ValidMask.
[Ny, Nx] = size(sx);
ValidMask = isfinite(sx) & isfinite(sy);

% Expand in x-direction.
sxEx = [NaN(Ny,1), sx, NaN(Ny,2)];
xEx  = [NaN(Ny,1),  x, NaN(Ny,2)];
syEx = [NaN(Ny,1), sy, NaN(Ny,2)];
yEx  = [NaN(Ny,1),  y, NaN(Ny,2)];

ExpandMaskx = isnan(sxEx);
se = [1 1 0 1 0];
DilatedExpandMaskx = imdilate(ExpandMaskx,se);
Maskx = DilatedExpandMaskx(:,2:end-2) & ~ExpandMaskx(:,2:end-2);

% Expand in x-direction.
sxEy = [NaN(1,Nx); sx; NaN(2,Nx)];
xEy  = [NaN(1,Nx);  x; NaN(2,Nx)];
syEy = [NaN(1,Nx); sy; NaN(2,Nx)];
yEy  = [NaN(1,Nx);  y; NaN(2,Nx)];

ExpandMasky = isnan(syEy);
se = [1;1;0;1;0];
DilatedExpandMasky = imdilate(ExpandMasky,se);
Masky = DilatedExpandMasky(2:end-2,:) & ~ExpandMasky(2:end-2,:);

% Compose matrices Dx and Dy.
Num = Ny*Nx;
ee = ones(Num,1);
Dx = spdiags([-ee,ee],[0,Ny],Num,Num);
Dy = spdiags([-ee,ee],[0, 1],Num,Num);

% Compose matrices Gx and Gy.
% O(h^5)
gx_x = (-1/13*sxEx(:,1:end-3)+sxEx(:,2:end-2)+sxEx(:,3:end-1)-1/13*sxEx(:,4:end)) ...
    .*(xEx(:,3:end-1)-xEx(:,2:end-2))*13/24;
gx_y = (-1/13*syEx(:,1:end-3)+syEx(:,2:end-2)+syEx(:,3:end-1)-1/13*syEx(:,4:end)) ...
    .*(yEx(:,3:end-1)-yEx(:,2:end-2))*13/24;
gy_x = (-1/13*sxEy(1:end-3,:)+sxEy(2:end-2,:)+sxEy(3:end-1,:)-1/13*sxEy(4:end,:)) ...
    .*(xEy(3:end-1,:)-xEy(2:end-2,:))*13/24;
gy_y = (-1/13*syEy(1:end-3,:)+syEy(2:end-2,:)+syEy(3:end-1,:)-1/13*syEy(4:end,:)) ...
    .*(yEy(3:end-1,:)-yEy(2:end-2,:))*13/24;

% O(h^3)
gx3_x = (sxEx(:,2:end-2)+sxEx(:,3:end-1)).*(xEx(:,3:end-1)-xEx(:,2:end-2))/2;
gx3_y = (syEx(:,2:end-2)+syEx(:,3:end-1)).*(yEx(:,3:end-1)-yEx(:,2:end-2))/2;
gy3_x = (sxEy(2:end-2,:)+sxEy(3:end-1,:)).*(xEy(3:end-1,:)-xEy(2:end-2,:))/2;
gy3_y = (syEy(2:end-2,:)+syEy(3:end-1,:)).*(yEy(3:end-1,:)-yEy(2:end-2,:))/2;

% Use O(h^3) values, if O(h^5) is not available.
gx_x(Maskx) = gx3_x(Maskx);
gx_y(Maskx) = gx3_y(Maskx);
gy_x(Masky) = gy3_x(Masky);
gy_y(Masky) = gy3_y(Masky);

Gx = gx_x + gx_y;
Gy = gy_x + gy_y;

% Compose D.
D = [Dx(isfinite(Gx),:); Dy(isfinite(Gy),:)];
clear Dx Dy;

% Compose matrix SpGx.
spGxQ = ComposeSpGx(x,sx,ValidMask,Nx,Ny);
% Compose matrix SpGy.
spGyQ = ComposeSpGy(y,sy,ValidMask,Nx,Ny);
clear sx sy x y;

% Replace with spline values, if available.
gy_y(end,:)=[];
gy_x(end,:)=[];
gx_x(isfinite(spGxQ)) = spGxQ(isfinite(spGxQ));
gy_y(isfinite(spGyQ)) = spGyQ(isfinite(spGyQ));

Gx = gx_x + gx_y;
Gy = gy_x + gy_y;

% Compose G.
G = [Gx(isfinite(Gx)); Gy(isfinite(Gy))];
clear Gx Gy;

% Solve "Warning: Rank deficient" for complete data by assuming Z(Ind)=0.  
Ind = find(D(1,:)==-1,1);
D(:,Ind) = [];
Z = D\G;
Z = [Z(1:Ind-1);0;Z(Ind:end)];

% Reconstructed result.
z_sli2q = reshape(Z,Ny,Nx);
z_sli2q(~ValidMask)= nan;

end




%% Subfunctions.

% Compose matrix spGx
function SpGx = ComposeSpGx(x,sx,ValidMask,Nx,Ny)

SpGx = NaN(Ny,Nx-1);
for ny = 1:Ny
    xl = x(ny,:)';
    vl = sx(ny,:)';
    
    % Check the number of sections.
    ml = ValidMask(ny,:)';   
    [Ns, Indices] = CheckSection(ml, Nx);
    
    % Spline fitting section by section.
    gs = cell(Ns,1);
    for ns = 1:Ns
        idx = Indices{ns};
        xx = xl(idx);
        vv = vl(idx);
        if length(xx)>1
            pp = spline(xx,vv); % "not-a-knot end condition"
            c = pp.coefs;
            switch(size(c,2))
                case 4  % 4 points for piecewise cubic spline fitting.
                    dx = diff(xx);
                    if sign(mean(dx))==1
                        gs{ns} = dx.*(c(:,4) + dx.*(c(:,3)./2 + dx.*(c(:,2)./3 + dx.*c(:,1)./4)));
                    else
                        dx = -flipud(dx);
                        gs{ns} = dx.*(c(:,4) + dx.*(c(:,3)./2 + dx.*(c(:,2)./3 + dx.*c(:,1)./4)));
                        gs{ns} = -flipud(gs{ns});
                    end     
                    
                case 3  % 3 points for 2nd order polynominal fitting.
                    % Here we do not use the polynomials. 
                    % We are going to use the Southwell expression instead
                    gs{ns} = diff(xx)*NaN; 
                
                case 2  % 2 points for 1st order polynominal fitting.
                    % Here we do not use the polynomials. 
                    % We are going to use the Southwell expression instead
                    gs{ns} = diff(xx)*NaN;
                
                case 1
                    % Logically impossible.
                    error('Only one point for fitting in x direction!');
                
                otherwise
                    % Logically impossible.
                    error('Unexpected number of points for fitting in x direction!');
            end
        end
    end
    sg = cat(1,gs{:});
    Valid = ml(1:end-1) & ml(1+1:end);
    pt = 1;
    for nx = 1 : Nx-1
        if Valid(nx) == 1
            SpGx(ny, nx) = sg(pt);
            pt = pt + 1;
        end
    end
end
end


% Compose matrix spGy
function SpGy = ComposeSpGy(y,sy,ValidMask,Nx,Ny)
SpGy = NaN(Ny-1,Nx);
for nx = 1:Nx
    yl = y(:,nx);
    vl = sy(:,nx);
    
    % Check the number of sections.
    ml = ValidMask(:,nx);
    [Ns, Indices] = CheckSection(ml, Ny);
    
    % Spline fitting section by section.
    gs = cell(Ns,1);
    for ns = 1:Ns
        idx = Indices{ns};
        yy = yl(idx);
        vv = vl(idx);
        if length(yy)>1
            pp = spline(yy,vv); % "not-a-knot end condition"
            c = pp.coefs;
            switch(size(c,2))
                case 4  % 4 points for piecewise cubic spline fitting.
                    dy = diff(yy);
                    if sign(mean(dy))==1
                        gs{ns} = dy.*(c(:,4) + dy.*(c(:,3)./2 + dy.*(c(:,2)./3 + dy.*c(:,1)./4)));
                    else
                        dy = -flipud(diff(yy));
                        gs{ns} = dy.*(c(:,4) + dy.*(c(:,3)./2 + dy.*(c(:,2)./3 + dy.*c(:,1)./4)));
                        gs{ns} = -flipud(gs{ns});
                    end
                    
                case 3  % 3 points
                    % Here we do not use the polynomials. 
                    % We are going to use the Southwell expression instead
                    gs{ns} = diff(yy)*NaN;

                case 2  % 2 points
                    % Here we do not use the polynomials. 
                    % We are going to use the Southwell expression instead
                    gs{ns} = diff(yy)*NaN;
                
                case 1
                    % Logically impossible.
                    error('Only one point for fitting in y direction!');
                
                otherwise
                    % Logically impossible.
                    error('Unexpected number of points for fitting in y direction!');
            end
        end
    end
    sg = cat(1,gs{:});
    Valid = ml(1:end-1) & ml(1+1:end);
    pt = 1;
    for ny = 1 : Ny-1
        if Valid(ny) == 1
            SpGy(ny, nx) = sg(pt);
            pt = pt + 1;
        end
    end
end
end


% Check Sections.
function [Ns, Indices] = CheckSection(ml, N)
if all(ml)==true      
    Ns = 1;
    Indices{Ns} = 1:N;
else
    Indices = cell(N,1);
    first = nan;
    last = nan;
    Ns = 0;
    for n = 1:N
        % Find the first.
        if n==1
            if ml(n)==true
                first = n;
            end
        else
            if ml(n)==true && ml(n-1)==false
                first = n;
            end
        end

        % Find the last.
        if n==N
            if ml(n)==true
                last = n;
            end
        else
            if ml(n)==true && ml(n+1)==false
                last = n;
            end
        end

        % Sum up the total number of sections and compose the Indices.
        if isfinite(first) && isfinite(last)
            Ns = Ns + 1;
            Indices{Ns} = first:last;
            first = nan;
            last = nan;
        end
    end
end
end
