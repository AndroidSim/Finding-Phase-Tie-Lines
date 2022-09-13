function bdy = interp_boundary(bdy,np,varargin)
% bdy = interp_boundary(bdy,np,varargin)
% inputs: bdy = matrix [x(i) y(i)], boundary represented as a set of points
%         np = scalar, number of points, >= 0.
%              if np > 1, output = interpolated boundary with np points
%              if np < 1, output = boundary point with parameter np
%         varargin = string, boundary representation, either linear or spline
% output: bdy = interpolated boundary or boundary point

if isempty(varargin)
    imethod = 'linear';
else
    if nargin > 3
        error('too many input arguments');
    else
        if ischar(varargin{1})
            imethod = varargin{1};
        else
            error('third argument must be a string specifying interpolation method');
        end
    end
end

if all(size(bdy) > 1) && ndims(bdy) == 2 % ie if matrix
    [nbp,nd] = size(bdy); % nbp = number of boundary pts, nd = number of dimensions
    
    if np > 1
        b = linspace(0,1,np)'; 
    else
        b = np;
    end
        
    if nd == 3
        bdy = tern2cart(bdy,1);
        [bdy_length,bdy_intervals] = bdy_fxn(bdy,'linear');
        bdy = interp1(bdy_intervals(:,1),bdy,b,imethod);
        bdy = cart2tern(bdy,1);
    elseif nd == 2
        [bdy_length,bdy_intervals] = bdy_fxn(bdy,'linear');
        bdy = interp1(bdy_intervals(:,1),bdy,b,imethod);
    else
        error('bdy must be an nx3 matrix (Xsm Xdopc Xchol) or nx2 matrix (X Y)');
    end  
else
    error('1st argument must be a matrix specifying boundary points');
end

return