function pt = param2bdypt(bdy,b,bdy_rep,varargin)
% pt = param2bdypt(bdy,b,bdy_rep)
% inputs: bdy = matrix, boundary in cartesian coordinates [x(i) y(i)]
%         b = parameter
%         bdy_rep = boundary representation, either 'polygon' or 'spline'
%         varargin = if bdy_rep = 'polygon', parameter values of bdy pts
%                    if bdy_rep = 'spline', spline structure of x and y
%                    values of boundary

switch bdy_rep
    case 'polygon'
        if isempty(varargin)
            bb = bdy_parameterization(bdy,bdy_rep);
            pt = interp1(bb,bdy,b,'linear');
        else
            bb = varargin{1};
            pt = interp1(bb,bdy,b,'linear');
        end
    case 'spline'
        if isempty(varargin)
            [bb,spx,spy] = bdy_parameterization(bdy,bdy_rep);
            x = fnval(spx,b);
            y = fnval(spy,b);
            pt = [x y];
        else
            spx = varargin{1};
            spy = varargin{2};
            x = fnval(spx,b);
            y = fnval(spy,b);
            pt = [x y];
        end
    otherwise
        error('boundary representation unknown');
end

return