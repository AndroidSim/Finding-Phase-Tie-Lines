function varargout = bdy_parameterization(bdy,bdy_rep)
% b = bdy_parameterization(bdy,bdy_rep)
% inputs: bdy = matrix, bdy in cartesian coordinates [x(i) y(i)]
%         bdy_rep = boundary representation, either 'polygon' or 'spline'
% outputs: if bdy_rep = 'polygon', varargout{1} = parameter vector
%          if bdy_rep = 'spline', varargout{1} = parameter vector
%                                 varargout{2} = spline structure for x
%                                 varargout{3} = spline structure for y

[nbpts,bxy] = size(bdy);

if bxy > 3 | bxy < 2
	error('the coordinate dimension of the boundary pts must be 3 (ternary) or 2 (cartesian)');
else
    if bxy == 3
        bdy = tern2cart(bdy,1);
    end
end

switch bdy_rep
    case 'polygon'
        % bdy_length = sum(sqrt(diff(bdy_pts(:,1)).^2+diff(bdy_pts(:,2)).^2));
        bdy_length = sum(sqrt(sum(diff(bdy).^2,2)));
        idist(1,1) = 0;
        for i = 1:nbpts-1
            pt2pt_dist = norm(bdy(i+1,:)-bdy(i,:));
            idist(i+1,1) = idist(i,1) + pt2pt_dist; 
        end
        b = idist./bdy_length;
        varargout{1} = b;
    case 'spline'
        b = bdy_parameterization(bdy,'polygon');
%         spx = spap2(nbpts-3+1,3,b,bdy(:,1));
%         spy = spap2(nbpts-3+1,3,b,bdy(:,2));
        w = ones(size(b)); 
        w([1 end]) = 100;
        spx = spap2(10,4,b,bdy(:,1),w);
        spy = spap2(10,4,b,bdy(:,2),w);
%         spx = spap2(newknt(spx),4,b,bdy(:,1),w);
%         spy = spap2(newknt(spy),4,b,bdy(:,2),w);
%         x = fnval(spx,b);
%         y = fnval(spy,b);
%         x_prime = fnval(fnder(spx),b);
%         y_prime = fnval(fnder(spy),b);
%         x_2prime = fnval(fnder(spx,2),b);
%         y_2prime = fnval(fnder(spy,2),b); 
%         k = (x_prime.*y_2prime - y_prime.*x_2prime)./(((x_prime.^2) + (y_prime.^2)).^(3/2));
        varargout{1} = b;
        varargout{2} = spx;
        varargout{3} = spy;
    case 'ruled'
        
    otherwise
        error('boundary representation unknown');
end

return