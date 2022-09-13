function varargout = chord_length(curve)
% clength = chord_length(curve)
% input: curve = matrix of points
% output: 1st vararg = clength = chord length
%         2nd vararg = parameterization = vector (nx1, n = #points)
% fxn calculates the chord length of curve and parameterization from first
% point to last point in matrix (nxd, n = #points, d = dimension)

% % bdy_length = sum(sqrt(diff(bdy_pts(:,1)).^2+diff(bdy_pts(:,2)).^2));
% bdy_length = sum(sqrt(sum(diff(bdy).^2,2)));
% idist(1,1) = 0;
% for i = 1:nbpts-1
%     pt2pt_dist = norm(bdy(i+1,:)-bdy(i,:));
%     idist(i+1,1) = idist(i,1) + pt2pt_dist; 
% end
% b = idist./bdy_length;
% varargout{1} = b;

[np,nd] = size(curve);

if np == 0
    error('curve has no points');
elseif np == 1
    error('curve has only 1 point');
elseif np == 2
    if all(diff(curve) == 0)
        clength = 0;
        c = [0;1];
    else
        v = sqrt(sum(diff(curve).^2,2));
        % idist(1,1) = 0;
        % idist(2:np-1,1) = cumsum(v);
        idist = [0;cumsum(v)];
        clength = idist(end);
        c = idist./clength;
    end
else
    v = sqrt(sum(diff(curve).^2,2));
    % idist(1,1) = 0;
    % idist(2:np-1,1) = cumsum(v);
    idist = [0;cumsum(v)];
    clength = idist(end);
    c = idist./clength;
end

if nargout == 1
    varargout{1} = clength;
elseif nargout == 2
    varargout{1} = clength;
    varargout{2} = c;
else
    error(nargoutchk(1,2,nargout));
end

return