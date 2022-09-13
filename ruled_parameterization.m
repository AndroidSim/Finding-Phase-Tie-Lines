function varargout = ruled_parameterization(coexpts,directrices,varargin)
% creates a ruled surface parameterization of the coexistence region.
% input: coexpts = matrix of points within coexistence region
%        directrices = cell array containing the directrices of ruled surface
%        gammas = cell array giving the gamma parameters for directrices
%        2 until the last directrix.
% output: coords = coordinates of parameterization

[np,nc] = size(coexpts);

if iscell(directrices)
    nd = length(directrices);
else
    error('1st argument must be a cell array of directrices');
end

if iscell(gammas)
    ng = length(gammas);
else
    error('2nd argument must be a cell array of gamma fxns');
end

if ng ~= nd-1
    error('the number of gamma fxns must equal the number of directrices minus 1');
end

directrix = directrices{1};
[ndp,ndc] = size(directrix);
v = sqrt(sum(diff(directrix).^2,2));
idist(1,1) = 0;
idist(2:ndp-1,1) = cumsum(v);
dlength = idist(end);
u = idist./dlength;
for d = 2:nd
    directrix = directrices{d};
    [ndp,ndc] = size(directrix);
    v = sqrt(sum(diff(directrix).^2,2));
    idist(1,1) = 0;
    idist(2:ndp-1,1) = cumsum(v);
    dlength = idist(end);       
%     bdy_length = sum(sqrt(sum(diff(bdy).^2,2)));
%     idist(1,1) = 0;
%     for i = 1:nbpts-1
%         pt2pt_dist = norm(bdy(i+1,:)-bdy(i,:));
%         idist(i+1,1) = idist(i,1) + pt2pt_dist; 
%     end
%     b = idist./bdy_length;
end

return