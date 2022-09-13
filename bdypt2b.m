function varargout = bdypt2b(bdypt,b,bdy)
% b = bdypt2b(bdypt,bdy)
% finds the parameter value b of boundary point bdypt on boundary bdy
% bi = bdypt2b(b2bdypt(b,BDY,bi),b,BDY)

[nbp,nbd] = size(bdy);
[np,npd] = size(bdypt);

if length(b) ~= nbp
    error('length of b and # boundary points are not equal');
end
if nbd ~= npd
    error('number of columns of 1st and 2nd argument of bdypt2b must be equal');
end

warning off;
for i = 1:np
    %initialization
    bi(i,1) = NaN;
    matrix = bdy-repmat(bdypt(i,:),nbp,1);
    for p = 1:nbp-1 
        if all(matrix(p,:) == 0)
            bi(i,1) = b(p);
            bracket{i} = p;
            break;
        elseif all(matrix(p+1,:) == 0)
            bi(i,1) = b(p+1);
            bracket{i} = p+1;
            break;
        else
            for d = 1:nbd
                if sign(matrix(p,d)) ~= sign(matrix(p+1,d)) 
                    z(d) = 1;
                elseif sign(matrix(p,d)) == 0 & sign(matrix(p+1,d)) == 0
                    z(d) = 1;
                else
                    z(d) = 0;
                end
            end
            if all(z)
                for d = 1:nbd
                    ip(d) = (bdypt(i,d)-bdy(p,d))/(bdy(p+1,d)-bdy(p,d));
                end
                ip = nanmean(ip);
                bi(i,1) = b(p)+ip*(b(p+1)-b(p));
                bracket{i} = [p p+1];
                break; %continue
            end
        end
    end
    clear matrix;
end
warning on;

if nargout < 2
    varargout{1} = bi;
else
    varargout{1} = bi;
    varargout{2} = bracket;
end

return

% bdy_length = sum(sqrt(sum(diff(bdy).^2,2)));
% [IN ON] = inpolygon(bdypt(:,1),bdypt(:,2),bdy(:,1),bdy(:,2));
% ion = find(ON);
% inon = find(~ON);
% vbdy = diff(bdy);
% for i = ion'
%     v2bdy = repmat(bdypt(i,:),nbp,1)-bdy;
%     for k = 1:nbp-1
% %         d(k) = norm(v2bdy(k,:));
%         d2lin(k) = pt2lin_dist(bdypt(i,:),[bdy(k,:);bdy(k+1,:)]);
%         u(k) = dot(v2bdy(k,:),vbdy(k,:))/dot(vbdy(k,:),vbdy(k,:));
%     end
% %     ionk = find((d == 0 | d2lin == 0) & (u >= 0 & u < 1));
%     ionk = find(round(d2lin.*1000)./1000 == 0 & (u >= 0 & u < 1));
%     length2i = sum(sqrt(sum(diff([bdy(1:ionk,:);bdypt(i,:)]).^2,2)));
%     bi(i,1) = length2i/bdy_length;
%     clear d u ionk;
% end
% for i = inon'
%     v2bdy = repmat(bdypt(i,:),nbp,1)-bdy;
%     for k = 1:nbp-1
%         u(k) = dot(v2bdy(k,:),vbdy(k,:))/dot(vbdy(k,:),vbdy(k,:));
%     end
%     iu = find(u > 0 & u < 1);
%     if isempty(iu)
%         error(sprintf('point %d does not lie close enough to boundary',i));
%     elseif length(iu) == 1
%         ptk(i,:) = bdy(iu,:) + u(iu).*vbdy(iu,:);
%     else
%         for k = 1:length(iu)
%             w = v2bdy(iu(k),:);
%             nv = vbdy(iu(k),:)./norm(vbdy(iu(k),:));
%             d2bdy(k) = norm(w - dot(w,nv)*nv);
%         end
%         [temp,isort] = sort(d2bdy);
%         iu = iu(isort(1));
%         ptk(i,:) = bdy(iu,:) + u(iu).*vbdy(iu,:);
%     end
%     length2i = sum(sqrt(sum(diff([bdy(1:iu,:);ptk(i,:)]).^2,2)));
%     bi(i,1) = length2i/bdy_length;
%     clear v2bdy u iu w nv d2bdy temp isort;
% end

% options = optimset('Display', 'off','NonlEqnAlgorithm','gn');
% bi = fsolve(@eqn2zero,0.5,options);   
%     function v = eqn2zero(x)
%         if x < 0 
%             x = -x;
%         end
%         if x > 1
%             x = 1-x;
%         end    
%         v = b2bdypt(b,bdy,x)' - bdypt';
%     end
% end

% function y = findzero(b, c, x0)
% 
% options = optimset('Display', 'off'); % Turn off Display
% y = fsolve(@poly, x0, options);
% 
%     function y = poly(x) % Compute the polynomial.
%     y = x^3 + b*x + c;
%     end
% end