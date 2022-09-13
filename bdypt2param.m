function b = bdypt2param(bdy,pt,bdy_rep,varargin)
% b = bdypt2param(bdy,pt,bdy_rep,varargin)
% finds the parameter value b of boundary point pt on boundary bdy
% inputs: bdy = matrix, boundary in cartesian coordinates [x(i) y(i)]
%         pt = point(s) that lies somewhere on boundary or close to it.
%         bdy_rep = boundary representation, either 'polygon' or 'spline'
%         varargin = if bdy_rep = 'polygon', parameter values of bdy pts
%                    if bdy_rep = 'spline', spline structure of x and y
%                    values of boundary

[nbp,nbd] = size(bdy);
[np,npd] = size(pt);

if nbd ~= npd
    error('number of columns of 1st and 2nd argument must be equal');
end

x = pt(:,1);
y = pt(:,2);
switch bdy_rep
    case 'polygon'
        if isempty(varargin)
            bb = bdy_parameterization(bdy,bdy_rep);
        else
            bb = varargin{1};
        end
        maxx = max(bdy(:,1));
        minx = min(bdy(:,1));
        maxy = max(bdy(:,2));
        miny = min(bdy(:,2));
        x(x > maxx) = maxx;
        x(x < minx) = minx;
        y(y > maxy) = maxy;
        y(y < miny) = miny;
%         options = optimset('display','off');
%         zx = [];
%         n = 1;
%         for ui = 0:0.1:0.9
%             if sign(eqn2zero4x(ui)) == 0
%                 zx(n) = ui;
%                 n = n+1;
%                 continue;
%             elseif sign(eqn2zero4x(ui+0.1)) == 0
%                 zx(n) = ui+0.1;
%                 n = n+1;
%                 continue;
%             else
%                 if sign(eqn2zero4x(ui)) ~= sign(eqn2zero4x(ui+0.1))
%                     zx(n) = fzero(@eqn2zero4x,[ui ui+0.1],options);
%                     n = n+1;
%                     continue;
%                 else   
%                     continue;
%                 end
%             end
%         end
%         zy = [];
%         n = 1;
%         for ui = 0:0.1:0.9
%             if sign(eqn2zero4y(ui)) == 0
%                 zy(n) = ui;
%                 n = n+1;
%                 continue;
%             elseif sign(eqn2zero4y(ui+0.1)) == 0
%                 zy(n) = ui+0.1;
%                 n = n+1;
%                 continue;
%             else
%                 if sign(eqn2zero4y(ui)) ~= sign(eqn2zero4y(ui+0.1))
%                     zy(n) = fzero(@eqn2zero4y,[ui ui+0.1],options);
%                     n = n+1;
%                     continue;
%                 else   
%                     continue;
%                 end
%             end
%         end
%         if isempty(zx) | length(zx) > 2
%             error('error in finding zeros of x coordinate');
%         end
%         if isempty(zy) | length(zy) > 2
%             error('error in finding zeros of y coordinate');
%         end
%         
%         function v = eqn2zero4x(u) 
%             v = interp1(bb,bdy(:,1),u,'linear') - x;
%         end
%         function v = eqn2zero4y(u) 
%             v = interp1(bb,bdy(:,2),u,'linear') - y;
%         end
        
        bdy_length = sum(sqrt(sum(diff(bdy).^2,2)));
        [IN ON] = inpolygon(pt(:,1),pt(:,2),bdy(:,1),bdy(:,2));
        ion = find(ON);
        inon = find(~ON);
        vbdy = diff(bdy);
        for i = ion'
            v2bdy = repmat(pt(i,:),nbp,1)-bdy;
            for k = 1:nbp-1
        %         d(k) = norm(v2bdy(k,:));
                d2lin(k) = pt2lin_dist(pt(i,:),[bdy(k,:);bdy(k+1,:)]);
                u(k) = dot(v2bdy(k,:),vbdy(k,:))/dot(vbdy(k,:),vbdy(k,:));
            end
        %     ionk = find((d == 0 | d2lin == 0) & (u >= 0 & u < 1));
            ionk = find(round(d2lin.*1000)./1000 == 0 & (u >= 0 & u < 1));
            length2i = sum(sqrt(sum(diff([bdy(1:ionk,:);pt(i,:)]).^2,2)));
            b(i,1) = length2i/bdy_length;
            clear d u ionk;
        end
        for i = inon'
            v2bdy = repmat(pt(i,:),nbp,1)-bdy;
            for k = 1:nbp-1
                u(k) = dot(v2bdy(k,:),vbdy(k,:))/dot(vbdy(k,:),vbdy(k,:));
            end
            iu = find(u > 0 & u < 1);
            if isempty(iu)
                error(sprintf('point %d does not lie close enough to boundary',i));
            elseif length(iu) == 1
                ptk(i,:) = bdy(iu,:) + u(iu).*vbdy(iu,:);
            else
                for k = 1:length(iu)
                    w = v2bdy(iu(k),:);
                    nv = vbdy(iu(k),:)./norm(vbdy(iu(k),:));
                    d2bdy(k) = norm(w - dot(w,nv)*nv);
                end
                [temp,isort] = sort(d2bdy);
                iu = iu(isort(1));
                ptk(i,:) = bdy(iu,:) + u(iu).*vbdy(iu,:);
            end
            length2i = sum(sqrt(sum(diff([bdy(1:iu,:);ptk(i,:)]).^2,2)));
            b(i,1) = length2i/bdy_length;
            clear v2bdy u iu w nv d2bdy temp isort;
        end
    case 'spline'
        if isempty(varargin)
            [bb,spx,spy] = bdy_parameterization(bdy,bdy_rep);
        else
            spx = varargin{1};
            spy = varargin{2};
        end
        fnx = fncmb(spx,'*',-1);
        maxx = abs(fnmin(fnx));
        minx = fnmin(spx);
        fny = fncmb(spy,'*',-1);
        maxy = abs(fnmin(fny));
        miny = fnmin(spy);
        x(x > maxx) = maxx;
        x(x < minx) = minx;
        y(y > maxy) = maxy;
        y(y < miny) = miny; 
        if isscalar(x)
            fnx = fncmb(spx,'-',x);
            zx = fnzeros(fnx);
            zx = mean(zx);
            nzx = length(zx);    
        elseif isvector(x)
            for i = 1:length(x)
                fnx = fncmb(spx,'-',x(i));
                zx = fnzeros(fnx);
                vzx(i,:) = mean(zx);
                nzx(i,1) = length(vzx(i,:));
            end
        else
            error('x in bdypt2param must be a vector or scalar');
        end
        if isscalar(y)
            fny = fncmb(spy,'-',y);
            zy = fnzeros(fny);
            zy = mean(zy);
            nzy = length(zy);
        elseif isvector(y)
            for i = 1:length(y)
                fny = fncmb(spy,'-',y(i));
                zy = fnzeros(fny);
                vzy(i,:) = mean(zy);
                nzy(i,1) = length(vzy(i,:));
            end
        else
            error('y in bdypt2param must be a vector or scalar');
        end
        if isscalar(x) && isscalar(y) 
            for i = 8:-1:0
                tzx = round(zx.*10^i)./10^i;
                tzy = round(zy.*10^i)./10^i;
                [z,izx,izy] = intersect(tzx,tzy);
                if isempty(z)
                    continue;
                else
                    b = mean([zx(izx) zy(izy)]);
                    break;
                end
            end
        elseif isvector(x) && isvector(y)
            % length(x) should equal length(y), and both should equal np
            for k = 1:np
                for i = 8:-1:0
                    tzx = round(vzx(k,:).*10^i)./10^i;
                    tzy = round(vzy(k,:).*10^i)./10^i;
                    [z,izx,izy] = intersect(tzx,tzy);
                    if isempty(z)
                        continue;
                    else
                        b(k,1) = mean([vzx(k,izx) vzy(k,izy)]);
                        break;
                    end
                end
            end
        else
            error('x and y must both be a scalar or vector');
        end
    otherwise
        error('boundary representation unknown');
end

end