function varargout = b2bdypt(b,bdy,bi)
% bdypt = b2bdypt(b,bdy,bi)
% linearly interpolates a pt on the boundary with parameter bi

bdypt = interp1(b,bdy,bi);
if isnan(bdypt)
    if bi <= b(1)
        bdypt = interp1(b,bdy,b(1));
        bracket{1} = 1;
    end
    if bi >= b(end)
        bdypt = interp1(b,bdy,b(end));
        bracket{1} = length(b);
    end
    if nargout < 2
        varargout{1} = bdypt;
    else
        varargout{1} = bdypt;
        varargout{2} = bracket;
    end
else
    if nargout < 2
        varargout{1} = bdypt;
    else
        for k = 1:length(bi)
            i = find(bi(k) == b);
            if isempty(i)
        %         il = find(b < bi);
        %         iu = find(b > bi);
        %         bracket = [il(end) iu(1)];
                il = 1;
                iu = length(b);
                while iu-il > 1
                    im = round((iu+il)/2);
                    if bi(k) > b(im)
                        il = im;
                    else
                        iu = im;
                    end
                end
                bracket{k} = [il iu];
            else
                bracket{k} = i;
            end
        end
        varargout{1} = bdypt;
        varargout{2} = bracket;
    end
end
    
return