function [parlines]=gen_par_lines(xmin,xmax,ymin,ymax,slope,div,npts)
% the function gen_par_lines generates a set of parallel lines with slope 'slope'
% within a box defined by xmin,xmax,ymin,ymax
% default values
vert_inc=1;
horz_inc=1;

if (slope == 0)
    ipar=1;
    for yi=ymin:vert_inc:ymax
        xline=linspace(xmin,xmax,npts);
        yline=yi;
        parlines(ipar,1:npts,1)=xline;
        parlines(ipar,1:npts,2)=yline;
        ipar=ipar+1;
    end
elseif (slope == Inf)
    ipar=1;
    for xi=xmin:horz_inc:xmax
        xline=xi;
        yline=linspace(ymin,ymax,npts);
        parlines(ipar,1:npts,1)=xline;
        parlines(ipar,1:npts,2)=yline;
        ipar=ipar+1;
    end
elseif (slope > 0)
    if (slope == div)%num2str(slope,10) == num2str(div,10)
        ipar=1;
        xline=linspace(xmin,xmax,npts);
        yline=linspace(ymin,ymax,npts);
        parlines(ipar,1:npts,1)=xline;
        parlines(ipar,1:npts,2)=yline;
        ipar=ipar+1;
        for xi=(xmin+horz_inc):horz_inc:xmax
            yi=slope.*(xmax-xi)+ymin;
            xline=linspace(xi,xmax,npts);
            yline=linspace(ymin,yi,npts);
            parlines(ipar,1:npts,1)=xline;
            parlines(ipar,1:npts,2)=yline;
            ipar=ipar+1;
        end
        for yi=(ymin+vert_inc):vert_inc:ymax
            xi=(ymax-yi+(slope.*xmin))./slope;
            xline=linspace(xmin,xi,npts);
            yline=linspace(yi,ymax,npts);
            parlines(ipar,1:npts,1)=xline;
            parlines(ipar,1:npts,2)=yline;
            ipar=ipar+1;
        end
    elseif (slope < div)
        ipar=1;
        ya=slope.*(xmax-xmin)+ymin;
        yb=slope.*(xmin-xmax)+ymax;
        for yi=ymin:vert_inc:ymax
            if (yi <= ya)
                xi=(ymin-yi+(slope.*xmax))./slope;
                xline=linspace(xi,xmax,npts);
                yline=linspace(ymin,yi,npts);
                parlines(ipar,1:npts,1)=xline;
                parlines(ipar,1:npts,2)=yline;
                ipar=ipar+1;    
            else
                yj=slope.*(xmin-xmax)+yi;
                xline=linspace(xmin,xmax,npts);
                yline=linspace(yj,yi,npts);
                parlines(ipar,1:npts,1)=xline;
                parlines(ipar,1:npts,2)=yline;
                ipar=ipar+1;
            end
        end
        for yi=yb:vert_inc:ymax
            xi=(ymax-yi+(slope.*xmin))./slope;
            xline=linspace(xmin,xi,npts);
            yline=linspace(yi,ymax,npts);
            parlines(ipar,1:npts,1)=xline;
            parlines(ipar,1:npts,2)=yline;
            ipar=ipar+1;
        end
    elseif (slope > div)
        ipar=1;
        xa=(ymax-ymin+(slope.*xmin))./slope;
        xb=(ymin-ymax+(slope.*xmax))./slope;
        for xi=xmin:horz_inc:xmax
            if (xi <= xa)
                yi=slope.*(xmin-xi)+ymax;
                xline=linspace(xmin,xi,npts);
                yline=linspace(yi,ymax,npts);
                parlines(ipar,1:npts,1)=xline;
                parlines(ipar,1:npts,2)=yline;
                ipar=ipar+1;
            else
                xj=(ymin-ymax+(slope.*xi))./slope;
                xline=linspace(xj,xi,npts);
                yline=linspace(ymin,ymax,npts);
                parlines(ipar,1:npts,1)=xline;
                parlines(ipar,1:npts,2)=yline;
                ipar=ipar+1;
            end
        end
        for xi=xb:horz_inc:xmax
            yi=slope.*(xmax-xi)+ymin;
            xline=linspace(xi,xmax,npts);
            yline=linspace(ymin,yi,npts);
            parlines(ipar,1:npts,1)=xline;
            parlines(ipar,1:npts,2)=yline;
            ipar=ipar+1;
        end
    end
elseif (slope < 0)
    % this following code is probably the hard way to draw these parallel lines.
    % the easier way would be the flipping or rotating across the appropriate
    % axis for the similar parallel lines when the slope is > 0.  I did it this
    % way because i did not want to take the time to figure out the right
    % transformation of the slope > 0 case that would give the correct
    % slope < 0 case.
    if (slope == -div)%num2str(abs(slope),10) == num2str(abs(-div),10)
        ipar=1;
        xline=linspace(xmin,xmax,npts);
        yline=linspace(ymax,ymin,npts);
        parlines(ipar,1:npts,1)=xline;
        parlines(ipar,1:npts,2)=yline;
        ipar=ipar+1;
        for xi=xmin:horz_inc:xmax
            yi=slope.*(xmin-xi)+ymin;
            xline=linspace(xmin,xi,npts);
            yline=linspace(yi,ymin,npts);
            parlines(ipar,1:npts,1)=xline;
            parlines(ipar,1:npts,2)=yline;
            ipar=ipar+1;
        end
        for yi=(ymin+vert_inc):vert_inc:ymax
            xi=(ymax-yi+(slope.*xmax))./slope;
            xline=linspace(xi,xmax,npts);
            yline=linspace(ymax,yi,npts);
            parlines(ipar,1:npts,1)=xline;
            parlines(ipar,1:npts,2)=yline;
            ipar=ipar+1;
        end
    elseif (slope > -div)
        ipar=1;
        ya=slope.*(xmin-xmax)+ymin;
        yb=slope.*(xmax-xmin)+ymax;
        for yi=ymin:vert_inc:ymax
            if (yi <= ya)
                xi=(ymin-yi+(slope.*xmin))./slope;
                xline=linspace(xmin,xi,npts);
                yline=linspace(yi,ymin,npts);
                parlines(ipar,1:npts,1)=xline;
                parlines(ipar,1:npts,2)=yline;
                ipar=ipar+1;
            else
                yj=slope.*(xmax-xmin)+yi;
                xline=linspace(xmin,xmax,npts);
                yline=linspace(yi,yj,npts);
                parlines(ipar,1:npts,1)=xline;
                parlines(ipar,1:npts,2)=yline;
                ipar=ipar+1;
            end
        end
        for yi=yb:vert_inc:ymax
            xi=(ymax-yi+(slope.*xmax))./slope;
            xline=linspace(xi,xmax,npts);
            yline=linspace(ymax,yi,npts);
            parlines(ipar,1:npts,1)=xline;
            parlines(ipar,1:npts,2)=yline;
            ipar=ipar+1;
        end
    elseif (slope < -div)
        ipar=1;
        xa=(ymin-ymax+(slope.*xmin))./slope;
        xb=(ymax-ymin+(slope.*xmax))./slope;
        for xi=xmin:horz_inc:xmax
            if (xi <= xa)
                yi=slope.*(xmin-xi)+ymin;
                xline=linspace(xmin,xi,npts);
                yline=linspace(yi,ymin,npts);
                parlines(ipar,1:npts,1)=xline;
                parlines(ipar,1:npts,2)=yline;
                ipar=ipar+1;
            else
                xj=(ymax-ymin+(slope.*xi))./slope;
                xline=linspace(xj,xi,npts);
                yline=linspace(ymax,ymin,npts);
                parlines(ipar,1:npts,1)=xline;
                parlines(ipar,1:npts,2)=yline;
                ipar=ipar+1;
            end
        end
        for xi=xb:horz_inc:xmax
            yi=slope.*(xmax-xi)+ymax;
            xline=linspace(xi,xmax,npts);
            yline=linspace(ymax,yi,npts);
            parlines(ipar,1:npts,1)=xline;
            parlines(ipar,1:npts,2)=yline;
            ipar=ipar+1;
        end
    end
end