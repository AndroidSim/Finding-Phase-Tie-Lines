function [line_pts]=vec2line(vector,point,xmin,xmax,ymin,ymax,npts)
% the function vec2line takes a vector from one point to another and calculates the line from
% xmin to xmax or ymin to ymax that the vector corresponds to and passing through the point point
slope=vector(2)./vector(1);
if (isnan(slope))
    error('vector = [0 0], in function vec2line')
elseif (isinf(slope))
    yline=linspace(ymin,ymax,npts);
    xline=ones(1,npts).*point(1);
    %xline=repmat(point(1),[1 npts])
    line_pts=cat(2,xline',yline');
    %line_pts=[xline',yline'];
elseif (slope == 0)
    xline=linspace(xmin,xmax,npts);
    yline=slope.*(xline-point(1))+point(2);
    line_pts=cat(2,xline',yline');
    %line_pts=[xline',yline'];
else
    xline=linspace(xmin,xmax,npts);
    yline=slope.*(xline-point(1))+point(2);
    line_pts=cat(2,xline',yline');
    %line_pts=[xline',yline'];
end