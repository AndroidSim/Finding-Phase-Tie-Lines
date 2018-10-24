function yesorno = is_above_line(qpt,line_pt1,line_pt2)

% is_above_line determines whether or not a point of interest (qpt) is above or below a line segment joining
% line_pt1 and line_pt2.  in case line segment is vertical (slope = m = Inf) above is to the right and below
% is to the left

% get equation of line segment
warning off
m=(line_pt1(2)-line_pt2(2))/(line_pt1(1)-line_pt2(1));
warning on
% equation == y = m*(x-line_pt1(1))+line_pt1(2) if m ~= Inf
%             x = line_pt1(1)                   if m == Inf
if (isinf(m))
    %keyboard
    if (qpt(1) > line_pt1(1)) % or line_pt2
        yesorno=logical(1);
    else
        yesorno=logical(0);
    end
elseif (iseqwithn(abs(m),0))
    if (qpt(2) > line_pt1(2)) % or line_pt2
        yesorno=logical(1);
    else
        yesorno=logical(0);
    end
else
    if (qpt(2) > m*(qpt(1)-line_pt1(1))+line_pt1(2))
        yesorno=logical(1);
    else
        yesorno=logical(0);
    end
end     