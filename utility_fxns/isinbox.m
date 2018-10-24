function yesorno = isinbox(pt1,pt2,qpt)

% isinbox(pt1,pt2) determine whether a query pt (qpt) lie within the box formed by the x and y coords of pt1
% and pt2 (length == deltax and width == deltay)

% all points of the form [x y]
xwn=0;
ywn=0;
if (pt1 == pt2)
    if (qpt == pt1) % or pt2
        yesorno=logical(1);
    end
else
    if (pt1(1) < pt2(1))
        if (pt1(1) < qpt(1) & qpt(1) < pt2(1))
            xwn=1;
        end
        if (pt1(2) < pt2(2))
            if (pt1(2) < qpt(2) & qpt(2) < pt2(2))
                ywn=1;
            end
        else
            if (pt2(2) < qpt(2) & qpt(2) < pt1(2))
                ywn=1;
            end
        end
    else
        if (pt2(1) < qpt(1) & qpt(1) < pt1(1))
            xwn=1;
        end
        if (pt1(2) < pt2(2))
            if (pt1(2) < qpt(2) & qpt(2) < pt2(2))
                ywn=1;
            end
        else
            if (pt2(2) < qpt(2) & qpt(2) < pt1(2))
                ywn=1;
            end
        end
    end
end
if (xwn == 1 & ywn == 1)
    yesorno=logical(1);
else
    yesorno=logical(0);
end