function plot_points(pts)
if (iscell(pts))
    [nrsets,ncsets]=size(pts);
    if (nrsets == 1 | ncsets == 1)
        figure;
        for i=1:nsets
            cell=i
            plot_points(pts{i});
            pause;
            clf reset;
        end
    else
        figure;
        for r=1:nrsets
            row_index=r
            for c=1:ncsets
                column_index=c
                plot_points(pts{r,c});
                pause;
                clf reset;
            end
        end
    end
end
if (ndims(pts) == 3)
    nsets=size(pts,3);
    if (nsets == 2)
        hold on;
        plot(pts(:,1,1),pts(:,2,1),'or');
        plot(pts(:,1,2),pts(:,2,2),'sb');
        hold off;
    end
    if (nsets == 3)
        hold on;
        plot(pts(:,1,1),pts(:,2,1),'or');
        plot(pts(:,1,2),pts(:,2,2),'sb');
        plot(pts(:,1,3),pts(:,2,3),'dm');
        hold off;
    end
end
if (ndims(pts) == 2)
    hold on;
    plot(pts(:,1),pts(:,2),'-k.');
    hold off;
end