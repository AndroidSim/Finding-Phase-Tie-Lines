function plot_tielinesandpoints(tielines,ab_pts)

if (iscell(tielines))
    [nrsets,ncsets]=size(tielines);
    if (nrsets == 1 | ncsets == 1)
        figure;
        for i=1:nrsets
            cell=i
            plot_tielines(tielines{i});
            hold on;
            plot_points(ab_pts{i});
            hold off;
            pause;
            clf reset;
        end
    else
        figure;
        for r=1:nrsets
            row_index=r
            for c=1:ncsets
                column_index=c
                plot_tielines(tielines{r,c});
                hold on;
                plot_points(ab_pts{r,c});
                hold off;
                pause;
                clf reset;
            end
        end
    end
end
if (ndims(tielines) == 4)
    [ntielines,ntpts,xy,nbpts]=size(tielines);
    figure;
    for i=1:nbpts
        tylines=tielines(:,:,:,i);
        plot_tielines(tylines);
        i
        pause;
        clf reset;
    end
end
if (ndims(tielines) == 3)
    [ntielines,ntpts,xy]=size(tielines);
    hold on;
    for i=1:ntielines
	    plot(tielines(i,:,1),tielines(i,:,2),':g');
    end
    plot_points(ab_pts);
    hold off;
end