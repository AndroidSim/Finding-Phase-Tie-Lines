function plot_tielines(tielines)
if (iscell(tielines))
    % nsets=size(tielines,2);
    [nrsets,ncsets]=size(tielines);
    if (nrsets == 1 | ncsets == 1)
        figure;
        for i=1:nrsets
            cell=i
            plot_tielines(tielines{i});
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
        % tylines=tielines(:,:,:,i);
        tylines=tielines(1:ntielines,1:ntpts,1:xy,i);
        plot_tielines(tylines);
        i
%        hold on;
%        for j=1:ntielines
%            plot(tylines(j,:,1),tylines(i,:,2),':g');
            %         end
%        hold off;
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
    hold off;
end