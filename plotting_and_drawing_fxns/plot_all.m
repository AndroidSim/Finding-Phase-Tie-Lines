function plot_all(tielines,ab_pts,all_surfaces)
if (~iscell(tielines) & ~iscell(ab_pts) & ~iscell(all_surfaces))
    scnsize=get(0,'ScreenSize');
    figure('Position',[scnsize(1)*20 scnsize(4)/2 scnsize(3)-scnsize(1)*40 scnsize(4)/2-30]);
    subplot(1,2,1);
    plot_tielines(tielines);
    hold on;
    plot_points(ab_pts);
    hold off;
    subplot(1,2,2);
    plot_surface(all_surfaces);
    %error('one or all inputs are not cell arrays')
else
    [nrsets,ncsets]=size(tielines); % dimensions should be same for ab_pts and all_surfaces
    scnsize=get(0,'ScreenSize');
    if (nrsets == 1 | ncsets == 1)
        figure('Position',[scnsize(1)*20 scnsize(4)/2 scnsize(3)-scnsize(1)*40 scnsize(4)/2-30]);
        for i=1:nrsets
            cell=i
            subplot(1,2,1);
            plot_tielines(tielines{i});
            hold on;
            plot_points(ab_pts{i});
            hold off;
            subplot(1,2,2);
            plot_surface(all_surfaces{i});
            pause;
            clf reset;
        end
    else
        figure('Position',[scnsize(1)*20 scnsize(4)/2 scnsize(3)-scnsize(1)*40 scnsize(4)/2-30]);
        for r=1:nrsets
            row_index=r
            for c=1:ncsets
                column_index=c
                subplot(1,2,1);
                plot_tielines(tielines{r,c});
                hold on;
                plot_points(ab_pts{r,c});
                hold off;
                subplot(1,2,2);
                plot_surface(all_surfaces{r,c});
                pause;
                clf reset;
            end
        end
    end
end