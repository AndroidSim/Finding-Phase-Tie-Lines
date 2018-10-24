function plot_surface(calc_surf)
if (iscell(calc_surf))
    [nrsets,ncsets]=size(calc_surf);
    if (nrsets == 1 | ncsets == 1)
        figure;
        for i=1:nrsets
            cell=i
            plot_surface(calc_surf{i});
            pause;
            clf reset;
        end
    else
        figure;
        for r=1:nrsets
            row_index=r
            for c=1:ncsets
                column_index=c
                plot_surface(calc_surf{r,c});
                pause;
                clf reset;
            end
        end
    end
end
if (ndims(calc_surf)==3)
    [nx,ny,nsurfs]=size(calc_surf);
    for i=1:nsurfs
        % c_surf=calc_surf(:,:,i);
        c_surf=calc_surf(1:nx,1:ny,i);
        surf(c_surf);
        i
        pause
    end
end
if (ndims(calc_surf)==2)
    surf(calc_surf);
end