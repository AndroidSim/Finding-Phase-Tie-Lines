for i=1:length(p2_output(11).data)
    figure('Position',[scnsize(1)*20 scnsize(4)/2 scnsize(3)-scnsize(1)*40 scnsize(4)/2-30]);
    subplot(1,2,1);
    plot_tielines(p2_output(3).data{p2_output(11).data(1,i),p2_output(11).data(2,i)});
    hold on;
    plot_points(p2_output(4).data{p2_output(11).data(1,i),p2_output(11).data(2,i)});
    hold off;
    subplot(1,2,2);
    if (p2_output(7).data(p2_output(11).data(1,i),p2_output(11).data(2,i))=='a')
        plot_surface(p2_output(1).data{p2_output(11).data(1,i),p2_output(11).data(2,i)});
    elseif (p2_output(7).data(p2_output(11).data(1,i),p2_output(11).data(2,i))=='b')
        plot_surface(p2_output(2).data{p2_output(11).data(1,i),p2_output(11).data(2,i)});
    else
        plot_surface(p2_output(1).data{p2_output(11).data(1,i),p2_output(11).data(2,i)});
    end
    fit_parameter=p2_output(9).data(ans(1,i),ans(2,i))
    pause
    clf reset
end