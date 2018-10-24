function plot_Kps(op_str,p2_output)

if (strcmp(op_str,'top'))
    for i=1:length(p2_output(11).data)
        p2_output(12).data{p2_output(11).data(1,i),p2_output(11).data(2,i)}
        plot(ans(:,1),'og'); % donor Kp
        hold on;
        plot(ans(:,2),'sr'); % acceptor Kp
        hold off;
        pause
        clf reset
    end
    close
end

if (strcmp(op_str,'all'))
    [nrsets,ncsets]=size(p2_output(12).data);
    for r=1:nrsets
        row_index=r
        for c=1:ncsets
            column_index=c
            p2_output(12).data{r,c};
            plot(ans(:,1),'og'); % donor Kp
            hold on;
            plot(ans(:,2),'sr'); % acceptor Kp
            hold off;
            pause
            clf reset
        end
    end
    close
end