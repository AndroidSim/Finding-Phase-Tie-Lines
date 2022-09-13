function output = ask4input(request,boundary,tieline)

switch request
    case 'bdy_confing'
        reply = input(sprintf('"2critpts", "1critpt/1endtl", or "2endtls" boundary configuration?\n'),'s');
        while isempty(reply) || ~any(strcmp(reply,{'2critpts','1critpt/1endtl','2endtls'}))
            reply = input(sprintf('"2critpts", "1critpt/1endtl", or "2endtls" boundary configuration?\n'),'s');
        end
        output = reply;
    case 'tieline_config'
        reply = input(sprintf('"parallel", "tangent", or "ruled" tieline configuration?\n'),'s');
        while isempty(reply) || ~any(strcmp(reply,{'parallel','tangent','ruled'}))
            reply = input(sprintf('"parallel", "tangent", or "ruled" tieline configuration?\n'),'s');
        end
        output = reply;
    case 'constraints'
        switch tieline.config
            case 'parallel'
                switch boundary.config
                    case '2critpts'
                        reply = input(sprintf('is the slope constant? y/n\n'),'s');
                        while isempty(reply) || ~any(strcmp(reply,{'y','n'}))
                            reply = input(sprintf('is the slope constant? y/n\n'),'s');
                        end
                        if strcmp(reply,'y') % reply == 'y'
                            slope_position = 'located';
                            reply = input(sprintf('what is the slope in angles? [0:180]\n'));
                            while isempty(reply) || ~isscalar(reply) || (reply < 0 | reply > 180)
                                reply = input(sprintf('what is the slope in angles? [0:180]\n'));
                            end
                            m = reply.*(1/180);  
                        else % reply == 'n'
                            reply = input(sprintf('is there a range of slope angles to search? y/n\n'),'s');
                            while isempty(reply) || ~any(strcmp(reply,{'y','n'}))
                                reply = input(sprintf('is there a range of slope angles to search? y/n\n'),'s');
                            end
                            if strcmp(reply,'y') % reply == 'y'
                                slope_position = 'range';
                                reply = input(sprintf('what is the range of slope angles to search? [0:180] (if unknown, press enter)\n'));
                                while ~isempty(reply) || ~isvector(reply) || (any(reply < 0) | any(reply > 180))
                                    reply = input(sprintf('what is the range of slope angles to search? [0:180] (if unknown, press enter)\n'));
                                end
                                if isempty(reply)
                                    slope_position = 'unknown';
                                    m = [0 1];
                                end
                                if isvector(reply) & all(reply >= 0 & reply <= 180)
                                    m = reply.*(1/180);
                                end
                            else % reply == 'n'
                                slope_position = 'unknown';
                                m = [0 1];
                            end      
                        end 
                        output.reply1 = slope_position;
                        output.reply2 = m;
                    case '1critpt/1endtl'
                        % loop over end points of end tielines and get their coordinates
                        for i = 1:2
                            reply = input(sprintf('do you know the location of end point %d? y/n\n',i),'s'); 
                            if reply == 'y'
                                boundary.(['endpt' num2str(i) 'position']) = 'located';
                                reply = input(sprintf('what is the boundary parameter for end point %d? [0:1] (if unknown press enter)\n',i));
                                if isempty(reply)
                                    disp(sprintf('click the point on the boundary where end point %d is located:\n',i));
                                    ternary_plot(C_bdy,'-k','linewidth',3);
                                    [x,y] = ginput(1);
                                    pause(3);
                                    close;
                                    p = bdypt2b([x y],b,C_bdy);
                                    endpt = b2bdypt(b,C_bdy,p);
                                    boundary.(['endpt' num2str(i)]) = endpt;
                                    boundary.(['endpt' num2str(i) 'parameter']) = p;
                                elseif isscalar(reply)
                                    p = reply;
                                    endpt = b2bdypt(b,C_bdy,p);
                                    boundary.(['endpt' num2str(i) 'parameter']) = p;
                                    boundary.(['endpt' num2str(i)]) = endpt;
                                else
                                    error('invalid response');
                                end             
                            elseif reply == 'n'
                                reply = input(sprintf('is there a range of possible locations for end point %d? y/n\n',i),'s');
                                if reply == 'y'
                                    boundary.(['endpt' num2str(i) 'position']) = 'range';
                                    reply = input(sprintf('what is the boundary parameter range to search for end point %d? [0:1] (if unknown press enter)\n',i));
                                    if isempty(reply)
                                        disp(sprintf('click 2 points on the boundary to indicate the possible range:\n'));
                                        ternary_plot(C_bdy,'-k','linewidth',3);
                                        [x,y] = ginput(2);
                                        pause(3);
                                        close;
                                        p = bdypt2b([x y],b,C_bdy);
                                        endpt = b2bdypt(b,C_bdy,p);
                                        boundary.(['endpt' num2str(i)]) = endpt;
                                        boundary.(['endpt' num2str(i) 'parameter']) = p;
                                    elseif isvector(reply)
                                        p = reply;
                                        endpt = b2bdypt(b,C_bdy,p);
                                        boundary.(['endpt' num2str(i) 'parameter']) = p;
                                        boundary.(['endpt' num2str(i)]) = endpt;
                                    else
                                        error('invalid response');
                                    end    
                                elseif reply == 'n'
                                    boundary.(['endpt' num2str(i) 'position']) = 'unknown';
                                    boundary.(['endpt' num2str(i) 'parameter']) = [0 1];
                                else
                                    error('invalid response');
                                end      
                            else
                                error('invalid response');
                            end  
                        end % end of loop over end points
                        output.reply1 = ;
                        output.reply2 = ;
                        output.reply3 = ;
                        output.reply4 = ;
                    case '2endtls'
                    otherwise
                        error('invalid boundary configuration');
                end
            case 'tangent'
                switch boundary.config
                    case '2critpts'
                        % loop over number of critical points and get their
                        % coordinates.
                        for i = 1:2
                            reply = input(sprintf('do you know the location for critical point %d? y/n\n',i),'s'); 
                            if reply == 'y'
                                boundary.(['critpt' num2str(i) 'position']) = 'located';
                                reply = input(sprintf('what is the boundary parameter for critical point %d? [0:1] (if unknown press enter)\n',i));
                                if isempty(reply)
                                    disp(sprintf('click the point on the boundary where critical point %d is located:\n',i));
                                    ternary_plot(C_bdy,'-k','linewidth',3);
                                    [x,y] = ginput(1);
                                    pause(3);
                                    close;
        %                             p = bdypt2b(cart2tern([x y],1),b,C_bdy);
                                    p = bdypt2b([x y],b,C_bdy);
                                    critpt = b2bdypt(b,C_bdy,p);
                                    boundary.(['critpt' num2str(i)]) = critpt;
                                    boundary.(['critpt' num2str(i) 'parameter']) = p;
                                elseif isscalar(reply)
                                    p = reply;
                                    critpt = b2bdypt(b,C_bdy,p);
                                    boundary.(['critpt' num2str(i) 'parameter']) = p;
                                    boundary.(['critpt' num2str(i)]) = critpt;
                                else
                                    error('invalid response');
                                end             
                            elseif reply == 'n'
                                reply = input(sprintf('is there a range of possible locations for critical point %d? y/n\n',i),'s');
                                if reply == 'y'
                                    boundary.(['critpt' num2str(i) 'position']) = 'range';
                                    reply = input(sprintf('what is the boundary parameter range to search for critical point %d? [0:1] (if unknown press enter)\n',i));
                                    if isempty(reply)
                                        disp(sprintf('click 2 points on the boundary to indicate the possible range:\n'));
                                        ternary_plot(C_bdy,'-k','linewidth',3);
                                        [x,y] = ginput(2);
                                        pause(3);
                                        close;
        %                                 p = bdypt2b(cart2tern([x y],1),b,C_bdy);
                                        p = bdypt2b([x y],b,C_bdy);
                                        critpt = b2bdypt(b,C_bdy,p);
                                        boundary.(['critpt' num2str(i)]) = critpt;
                                        boundary.(['critpt' num2str(i) 'parameter']) = p;
                                    elseif isvector(reply)
                                        p = reply;
                                        critpt = b2bdypt(b,C_bdy,p);
                                        boundary.(['critpt' num2str(i) 'parameter']) = p;
                                        boundary.(['critpt' num2str(i)]) = critpt;
                                    else
                                        error('invalid response');
                                    end    
                                elseif reply == 'n'
                                    boundary.(['critpt' num2str(i) 'position']) = 'unknown';
                                    boundary.(['critpt' num2str(i) 'parameter']) = [0 1];
                                else
                                    error('invalid response');
                                end      
                            else
                                error('invalid response');
                            end  
                        end % end of loop over number critical points
                    case '1critpt/1endtl'
                        reply = input(sprintf('do you know the location of the critical point? y/n\n'),'s'); 
                        if reply == 'y'
                            boundary.('critptposition') = 'located';
                            reply = input(sprintf('what is the boundary parameter for the critical point? [0:1] (if unknown press enter)\n'));
                            if isempty(reply)
                                disp(sprintf('click the point on the boundary where the critical point is located:\n'));
                                ternary_plot(C_bdy,'-k','linewidth',3);
                                [x,y] = ginput(1);
                                pause(3);
                                close;
                                p = bdypt2b([x y],b,C_bdy);
                                critpt = b2bdypt(b,C_bdy,p);
                                boundary.('critpt') = critpt;
                                boundary.('critptparameter') = p;
                            elseif isscalar(reply)
                                p = reply;
                                critpt = b2bdypt(b,C_bdy,p);
                                boundary.('critptparameter') = p;
                                boundary.('critpt') = critpt;
                            else
                                error('invalid response');
                            end             
                        elseif reply == 'n'
                            reply = input(sprintf('is there a range of possible locations for the critical point? y/n\n'),'s');
                            if reply == 'y'
                                boundary.('critptposition') = 'range';
                                reply = input(sprintf('what is the boundary parameter range to search for the critical point? [0:1] (if unknown press enter)\n'));
                                if isempty(reply)
                                    disp(sprintf('click 2 points on the boundary to indicate the possible range:\n'));
                                    ternary_plot(C_bdy,'-k','linewidth',3);
                                    [x,y] = ginput(2);
                                    pause(3);
                                    close;
                                    p = bdypt2b([x y],b,C_bdy);
                                    critpt = b2bdypt(b,C_bdy,p);
                                    boundary.('critpt') = critpt;
                                    boundary.('critptparameter') = p;
                                elseif isvector(reply)
                                    p = reply;
                                    critpt = b2bdypt(b,C_bdy,p);
                                    boundary.('critptparameter') = p;
                                    boundary.('critpt') = critpt;
                                else
                                    error('invalid response');
                                end    
                            elseif reply == 'n'
                                boundary.('critptposition') = 'unknown';
                                boundary.('critptparameter') = [0 1];
                            else
                                error('invalid response');
                            end      
                        else
                            error('invalid response');
                        end  
                        % loop over end points of end tielines and get their coordinates
                        for i = 1:2
                            reply = input(sprintf('do you know the location of end point %d? y/n\n',i),'s'); 
                            if reply == 'y'
                                boundary.(['endpt' num2str(i) 'position']) = 'located';
                                reply = input(sprintf('what is the boundary parameter for end point %d? [0:1] (if unknown press enter)\n',i));
                                if isempty(reply)
                                    disp(sprintf('click the point on the boundary where end point %d is located:\n',i));
                                    ternary_plot(C_bdy,'-k','linewidth',3);
                                    [x,y] = ginput(1);
                                    pause(3);
                                    close;
                                    p = bdypt2b([x y],b,C_bdy);
                                    endpt = b2bdypt(b,C_bdy,p);
                                    boundary.(['endpt' num2str(i)]) = endpt;
                                    boundary.(['endpt' num2str(i) 'parameter']) = p;
                                elseif isscalar(reply)
                                    p = reply;
                                    endpt = b2bdypt(b,C_bdy,p);
                                    boundary.(['endpt' num2str(i) 'parameter']) = p;
                                    boundary.(['endpt' num2str(i)]) = endpt;
                                else
                                    error('invalid response');
                                end             
                            elseif reply == 'n'
                                reply = input(sprintf('is there a range of possible locations for end point %d? y/n\n',i),'s');
                                if reply == 'y'
                                    boundary.(['endpt' num2str(i) 'position']) = 'range';
                                    reply = input(sprintf('what is the boundary parameter range to search for end point %d? [0:1] (if unknown press enter)\n',i));
                                    if isempty(reply)
                                        disp(sprintf('click 2 points on the boundary to indicate the possible range:\n'));
                                        ternary_plot(C_bdy,'-k','linewidth',3);
                                        [x,y] = ginput(2);
                                        pause(3);
                                        close;
                                        p = bdypt2b([x y],b,C_bdy);
                                        endpt = b2bdypt(b,C_bdy,p);
                                        boundary.(['endpt' num2str(i)]) = endpt;
                                        boundary.(['endpt' num2str(i) 'parameter']) = p;
                                    elseif isvector(reply)
                                        p = reply;
                                        endpt = b2bdypt(b,C_bdy,p);
                                        boundary.(['endpt' num2str(i) 'parameter']) = p;
                                        boundary.(['endpt' num2str(i)]) = endpt;
                                    else
                                        error('invalid response');
                                    end    
                                elseif reply == 'n'
                                    boundary.(['endpt' num2str(i) 'position']) = 'unknown';
                                    boundary.(['endpt' num2str(i) 'parameter']) = [0 1];
                                else
                                    error('invalid response');
                                end      
                            else
                                error('invalid response');
                            end  
                        end % end of loop over end points
                    case '2endtls'
                    otherwise
                        error('invalid boundary configuration');
                end
            case 'ruled'
                switch boundary.config
                    case '2critpts'
                        % loop over number of critical points and get their
                        % coordinates.
                        for i = 1:2
                            reply = input(sprintf('do you know the location for critical point %d? y/n\n',i),'s'); 
                            if reply == 'y'
                                boundary.(['critpt' num2str(i) 'position']) = 'located';
                                reply = input(sprintf('what is the boundary parameter for critical point %d? [0:1] (if unknown press enter)\n',i));
                                if isempty(reply)
                                    disp(sprintf('click the point on the boundary where critical point %d is located:\n',i));
                                    ternary_plot(C_bdy,'-k','linewidth',3);
                                    [x,y] = ginput(1);
                                    pause(3);
                                    close;
        %                             p = bdypt2b(cart2tern([x y],1),b,C_bdy);
                                    p = bdypt2b([x y],b,C_bdy);
                                    critpt = b2bdypt(b,C_bdy,p);
                                    boundary.(['critpt' num2str(i)]) = critpt;
                                    boundary.(['critpt' num2str(i) 'parameter']) = p;
                                elseif isscalar(reply)
                                    p = reply;
                                    critpt = b2bdypt(b,C_bdy,p);
                                    boundary.(['critpt' num2str(i) 'parameter']) = p;
                                    boundary.(['critpt' num2str(i)]) = critpt;
                                else
                                    error('invalid response');
                                end             
                            elseif reply == 'n'
                                reply = input(sprintf('is there a range of possible locations for critical point %d? y/n\n',i),'s');
                                if reply == 'y'
                                    boundary.(['critpt' num2str(i) 'position']) = 'range';
                                    reply = input(sprintf('what is the boundary parameter range to search for critical point %d? [0:1] (if unknown press enter)\n',i));
                                    if isempty(reply)
                                        disp(sprintf('click 2 points on the boundary to indicate the possible range:\n'));
                                        ternary_plot(C_bdy,'-k','linewidth',3);
                                        [x,y] = ginput(2);
                                        pause(3);
                                        close;
        %                                 p = bdypt2b(cart2tern([x y],1),b,C_bdy);
                                        p = bdypt2b([x y],b,C_bdy);
                                        critpt = b2bdypt(b,C_bdy,p);
                                        boundary.(['critpt' num2str(i)]) = critpt;
                                        boundary.(['critpt' num2str(i) 'parameter']) = p;
                                    elseif isvector(reply)
                                        p = reply;
                                        critpt = b2bdypt(b,C_bdy,p);
                                        boundary.(['critpt' num2str(i) 'parameter']) = p;
                                        boundary.(['critpt' num2str(i)]) = critpt;
                                    else
                                        error('invalid response');
                                    end    
                                elseif reply == 'n'
                                    boundary.(['critpt' num2str(i) 'position']) = 'unknown';
                                    boundary.(['critpt' num2str(i) 'parameter']) = [0 1];
                                else
                                    error('invalid response');
                                end      
                            else
                                error('invalid response');
                            end  
                        end % end of loop over number critical points
                    case '1critpt/1endtl'
                        reply = input(sprintf('do you know the location of the critical point? y/n\n'),'s'); 
                        if reply == 'y'
                            boundary.('critptposition') = 'located';
                            reply = input(sprintf('what is the boundary parameter for the critical point? [0:1] (if unknown press enter)\n'));
                            if isempty(reply)
                                disp(sprintf('click the point on the boundary where the critical point is located:\n'));
                                ternary_plot(C_bdy,'-k','linewidth',3);
                                [x,y] = ginput(1);
                                pause(3);
                                close;
                                p = bdypt2b([x y],b,C_bdy);
                                critpt = b2bdypt(b,C_bdy,p);
                                boundary.('critpt') = critpt;
                                boundary.('critptparameter') = p;
                            elseif isscalar(reply)
                                p = reply;
                                critpt = b2bdypt(b,C_bdy,p);
                                boundary.('critptparameter') = p;
                                boundary.('critpt') = critpt;
                            else
                                error('invalid response');
                            end             
                        elseif reply == 'n'
                            reply = input(sprintf('is there a range of possible locations for the critical point? y/n\n'),'s');
                            if reply == 'y'
                                boundary.('critptposition') = 'range';
                                reply = input(sprintf('what is the boundary parameter range to search for the critical point? [0:1] (if unknown press enter)\n'));
                                if isempty(reply)
                                    disp(sprintf('click 2 points on the boundary to indicate the possible range:\n'));
                                    ternary_plot(C_bdy,'-k','linewidth',3);
                                    [x,y] = ginput(2);
                                    pause(3);
                                    close;
                                    p = bdypt2b([x y],b,C_bdy);
                                    critpt = b2bdypt(b,C_bdy,p);
                                    boundary.('critpt') = critpt;
                                    boundary.('critptparameter') = p;
                                elseif isvector(reply)
                                    p = reply;
                                    critpt = b2bdypt(b,C_bdy,p);
                                    boundary.('critptparameter') = p;
                                    boundary.('critpt') = critpt;
                                else
                                    error('invalid response');
                                end    
                            elseif reply == 'n'
                                boundary.('critptposition') = 'unknown';
                                boundary.('critptparameter') = [0 1];
                            else
                                error('invalid response');
                            end      
                        else
                            error('invalid response');
                        end  
                        % loop over end points of end tielines and get their coordinates
                        for i = 1:2
                            reply = input(sprintf('do you know the location of end point %d? y/n\n',i),'s'); 
                            if reply == 'y'
                                boundary.(['endpt' num2str(i) 'position']) = 'located';
                                reply = input(sprintf('what is the boundary parameter for end point %d? [0:1] (if unknown press enter)\n',i));
                                if isempty(reply)
                                    disp(sprintf('click the point on the boundary where end point %d is located:\n',i));
                                    ternary_plot(C_bdy,'-k','linewidth',3);
                                    [x,y] = ginput(1);
                                    pause(3);
                                    close;
                                    p = bdypt2b([x y],b,C_bdy);
                                    endpt = b2bdypt(b,C_bdy,p);
                                    boundary.(['endpt' num2str(i)]) = endpt;
                                    boundary.(['endpt' num2str(i) 'parameter']) = p;
                                elseif isscalar(reply)
                                    p = reply;
                                    endpt = b2bdypt(b,C_bdy,p);
                                    boundary.(['endpt' num2str(i) 'parameter']) = p;
                                    boundary.(['endpt' num2str(i)]) = endpt;
                                else
                                    error('invalid response');
                                end             
                            elseif reply == 'n'
                                reply = input(sprintf('is there a range of possible locations for end point %d? y/n\n',i),'s');
                                if reply == 'y'
                                    boundary.(['endpt' num2str(i) 'position']) = 'range';
                                    reply = input(sprintf('what is the boundary parameter range to search for end point %d? [0:1] (if unknown press enter)\n',i));
                                    if isempty(reply)
                                        disp(sprintf('click 2 points on the boundary to indicate the possible range:\n'));
                                        ternary_plot(C_bdy,'-k','linewidth',3);
                                        [x,y] = ginput(2);
                                        pause(3);
                                        close;
                                        p = bdypt2b([x y],b,C_bdy);
                                        endpt = b2bdypt(b,C_bdy,p);
                                        boundary.(['endpt' num2str(i)]) = endpt;
                                        boundary.(['endpt' num2str(i) 'parameter']) = p;
                                    elseif isvector(reply)
                                        p = reply;
                                        endpt = b2bdypt(b,C_bdy,p);
                                        boundary.(['endpt' num2str(i) 'parameter']) = p;
                                        boundary.(['endpt' num2str(i)]) = endpt;
                                    else
                                        error('invalid response');
                                    end    
                                elseif reply == 'n'
                                    boundary.(['endpt' num2str(i) 'position']) = 'unknown';
                                    boundary.(['endpt' num2str(i) 'parameter']) = [0 1];
                                else
                                    error('invalid response');
                                end      
                            else
                                error('invalid response');
                            end  
                        end % end of loop over end points
                    case '2endtls'
                    otherwise
                        error('invalid boundary configuration');
                end
            otherwise
                error('invalid tieline configuration');
        end       
    case 'x0'
        switch tieline.config
            case 'parallel'
                switch boundary.config
                    case '2critpts'
                        m = tieline.slope;
                        switch tieline.slope_position
                            case 'located'
                                reply = input(sprintf('enter 2 parameter starting vector as [%4.3f c]:(enter for default)\n',m*180));
                                while ~isempty(reply) || ~isvector(reply) || (reply(1) ~= m*180)
                                    reply = input(sprintf('enter 2 parameter starting vector as [%4.3f c]:(enter for default)\n',m*180));
                                end
                                if isempty(reply)
                                    x0 = [m 0]';
                                elseif isvector(reply) & length(reply) == 2 & (reply(1) == m*180) 
                                    x0 = [reply(1).*(1/180) reply(2)]';
                                else 
                                    error('invalid response');
                                end
                                output = x0;
                            case {'range','unknown'}
                                reply = input(sprintf('enter 2 parameter starting vector as [%4.3f<=m<=%4.3f c]:(enter for default)\n',[m(1)*180;m(end)*180]));
                                if isempty(reply)
                                    x0 = [m(1)+(m(end)-m(1)).*rand(1) 0]';
                                elseif isvector(reply) & length(reply) == 2 & (reply(1) >= m(1)*180 & reply(1) <= m(end)*180)
                                    x0 = [reply(1).*(1/180) reply(2)]';
                                else 
                                    error('invalid response');
                                end 
                            otherwise
                                error('invalid slope position');
                        end
                    case '1critpt/1endtl'
                        switch boundary.endpt1position
                            case 'located'
                                switch boundary.endpt2position
                                    case 'located'
                                        reply = input(sprintf('enter 4 parameter starting vector as [%4.3f %4.3f c1 c2]:(enter for default)\n',[e1;e2]));
                                        if isempty(reply)
                                            x0 = [e1 e2 0 0]';
                                        elseif isvector(reply) & length(reply) == 4 & (reply(1) == e1) & (reply(2) == e2)
                                            x0 = reply';
                                        else 
                                            error('invalid response');
                                        end 
                                    case {'range','unknown'}
                                        reply = input(sprintf('enter 4 parameter starting vector as [%4.3f %4.3f<=e2<=%4.3f c1 c2]:(enter for default)\n',[e1;e2(1);e2(end)]));
                                        if isempty(reply)
                                            x0 = [e1 e2(1)+(e2(end)-e2(1)).*rand(1) 0 0]';
                                        elseif isvector(reply) & length(reply) == 4 & (reply(1) == e1) & (reply(2) >= e2(1) & reply(2) <= e2(end))
                                            x0 = reply';
                                        else 
                                            error('invalid response');
                                        end 
                                    otherwise 
                                        error('invalid endpt 2 position');
                                end
                            case {'range','unknown'}
                                switch boundary.endpt2position
                                    case 'located'
                                        reply = input(sprintf('enter 4 parameter starting vector as [%4.3f<=e1<=%4.3f %4.3f c1 c2]:(enter for default)\n',[e1(1);e1(end);e2]));
                                        if isempty(reply)
                                            x0 = [e1(1)+(e1(end)-e1(1)).*rand(1) e2 0 0]';
                                        elseif isvector(reply) & length(reply) == 4 & (reply(1) >= e1(1) & reply(1) <= e1(end)) & (reply(2) == e2)
                                            x0 = reply';
                                        else 
                                            error('invalid response');
                                        end 
                                    case {'range','unknown'}
                                        reply = input(sprintf('enter 4 parameter starting vector as [%4.3f<=e1<=%4.3f %4.3f<=e2<=%4.3f c1 c2]:(enter for default)\n',[e1(1);e1(end);e2(1);e2(end)]));
                                        if isempty(reply)
                                            x0 = [e1(1)+(e1(end)-e1(1)).*rand(1) e2(1)+(e2(end)-e2(1)).*rand(1) 0 0]';
                                        elseif isvector(reply) & length(reply) == 4 & (reply(1) >= e1(1) & reply(1) <= e1(end)) & (reply(2) >= e2(1) & reply(2) <= e2(end))
                                            x0 = reply';
                                        else 
                                            error('invalid response');
                                        end 
                                    otherwise 
                                        error('invalid endpt 2 position');
                                end
                            otherwise 
                                error('invalid endpt 1 position');
                        end
                    case '2endtls'
                    otherwise
                        error('invalid boundary configuration');
                end
            case 'tangent'
                switch boundary.config
                    case '2critpts'
                        switch boundary.critpt1position
                            case 'located'
                                switch boundary.critpt2position
                                    case 'located'
                                        reply = input(sprintf('enter 3 parameter starting vector as [%4.3f %4.3f c]:(enter for default)\n',[c1;c2]));
                                        if isempty(reply)
                                            x0 = [c1 c2 0]';
                                        elseif isvector(reply) & length(reply) == 3 & (reply(1) == c1) & (reply(2) == c2)
                                            x0 = reply';
                                        else 
                                            error('invalid response');
                                        end 
                                    case {'range','unknown'}
                                        reply = input(sprintf('enter 3 parameter starting vector as [%4.3f %4.3f<=c2<=%4.3f c]:(enter for default)\n',[c1;c2(1);c2(end)]));
                                        if isempty(reply)
                                            x0 = [c1 c2(1)+(c2(end)-c2(1)).*rand(1) 0]';
                                        elseif isvector(reply) & length(reply) == 3 & (reply(1) == c1) & (reply(2) >= c2(1) & reply(2) <= c2(end))
                                            x0 = reply';
                                        else 
                                            error('invalid response');
                                        end 
                                    otherwise
                                        error('invalid critpt 2 position');
                                end
                            case {'range','unknown'}
                                switch boundary.critpt2position
                                    case 'located'
                                        reply = input(sprintf('enter 3 parameter starting vector as [%4.3f<=c1<=%4.3f %4.3f c]:(enter for default)\n',[c1(1);c1(end);c2]));
                                        if isempty(reply)
                                            x0 = [c1(1)+(c1(end)-c1(1)).*rand(1) c2 0]';
                                        elseif isvector(reply) & length(reply) == 3 & (reply(1) >= c1(1) & reply(1) <= c1(end)) & (reply(2) == c2)
                                            x0 = reply';
                                        else 
                                            error('invalid response');
                                        end 
                                    case {'range','unknown'}
                                        reply = input(sprintf('enter 3 parameter starting vector as [%4.3f<=c1<=%4.3f %4.3f<=c2<=%4.3f c]:(enter for default)\n',[c1(1);c1(end);c2(1);c2(end)]));
                                        if isempty(reply)
                                            x0 = [c1(1)+(c1(end)-c1(1))*rand(1) c2(1)+(c2(end)-c2(1))*rand(1) 0]';
                                        elseif isvector(reply) & length(reply) == 3 & (reply(1) >= c1(1) & reply(1) <= c1(end)) & (reply(2) >= c2(1) & reply(2) <= c2(end))
                                            x0 = reply';
                                        else 
                                            error('invalid response');
                                        end 
                                    otherwise
                                        error('invalid critpt 2 position');
                                end
                            otherwise
                                error('invalid critpt 1 position');
                        end
                    case '1critpt/1endtl'
                        switch boundary.critptposition
                            case 'located'
                                switch boundary.endpt1position
                                    case 'located'
                                        switch boundary.endpt2position
                                            case 'located'
                                                reply = input(sprintf('enter 5 parameter starting vector as [%4.3f %4.3f %4.3f c1 c2]:(enter for default)\n',[cp;e1;e2]));
                                                if isempty(reply)
                                                    x0 = [cp e1 e2 0 0]';
                                                elseif isvector(reply) & length(reply) == 5 & (reply(1) == cp) & (reply(2) == e1) & (reply(3) == e2) 
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end 
                                            case {'range','unknown'}
                                                reply = input(sprintf('enter 5 parameter starting vector as [%4.3f %4.3f %4.3f<=e2<=%4.3f c1 c2]:(enter for default)\n',[cp;e1;e2(1);e2(end)]));
                                                if isempty(reply)
                                                    x0 = [cp e1 e2(1)+(e2(end)-e2(1)).*rand(1) 0 0]';
                                                elseif isvector(reply) & length(reply) == 5 & (reply(1) == cp) & (reply(2) == e1) & (reply(3) >= e2(1) & reply(3) <= e2(end))
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end 
                                            otherwise
                                                error('invalid endpt 2 position');
                                        end
                                    case {'range','unknown'}
                                        switch boundary.endpt2position
                                            case 'located'
                                                reply = input(sprintf('enter 5 parameter starting vector as [%4.3f %4.3f<=e1<=%4.3f %4.3f c1 c2]:(enter for default)\n',[cp;e1(1);e1(end);e2]));
                                                if isempty(reply)
                                                    x0 = [cp e1(1)+(e1(end)-e1(1)).*rand(1) e2 0 0]';
                                                elseif isvector(reply) & length(reply) == 5 & (reply(1) == cp) & (reply(2) >= e1(1) & reply(2) <= e1(end)) & (reply(3) == e2) 
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end
                                            case {'range','unknown'}
                                                reply = input(sprintf('enter 5 parameter starting vector as [%4.3f %4.3f<=e1<=%4.3f %4.3f<=e2<=%4.3f c1 c2]:(enter for default)\n',[cp;e1(1);e1(end);e2(1);e2(end)]));
                                                if isempty(reply)
                                                    x0 = [cp e1(1)+(e1(end)-e1(1)).*rand(1) e2(1)+(e2(end)-e2(1)).*rand(1) 0 0]';
                                                elseif isvector(reply) & length(reply) == 5 & (reply(1) == cp) & (reply(2) >= e1(1) & reply(2) <= e1(end)) & (reply(3) >= e2(1) & reply(3) <= e2(end))
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end 
                                            otherwise
                                                error('invalid endpt 2 position');
                                        end
                                    otherwise
                                        error('invalid endpt 1position');
                                end
                            case {'range','unknown'}
                                switch boundary.endpt1position
                                    case 'located'
                                        switch boundary.endpt2position
                                            case 'located'
                                                reply = input(sprintf('enter 5 parameter starting vector as [%4.3f<=cp<=%4.3f %4.3f %4.3f c1 c2]:(enter for default)\n',[cp(1);cp(end);e1;e2]));
                                                if isempty(reply)
                                                    x0 = [cp(1)+(cp(end)-cp(1)).*rand(1) e1 e2 0 0]';
                                                elseif isvector(reply) & length(reply) == 5 & (reply(1) >= cp(1) & reply(1) <= cp(end)) & (reply(2) == e1) & (reply(3) == e2)
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end 
                                            case {'range','unknown'}
                                                reply = input(sprintf('enter 5 parameter starting vector as [%4.3f<=cp<=%4.3f %4.3f %4.3f<=e2<=%4.3f a>0 c1 c2]:(enter for default)\n',[cp(1);cp(end);e1;e2(1);e2(end)]));
                                                if isempty(reply)
                                                    x0 = [cp(1)+(cp(end)-cp(1)).*rand(1) e1 e2(1)+(e2(end)-e2(1)).*rand(1) 0 0]';
                                                elseif isvector(reply) & length(reply) == 5 & (reply(1) >= cp(1) & reply(1) <= cp(end)) & (reply(2) == e1) & (reply(3) >= e2(1) & reply(3) <= e2(end))
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end 
                                            otherwise
                                                error('invalid endpt 2 position');
                                        end
                                    case {'range','unknown'}
                                        switch boundary.endpt2position
                                            case 'located'
                                                reply = input(sprintf('enter 5 parameter starting vector as [%4.3f<=cp<=%4.3f %4.3f<=e1<=%4.3f %4.3f c1 c2]:(enter for default)\n',[cp(1);cp(end);e1(1);e1(end);e2]));
                                                if isempty(reply)
                                                    x0 = [cp(1)+(cp(end)-cp(1)).*rand(1) e1(1)+(e1(end)-e1(1)).*rand(1) e2 0 0]';
                                                elseif isvector(reply) & length(reply) == 5 & (reply(1) >= cp(1) & reply(1) <= cp(end)) & (reply(2) >= e1(1) & reply(2) <= e1(end)) & (reply(3) == e2) 
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end 
                                            case {'range','unknown'}
                                                reply = input(sprintf('enter 5 parameter starting vector as [%4.3f<=cp<=%4.3f %4.3f<=e1<=%4.3f %4.3f<=e2<=%4.3f c1 c2]:(enter for default)\n',[cp(1);cp(end);e1(1);e1(end);e2(1);e2(end)]));
                                                if isempty(reply)
                                                    x0 = [cp(1)+(cp(end)-cp(1)).*rand(1) e1(1)+(e1(end)-e1(1)).*rand(1) e2(1)+(e2(end)-e2(1)).*rand(1) 0 0]';
                                                elseif isvector(reply) & length(reply) == 5 & (reply(1) >= cp(1) & reply(1) <= cp(end)) & (reply(2) >= e1(1) & reply(2) <= e1(end)) & (reply(3) >= e2(1) & reply(3) <= e2(end)) 
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end 
                                            otherwise
                                                error('invalid endpt 2 position');
                                        end
                                    otherwise
                                        error('invalid endpt 1 position');
                                end
                            otherwise
                                error('invalid critpt position');
                        end
                    case '2endtls'
                    otherwise
                        error('invalid boundary configuration');
                end
            case 'ruled'
                switch boundary.config
                    case '2critpts'
                        switch boundary.critpt1position
                            case 'located'
                                switch boundary.critpt2position
                                    case 'located'
                                        reply = input(sprintf('enter 4 parameter starting vector as [%4.3f %4.3f a>0 c]:(enter for default)\n',[c1;c2]));
                                        if isempty(reply)
                                            x0 = [c1 c2 1 0]';
                                        elseif isvector(reply) & length(reply) == 4 & (reply(1) == c1) & (reply(2) == c2) & (reply(3) > 0)
                                            x0 = reply';
                                        else 
                                            error('invalid response');
                                        end 
                                    case {'range','unknown'}
                                        reply = input(sprintf('enter 4 parameter starting vector as [%4.3f %4.3f<=c2<=%4.3f a>0 c]:(enter for default)\n',[c1;c2(1);c2(end)]));
                                        if isempty(reply)
                                            x0 = [c1 c2(1)+(c2(end)-c2(1)).*rand(1) 1 0]';
                                        elseif isvector(reply) & length(reply) == 4 & (reply(1) == c1) & (reply(2) >= c2(1) & reply(2) <= c2(end)) & (reply(3) > 0)
                                            x0 = reply';
                                        else 
                                            error('invalid response');
                                        end 
                                    otherwise
                                        error('invalid critpt 2 position');
                                end
                            case {'range','unknown'}
                                switch boundary.critpt2position
                                    case 'located'
                                        reply = input(sprintf('enter 4 parameter starting vector as [%4.3f<=c1<=%4.3f %4.3f a>0 c]:(enter for default)\n',[c1(1);c1(end);c2]));
                                        if isempty(reply)
                                            x0 = [c1(1)+(c1(end)-c1(1)).*rand(1) c2 1 0]';
                                        elseif isvector(reply) & length(reply) == 4 & (reply(1) >= c1(1) & reply(1) <= c1(end)) & (reply(2) == c2) & (reply(3) > 0)
                                            x0 = reply';
                                        else 
                                            error('invalid response');
                                        end 
                                    case {'range','unknown'}
                                        reply = input(sprintf('enter 4 parameter starting vector as [%4.3f<=c1<=%4.3f %4.3f<=c2<=%4.3f a>0 c]:(enter for default)\n',[c1(1);c1(end);c2(1);c2(end)]));
                                        if isempty(reply)
                                            x0 = [c1(1)+(c1(end)-c1(1)).*rand(1) c2(1)+(c2(end)-c2(1)).*rand(1) 1 0]';
                                        elseif isvector(reply) & length(reply) == 4 & (reply(1) >= c1(1) & reply(1) <= c1(end)) & (reply(2) >= c2(1) & reply(2) <= c2(end)) & (reply(3) > 0)
                                            x0 = reply';
                                        else 
                                            error('invalid response');
                                        end 
                                    otherwise
                                        error('invalid critpt 2 position');
                                end
                            otherwise
                                error('invalid critpt 1 position');
                        end
                    case '1critpt/1endtl'
                        switch boundary.critptposition
                            case 'located'
                                switch boundary.endpt1position
                                    case 'located'
                                        switch boundary.endpt2position
                                            case 'located'
                                                reply = input(sprintf('enter 6 parameter starting vector as [%4.3f %4.3f %4.3f a>0 c1 c2]:(enter for default)\n',[cp;e1;e2]));
                                                if isempty(reply)
                                                    x0 = [cp e1 e2 1 0 0]';
                                                elseif isvector(reply) & length(reply) == 6 & (reply(1) == cp) & (reply(2) == e1) & (reply(3) == e2) & (reply(4) > 0)
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end 
                                            case {'range','unknown'}
                                                reply = input(sprintf('enter 6 parameter starting vector as [%4.3f %4.3f %4.3f<=e2<=%4.3f a>0 c1 c2]:(enter for default)\n',[cp;e1;e2(1);e2(end)]));
                                                if isempty(reply)
                                                    x0 = [cp e1 e2(1)+(e2(end)-e2(1)).*rand(1) 1 0 0]';
                                                elseif isvector(reply) & length(reply) == 6 & (reply(1) == cp) & (reply(2) == e1) & (reply(3) >= e2(1) & reply(3) <= e2(end)) & (reply(4) > 0)
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end
                                            otherwise
                                                error('invalid endpt 2 position');
                                        end
                                    case {'range','unknown'}
                                        switch boundary.endpt2position
                                            case 'located'
                                                reply = input(sprintf('enter 6 parameter starting vector as [%4.3f %4.3f<=e1<=%4.3f %4.3f a>0 c1 c2]:(enter for default)\n',[cp;e1(1);e1(end);e2]));
                                                if isempty(reply)
                                                    x0 = [cp e1(1)+(e1(end)-e1(1)).*rand(1) e2 1 0 0]';
                                                elseif isvector(reply) & length(reply) == 6 & (reply(1) == cp) & (reply(2) >= e1(1) & reply(2) <= e1(end)) & (reply(3) == e2) & (reply(4) > 0)
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end 
                                            case {'range','unknown'}
                                                reply = input(sprintf('enter 6 parameter starting vector as [%4.3f %4.3f<=e1<=%4.3f %4.3f<=e2<=%4.3f a>0 c1 c2]:(enter for default)\n',[cp;e1(1);e1(end);e2(1);e2(end)]));
                                                if isempty(reply)
                                                    x0 = [cp e1(1)+(e1(end)-e1(1)).*rand(1) e2(1)+(e2(end)-e2(1)).*rand(1) 1 0 0]';
                                                elseif isvector(reply) & length(reply) == 6 & (reply(1) == cp) & (reply(2) >= e1(1) & reply(2) <= e1(end)) & (reply(3) >= e2(1) & reply(3) <= e2(end)) & (reply(4) > 0)
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end 
                                            otherwise
                                                error('invalid endpt 2 position');
                                        end
                                    otherwise
                                        error('invalid endpt 1position');
                                end
                            case {'range','unknown'}
                                switch boundary.endpt1position
                                    case 'located'
                                        switch boundary.endpt2position
                                            case 'located'
                                                reply = input(sprintf('enter 6 parameter starting vector as [%4.3f<=cp<=%4.3f %4.3f %4.3f a>0 c1 c2]:(enter for default)\n',[cp(1);cp(end);e1;e2]));
                                                if isempty(reply)
                                                    x0 = [cp(1)+(cp(end)-cp(1)).*rand(1) e1 e2 1 0 0]';
                                                elseif isvector(reply) & length(reply) == 6 & (reply(1) >= cp(1) & reply(1) <= cp(end)) & (reply(2) == e1) & (reply(3) == e2) & (reply(4) > 0)
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end 
                                            case {'range','unknown'}
                                                reply = input(sprintf('enter 6 parameter starting vector as [%4.3f<=cp<=%4.3f %4.3f %4.3f<=e2<=%4.3f a>0 c1 c2]:(enter for default)\n',[cp(1);cp(end);e1;e2(1);e2(end)]));
                                                if isempty(reply)
                                                    x0 = [cp(1)+(cp(end)-cp(1)).*rand(1) e1 e2(1)+(e2(end)-e2(1)).*rand(1) 1 0 0]';
                                                elseif isvector(reply) & length(reply) == 6 & (reply(1) >= cp(1) & reply(1) <= cp(end)) & (reply(2) == e1) & (reply(3) >= e2(1) & reply(3) <= e2(end)) & (reply(4) > 0)
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end 
                                            otherwise
                                                error('invalid endpt 2 position');
                                        end
                                    case {'range','unknown'}
                                        switch boundary.endpt2position
                                            case 'located'
                                                reply = input(sprintf('enter 6 parameter starting vector as [%4.3f<=cp<=%4.3f %4.3f<=e1<=%4.3f %4.3f a>0 c1 c2]:(enter for default)\n',[cp(1);cp(end);e1(1);e1(end);e2]));
                                                if isempty(reply)
                                                    x0 = [cp(1)+(cp(end)-cp(1)).*rand(1) e1(1)+(e1(end)-e1(1)).*rand(1) e2 1 0 0]';
                                                elseif isvector(reply) & length(reply) == 6 & (reply(1) >= cp(1) & reply(1) <= cp(end)) & (reply(2) >= e1(1) & reply(2) <= e1(end)) & (reply(3) == e2) & (reply(4) > 0)
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end 
                                            case {'range','unknown'}
                                                reply = input(sprintf('enter 6 parameter starting vector as [%4.3f<=cp<=%4.3f %4.3f<=e1<=%4.3f %4.3f<=e2<=%4.3f a>0 c1 c2]:(enter for default)\n',[cp(1);cp(end);e1(1);e1(end);e2(1);e2(end)]));
                                                if isempty(reply)
                                                    x0 = [cp(1)+(cp(end)-cp(1)).*rand(1) e1(1)+(e1(end)-e1(1)).*rand(1) e2(1)+(e2(end)-e2(1)).*rand(1) 1 0 0]';
                                                elseif isvector(reply) & length(reply) == 6 & (reply(1) >= cp(1) & reply(1) <= cp(end)) & (reply(2) >= e1(1) & reply(2) <= e1(end)) & (reply(3) >= e2(1) & reply(3) <= e2(end)) & (reply(4) > 0)
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end 
                                            otherwise
                                                error('invalid endpt 2 position');
                                        end
                                    otherwise
                                        error('invalid endpt 1 position');
                                end
                            otherwise
                                error('invalid critpt position');
                        end
                    case '2endtls'
                    otherwise
                        error('invalid boundary configuration');
                end
            otherwise 
                error('invalid tieline configuration');
        end
    otherwise
        error('invalid request for input');
end

return