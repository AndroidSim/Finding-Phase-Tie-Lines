function stop = outfun(x,optimValues,state,varargin)
stop=[];
switch state
    case 'init'
        switch length(x)
            case 2
                disp(sprintf('initial point = [%f %f]',x));
            case 3
                disp(sprintf('initial point = [%f %f %f]',x));
            case 4
                disp(sprintf('initial point = [%f %f %f %f]',x));
            case 5
                disp(sprintf('initial point = [%f %f %f %f %f]',x));
            case 6
                disp(sprintf('initial point = [%f %f %f %f %f %f]',x));
            otherwise
                error('too many or too little number of parameters');
        end
    case 'iter'
        if ~isempty(varargin)
            fit_method = varargin{8};
            if strcmp(fit_method,'lsqnonlin')
                switch length(x)
                    case 2
                        disp(sprintf('current point = [%f %f]',x));
                    case 3
                        disp(sprintf('current point = [%f %f %f]',x));
                    case 4
                        disp(sprintf('current point = [%f %f %f %f]',x));
                    case 5
                        disp(sprintf('current point = [%f %f %f %f %f]',x));
                    case 6
                        disp(sprintf('current point = [%f %f %f %f %f %f]',x));
                    otherwise
                        error('too many or too little number of parameters');
                end
            end
        end
    case 'done'
%         switch length(x)
%             case 2
%                 disp(sprintf('last point = [%f %f]',x));
%             case 3
%                 disp(sprintf('last point = [%f %f %f]',x));
%             case 4
%                 disp(sprintf('last point = [%f %f %f %f]',x));
%             case 5
%                 disp(sprintf('last point = [%f %f %f %f %f]',x));
%             case 6
%                 disp(sprintf('last point = [%f %f %f %f %f %f]',x));
%             otherwise
%                 error('too many or too little number of parameters');
%         end
otherwise
end

return