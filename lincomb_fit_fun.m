function [V,D] = lincomb_fit_fun(kp,f_beta,varargin)

switch varargin{1}
    case 'spectrum'
        switch varargin{2}
            case 'lsqcurvefit'
                B = varargin{3};
                B_alpha = B(:,1);
                B_beta = B(:,2);
                V = ((1-f_beta)/((1-f_beta) + kp*f_beta)).*B_alpha + ((kp*f_beta)/((1-f_beta) + kp*f_beta)).*B_beta;
                
                if nargout > 1
                    % calculate gradient D
                    fb = f_beta;
                    Sa = B(:,1);
                    Sb = B(:,2);
                    D = (fb*(fb-1)/((1-fb*(1-kp))^2)).*Sa + (fb*(1-fb)/((1-fb*(1-kp))^2)).*Sb;
                end
            case 'chisq_min'
                S = varargin{3};
                B = varargin{4};
                nb = varargin{5};
                std_bline = varargin{6};
                B_alpha = B(:,1);
                B_beta = B(:,2);
                T = ((1-f_beta)/((1-f_beta) + kp*f_beta)).*B_alpha + ((kp*f_beta)/((1-f_beta) + kp*f_beta)).*B_beta;
                V = (1/((nb-1)*(std_bline^2)))*((S-T)'*(S-T));
            otherwise
                error('invalid spectrum fitting algorithm in fit function');
        end
    case 'fpafpb'
        V = (kp.*f_beta)./((1-f_beta) + kp*f_beta); % == fraction probe in beta
    case 'tieline'
        B = varargin{2};
        B_alpha = B(:,1);
        B_beta = B(:,2);
        nb = length(B_alpha);% or B_beta
        ns = length(f_beta);
        for s = 1:ns
            v(:,s) = ((1-f_beta(s))/((1-f_beta(s)) + kp*f_beta(s))).*B_alpha + ((kp*f_beta(s))/((1-f_beta(s)) + kp*f_beta(s))).*B_beta;
        end
        V = reshape(v,nb*ns,1);
        if nargout > 1
            % calculate gradient D
            fb = f_beta;
            Sa = B(:,1);
            Sb = B(:,2);
            D = (fb*(fb-1)/((1-fb*(1-kp))^2)).*Sa + (fb*(1-fb)/((1-fb*(1-kp))^2)).*Sb;
        end
    otherwise
        error('invalid fitting algorithm for data type in fit function');
end

return