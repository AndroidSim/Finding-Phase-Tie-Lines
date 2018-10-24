function F = FRET_theory_calc(x,xdata,exp_surface,bd_pts,bd_rep,Kp_fxn_type,Xd,Xa,Ca,Cb,varargin)

% FRET_theory_calc calculates FRET based on Jeff Buboltz's model for the data points in xdata and is built to
% be used with Matlab's lsqcurvefit fxn

% x vector contains model parameters (or the gammas): x(i) = gamma(i)
if (length(x) == 2 & length(varargin) == 3)
    if (varargin{1} == 1)
        g1=varargin{2};
        g2=varargin{3};
        for i=1:size(xdata,1)
            %disp('finding tie line through data point: alpha and xt');
            [alpha,xt,alpha_pt,beta_pt]=tieline_pt_coords(bd_pts,bd_rep,g1,g2,xdata(i,1),xdata(i,2));
            if (~isfinite(alpha) | ~isfinite(xt) | ~isfinite(alpha_pt) | ~isfinite(beta_pt))
                error('one of the outputs of tieline_pt_coords is not finite')
            end
            %disp('evaluating model fret at data point');
            F(i)=calc_FRET_at_pt(exp_surface,alpha,xt,Kp_fxn_type,x(1),x(2),alpha_pt,beta_pt,Xd,Xa,Ca,Cb);
            if (~isfinite(F(i)))
                error('the calculated FRET value is not finite')
            end
        end
        F=F';
        return
    elseif (varargin{1} == 2)
        g3=varargin{2};
        g4=varargin{3};
        for i=1:size(xdata,1)
            %disp('finding tie line through data point: alpha and xt');
            [alpha,xt,alpha_pt,beta_pt]=tieline_pt_coords(bd_pts,bd_rep,x(1),x(2),xdata(i,1),xdata(i,2));
            if (~isfinite(alpha) | ~isfinite(xt) | ~isfinite(alpha_pt) | ~isfinite(beta_pt))
                error('one of the outputs of tieline_pt_coords is not finite')
            end
            %disp('evaluating model fret at data point');
            F(i)=calc_FRET_at_pt(exp_surface,alpha,xt,Kp_fxn_type,g3,g4,alpha_pt,beta_pt,Xd,Xa,Ca,Cb);
            if (~isfinite(F(i)))
                error('the calculated FRET value is not finite')
            end
        end
        F=F';
        return
    end
else
    for i=1:size(xdata,1)
        %disp('finding tie line through data point: alpha and xt');
        [alpha,xt,alpha_pt,beta_pt]=tieline_pt_coords(bd_pts,bd_rep,x(1),x(2),xdata(i,1),xdata(i,2));
        if (~isfinite(alpha) | ~isfinite(xt) | ~isfinite(alpha_pt) | ~isfinite(beta_pt))
            error('one of the outputs of tieline_pt_coords is not finite')
        end
        %disp('evaluating model fret at data point');
        F(i)=calc_FRET_at_pt(exp_surface,alpha,xt,Kp_fxn_type,x(3),x(4),alpha_pt,beta_pt,Xd,Xa,Ca,Cb);
        if (~isfinite(F(i)))
            error('the calculated FRET value is not finite')
        end
    end
    F=F';
    return
end