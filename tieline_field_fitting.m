function output = tieline_field_fitting(S_coex,C_coex,S_bdy,C_bdy,fit_method,varargin)
% output = global_tieline_fitting(S_coex,C_coex,S_bdy,C_bdy,varargin)
% global_tieline_fitting is a program to fit ESR spectra within a
% coexistence region at evenly or randomly distributed compositions using
% ESR spectra obtained on the coexistence phase boundary by searching
% different tieline configurations.
%
% related programs:  tie_line_search.m, ref_tie_line_search.m,
% tie_line_m_kp_fit_fxn.m. ref_tie_line_fit_fxn.m,
% ref_tie_line_pt_m_fit_fxn.m, etc ...
% 
% input arguments:
% S_coex = ESR spectra within the phase coexistence region.
% C_coex = the compositions of the spectra in S_coex.
% S_bdy = ESR spectra on the phase boundary.
% C_bdy = the composition of the spectra in S_bdy.
% np = number of coexisting phases.
% nc = number of chemically distinct components (not including the spin
% probe)
%
% Bdy = a different list of compositions on the boundary other than C_bdy.
% (points in C_bdy must lie on the boundary represented by 1st variable
% argument)
%
% variable input arguments:

disp(sprintf('\nchecking input arguments...\n'));

% checking coexistence spectra
if all(size(S_coex) > 1) & ndims(S_coex) == 2 % if spectra is a matrix
    [ncB,ncScol] = size(S_coex);

    if ncScol < 2 || rem(ncScol,2) ~= 0
        error('each coexistence region spectrum is two columns: [B-field intensity_values]');
    end
    
    ncS = ncScol/2;
else
    error('first argument must be a matrix of the coexistence region spectra');
end

% checking coexistence region compositions
if  any(size(C_coex) == 1) & any(size(C_coex) > 1) % if compositions is a vector
    [nrow,ncol] = size(C_coex);
    
    if ncol == 1 % if column vector
        C_coex = C_coex'; % convert to row vector
    end
    
    ncC = 1; 
    cCd = length(C_coex);
    clear nrow ncol;
elseif all(size(C_coex) > 1) & ndims(C_coex) == 2 % if compositions is a matrix 
    [ncC,cCd] = size(C_coex);
end

if ncC ~= ncS
    error('number of coexistence region compositions must equal number of coexistence region spectra');
end

% sum of coexistence region mole fractions must equal one
if all(sum(C_coex')) ~= 1
    error('sum of coexistence region mole fractions must equal one');
end

% checking boundary spectra
if all(size(S_bdy) > 1) & ndims(S_bdy) == 2 % if spectra is a matrix
    [nbB,nbScol] = size(S_bdy);

    if nbScol < 2 || rem(nbScol,2) ~= 0
        error('each boundary spectrum is two columns: [B-field intensity_values]');
    end
    
    nbS = nbScol/2;
else
    error('third argument must be a matrix of the boundary spectra');
end

% checking boundary compositions
if  any(size(C_bdy) == 1) & any(size(C_bdy) > 1) % if compositions is a vector
    [nrow,ncol] = size(C_bdy);
    
    if ncol == 1 % if column vector
        C_bdy = C_bdy'; % convert to row vector
    end
    
    nbC = 1; 
    bCd = length(C_bdy);
    clear nrow ncol; 
elseif all(size(C_bdy) > 1) & ndims(C_bdy) == 2 % if compositions is a matrix 
    [nbC,bCd] = size(C_bdy);
end

if nbC ~= nbS
    error('number of boundary compositions must equal number of boundary spectra');
end
    
% sum of boundary mole fractions must equal one
if all(sum(C_bdy')) ~= 1
    error('sum of boundary mole fractions with spectra must equal one');
end

% if ~isscalar(nc)
%     error('number of chemical components must be a scalar');
% end
% 
% if ~isscalar(np)
%     error('number of phases must be a scalar');
% end

if isempty(fit_method)
    fit_method = 'simplex';
elseif ischar(fit_method)
    if ~any(strcmp(fit_method,{'lsqnonlin';'levmar';'simplex'}))
        error('fit method must be either "lsqnonlin","levmar",or "simplex"');
    end 
else
    error('fit method must be either "lsqnonlin","levmar",or "simplex"');
end

% declare and initialize boundary and tieline data structures
boundary = struct('config',[]);
tieline = struct('config',[]);
part_coeff.config = 'function';

if isempty(varargin)
    interactive = true;
    % ask user for the phase boundary configuration
    reply = ask4input('bdy_config',boundary,tieline);
    boundary.config = reply;
    % ask user for the tieline configuration
    reply = ask4input('tieline_config',boundary,tieline);
    tieline.config = reply;
elseif isstruct(varargin{1})
    interactive = false;
    fit_conditions = varargin{1};
    if isfield(fit_conditions,'bdy_config')
        boundary.config = fit_conditions.bdy_config;
    else
        error('fit_conditions structure must have a field "bdy_config" containing the boundary configuration');
    end
    if isfield(fit_conditions,'tieline_config')
        tieline.config = fit_conditions.tieline_config;
    else
        error('fit_conditions structure must have a field "tieline_config" containing the tieline configuration');
    end
    switch tieline.config
        case 'parallel'
            switch boundary.config
                case '2critpts'
                    if isfield(fit_conditions,'slope')
                        if isscalar(fit_conditions.slope)
                            tieline.slope_position = 'located';
                            tieline.slope = fit_conditions.slope;
                        elseif isvector(fit_conditions.slope)
                            tieline.slope_position = 'range';
                            tieline.slope = fit_conditions.slope;
                        elseif isempty(fit_conditions.slope)
                            tieline.slope_position = 'unknown';
                            tieline.slope = [0 1];
                        else
                            error('slope must be a scalar, a vector specifying a range, or empty');
                        end
                    else
                        error('fit_conditions structure must have a field "slope" containing the parallel tieline slope');
                    end
                case '1critpt/1endtl'
                    if isfield(fit_conditions,'endpt1parameter')
                        if isscalar(fit_conditions.endpt1parameter)
                            boundary.endpt1position = 'located';
                            boundary.endpt1parameter = fit_conditions.endpt1parameter;
                        elseif isvector(fit_conditions.endpt1parameter)
                            boundary.endpt1position = 'range';
                            boundary.endpt1parameter = fit_conditions.endpt1parameter;
                        elseif isempty(fit_conditions.endpt1parameter)
                            boundary.endpt1position = 'range';
                            boundary.endpt1parameter = [0 1];
                        else
                            error('endpt 1 parameter must be a scalar, a vector specifying a range, or empty');
                        end
                    else
                        error('fit_conditions structure must have a field "endpt1parameter" containing the boundary parameter for endpt 1');
                    end
                    if isfield(fit_conditions,'endpt2parameter')
                        if isscalar(fit_conditions.endpt2parameter)
                            boundary.endpt2position = 'located';
                            boundary.endpt2parameter = fit_conditions.endpt2parameter;
                        elseif isvector(fit_conditions.endpt2parameter)
                            boundary.endpt2position = 'range';
                            boundary.endpt2parameter = fit_conditions.endpt2parameter;
                        elseif isempty(fit_conditions.endpt2parameter)
                            boundary.endpt2position = 'range';
                            boundary.endpt2parameter = [0 1];
                        else
                            error('endpt 2 parameter must be a scalar, a vector specifying a range, or empty');
                        end
                    else
                        error('fit_conditions structure must have a field "endpt2parameter" containing the boundary parameter for endpt 2');
                    end  
                case '2endtls'
                otherwise
                    error('invalid boundary configuration');
            end
        case {'tangent','ruled'}
            switch boundary.config
                case '2critpts'
                    if isfield(fit_conditions,'critpt1parameter')
                        if isscalar(fit_conditions.critpt1parameter)
                            boundary.critpt1position = 'located';
                            boundary.critpt1parameter = fit_conditions.critpt1parameter;
                        elseif isvector(fit_conditions.critpt1parameter)
                            boundary.critpt1position = 'range';
                            boundary.critpt1parameter = fit_conditions.critpt1parameter;
                        elseif isempty(fit_conditions.critpt1parameter)
                            boundary.critpt1position = 'range';
                            boundary.critpt1parameter = [0 1];
                        else
                            error('critpt 1 parameter must be a scalar, a vector specifying a range, or empty');
                        end
                    else
                        error('fit_conditions structure must have a field "critpt1parameter" containing the boundary parameter for critpt 1');
                    end
                    if isfield(fit_conditions,'critpt2parameter')
                        if isscalar(fit_conditions.critpt2parameter)
                            boundary.critpt2position = 'located';
                            boundary.critpt2parameter = fit_conditions.critpt2parameter;
                        elseif isvector(fit_conditions.critpt2parameter)
                            boundary.critpt2position = 'range';
                            boundary.critpt2parameter = fit_conditions.critpt2parameter;
                        elseif isempty(fit_conditions.critpt2parameter)
                            boundary.critpt2position = 'range';
                            boundary.critpt2parameter = [0 1];
                        else
                            error('critpt 2 parameter must be a scalar, a vector specifying a range, or empty');
                        end
                    else
                        error('fit_conditions structure must have a field "critpt2parameter" containing the boundary parameter for critpt 2');
                    end
                case '1critpt/1endtl'
                    if isfield(fit_conditions,'critptparameter')
                        if isscalar(fit_conditions.critptparameter)
                            boundary.critptposition = 'located';
                            boundary.critptparameter = fit_conditions.critptparameter;
                        elseif isvector(fit_conditions.critptparameter)
                            boundary.critptposition = 'range';
                            boundary.critptparameter = fit_conditions.critptparameter;
                        elseif isempty(fit_conditions.critptparameter)
                            boundary.critptposition = 'range';
                            boundary.critptparameter = [0 1];
                        else
                            error('critpt parameter must be a scalar, a vector specifying a range, or empty');
                        end
                    else
                        error('fit_conditions structure must have a field "bdy_config" containing the boundary parameter for the critpt');
                    end
                    if isfield(fit_conditions,'endpt1parameter')
                        if isscalar(fit_conditions.endpt1parameter)
                            boundary.endpt1position = 'located';
                            boundary.endpt1parameter = fit_conditions.endpt1parameter;
                        elseif isvector(fit_conditions.endpt1parameter)
                            boundary.endpt1position = 'range';
                            boundary.endpt1parameter = fit_conditions.endpt1parameter;
                        elseif isempty(fit_conditions.endpt1parameter)
                            boundary.endpt1position = 'range';
                            boundary.endpt1parameter = [0 1];
                        else
                            error('endpt 1 parameter must be a scalar, a vector specifying a range, or empty');
                        end
                    else
                        error('fit_conditions structure must have a field "tieline_config" containing the boundary parameter for endpt 1');
                    end
                    if isfield(fit_conditions,'endpt2parameter')
                        if isscalar(fit_conditions.endpt2parameter)
                            boundary.endpt2position = 'located';
                            boundary.endpt2parameter = fit_conditions.endpt2parameter;
                        elseif isvector(fit_conditions.endpt2parameter)
                            boundary.endpt2position = 'range';
                            boundary.endpt2parameter = fit_conditions.endpt2parameter;
                        elseif isempty(fit_conditions.endpt2parameter)
                            boundary.endpt2position = 'range';
                            boundary.endpt2parameter = [0 1];
                        else
                            error('endpt 2 parameter must be a scalar, a vector specifying a range, or empty');
                        end
                    else
                        error('fit_conditions structure must have a field "tieline_config" containing the boundary parameter for endpt 2');
                    end    
                case '2endtls'
                otherwise
                    error('invalid boundary configuration');
            end
        otherwise
            error('invalid tieline configuration');
    end
else
    error('varargin must be either empty or a structure specifying the non-interactive conditions for the fit');
end

% processing spectra (potential area of improvement)
disp(sprintf('processing spectra...\n'));
% disp(sprintf('shifting...\n'));
% allS = shift_spectra([S_coex S_bdy]); % this fxn can change results!!
% disp(sprintf('scaling...\n'));
% allS = scale_spectra(allS,500);
allS = [S_coex S_bdy];
disp(sprintf('aligning...\n'));
allS = align_spectra(allS); 
% normalize by converting to absorbance spectrum, transform so that the
% area of the absorbance spectrum equals one thousand, then convert back to the
% derivative spectrum.
disp(sprintf('normalizing...\n'));
allS = normalize_spectra(allS,'d','area',1000);
S_coex = allS(:,1:ncScol);
S_bdy = allS(:,ncScol+1:ncScol+nbScol);
clear allS;
magnetic_fields = S_coex(:,1); % or S_bdy, magnetic field values the same for all spectra because of alignment
% reduced all spectra to a vector of intensity values (magnetic fields not needed for fitting)
S_coex = S_coex(:,2:2:end);
S_bdy = S_bdy(:,2:2:end);

% linearly transform ternary(i.e. homogeneous barycentric coordinates) to 
% cartesian coordinates using 1st (Sm) and 3rd (chol) mole fractions
C_bdy = tern2cart(C_bdy,1);
C_coex = tern2cart(C_coex,1);

% determine sizes of inputs before fitting
[ncC,ncd] = size(C_coex);
[ncB,ncS] = size(S_coex);
[nbC,nbd] = size(C_bdy);
[nbB,nbS] = size(S_bdy);

% [u,w,v] = svd([S_coex S_bdy],0);
% S = u(:,1:5)*w(1:5,1:5)*v(:,1:5)';
% S_coex = S(:,1:ncS);
% S_bdy = S(:,ncS+1:end);

disp(sprintf('analyzing baseline...\n'));
% get baseline statistics
for s = 1:ncS
    cSbaseline(:,s) = S_coex([1:200 end-200:end]',s);
    avg_cSbline(s) = mean(cSbaseline(:,s));
    var_cSbline(s) = var(cSbaseline(:,s));
    std_cSbline(s) = std(cSbaseline(:,s));
end 
for s = 1:nbS
    bSbaseline(:,s) = S_bdy([1:200 end-200:end]',s);
    avg_bSbline(s) = mean(bSbaseline(:,s));
    var_bSbline(s) = var(bSbaseline(:,s));
    std_bSbline(s) = std(bSbaseline(:,s));
end 
stdbS = std(S_bdy')';
varbS = var(S_bdy')';

% weight_fxn_S = ones(size(stdbS));

weight_fxn_S = zeros(size(stdbS));
weight_fxn_S(300:ncB-200) = ones(length(300:ncB-200),1);

w = zeros(size(stdbS));
w(stdbS >= 10*mean(std_bSbline)) = 1;
weight_fxn_fp = w;
% weight_fxn_S = w;

% w = zeros(size(stdbS));
% w(weight_fxn >= 10*mean(std_bSbline)) = 1;
% w(weight_fxn >= 100*mean(std_bSbline)) = 2;
% w(weight_fxn >= 1000*mean(std_bSbline)) = 3;
% w(weight_fxn >= 10000*mean(std_bSbline)) = 4;
% weight_fxn = w;

% declare spectra data structure that contains information about the spectra
spectra = struct('magnetic_fields',magnetic_fields,'cSbline_var',var_cSbline,'cSbline_std',std_cSbline,...
    'bSC_std',stdbS,'bSC_var',varbS,'weight_fxn_S',weight_fxn_S,'weight_fxn_fp',weight_fxn_fp);

disp(sprintf('initializing random number generator...\n'));
rand('state',sum(100*clock));

disp(sprintf('parameterizing boundary...\n'));
% chord length parameterization of boundary.
[bdy_length,b] = chord_length(C_bdy); % b = boundary parameter
boundary.parameterization = b;

switch tieline.config
    case 'parallel'
        switch fit_method
            case 'lsqnonlin'
                options = optimset('OutputFcn',@outfun,'MaxFunEvals',6000,'MaxIter',1000,'TolFun',1e-6,'TolX',1e-6,'LargeScale','on','Display','iter','Jacobian','off');
            case 'levmar'
                options = optimset('OutputFcn',@outfun,'MaxFunEvals',6000,'MaxIter',1000,'TolFun',1e-6,'TolX',1e-6,'LargeScale','on','LevenbergMarquardt','on','Display','iter','Jacobian','off');
            case 'simplex'
                options = optimset('OutputFcn',@outfun,'MaxFunEvals',6000,'MaxIter',1000,'TolFun',1e-6,'TolX',1e-6,'Display','iter');
            otherwise
                error('invalid fit method');
        end
        switch boundary.config
            case '2critpts'
                % parallel tielines implicitly specify both critical points;
                % therefore, the position of the critical points will not be
                % searched.  instead the search parameter will be the slope
                % of the tielines.
                if interactive
                    reply = ask4input('constraints',boundary,tieline);
                    tieline.slope_position = reply.reply1;
                    tieline.slope = reply.reply2;
                    m = tieline.slope;
                else
                    m = tieline.slope./180;
                end
                switch tieline.slope_position
                    case 'located'
                        if interactive
                            reply = ask4input('x0',boundary,tieline);
                            x0 = reply;
                        else
                            x0 = [m 0]';
                        end

                        disp(sprintf('starting search and fitting...'));
                        tic;
                        switch fit_method
                            case 'lsqnonlin'
                                lb = [m-1e-8 -Inf]';
                                ub = [m+1e-8 Inf]';
                                [x,resnorm,residual,exitflag] = lsqnonlin(@parallel_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                            case 'levmar'
                            case 'simplex'
                                lb = [m -100]';
                                ub = [m 100]';
                                A = [];
                                b = [];
                                nonlcon = [];
                                [x,fval,exitflag] = fminsearchcon(@parallel_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                            otherwise
                                error('invalid fit method');
                        end

                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                        disp(sprintf('exitflag = %d',exitflag));
                        [fitmeasure,output] = parallel_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
%                         plot_tielines(output);
                        display = [fitmeasure;x];
                        disp(sprintf('best fit measure = %f, with parameters [%f %f]\n',display));
                        disp(sprintf('done!'));
                        toc;
                        return % fitting finished
                    case {'range','unknown'} % slope_position
                        if m(1) > m(end)
                            start = m(end);
                            finish = m(1);
                            m = [start finish];
                        elseif m(1) < m(end)
                            start = m(1);
                            finish = m(end);
                            m = [start finish];
                        else
                            error('slope range is same point');
                        end
                        
                        if interactive
                            reply = ask4input('x0',boundary,tieline);
                            x0 = reply;
                        else
                            x0 = [m(1)+(m(end)-m(1)).*rand(1) 0]';
                        end
                        
                        disp(sprintf('starting search and fitting...'));
                        tic;
                        switch fit_method
                            case 'lsqnonlin'
                                lb = [m(1) -Inf]';
                                ub = [m(end) Inf]';
                                [x,resnorm,residual,exitflag] = lsqnonlin(@parallel_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                            case 'levmar'
                            case 'simplex'
                                lb = [m(1) -100]';
                                ub = [m(end) 100]';
                                A = [];
                                b = [];
                                nonlcon = [];
                                [x,fval,exitflag] = fminsearchcon(@parallel_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                            otherwise
                                error('invalid fit method');
                        end

                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                        disp(sprintf('exitflag = %d',exitflag));
                        [fitmeasure,output] = parallel_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                        plot_tielines(output);
                        display = [fitmeasure;x];
                        disp(sprintf('best fit measure = %f, with parameters [%f %f]\n',display));
                        disp(sprintf('done!'));
                        toc;
                        return % fitting finished
                    otherwise
                        error('invalid slope position');
                end
            case '1critpt/1endtl'
                if interactive
                    reply = ask4input('constraints',boundary,tieline);
                    boundary.endpt1position = reply.reply1;
                    boundary.endpt1parameter = reply.reply2;
                    boundary.endpt2position = reply.reply3;
                    boundary.endpt2parameter = reply.reply4;
                end
                e1 = boundary.endpt1parameter;
                e2 = boundary.endpt2parameter;
                switch boundary.endpt1position
                    case 'located'
                        switch boundary.endpt2position
                            case 'located'
                                if interactive
                                    reply = ask4input('x0',boundary,tieline);
                                    x0 = reply;
                                else
                                    x0 = [e1 e2 0 0]';
                                end
                                
                                disp(sprintf('starting search and fitting...'));
                                tic;
                                switch fit_method
                                    case 'lsqnonlin'
                                        lb = [e1-1e-8 e2-1e-8 -Inf -Inf]';
                                        ub = [e1+1e-8 e2+1e-8 Inf Inf]';
                                        [x,resnorm,residual,exitflag] = lsqnonlin(@parallel_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                    case 'levmar'
                                    case 'simplex'
                                        lb = [e1 e2 -100 -100]';
                                        ub = [e1 e2 100 100]';
                                        A = [0 0 -1 -1];
                                        b = 1-(1e-8);
                                        nonlcon = [];
                                        [x,fval,exitflag] = fminsearchcon(@parallel_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                    otherwise
                                        error('invalid fit method');
                                end

                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                disp(sprintf('exitflag = %d',exitflag));
                                [fitmeasure,output] = parallel_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
%                                 plot_tielines(output);
                                display = [fitmeasure;x];
                                disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f]\n',display));
                                disp(sprintf('done!'));
                                toc;
                                return % fitting finished
                            case {'range','unknown'} % e1 = 'located'
                                if e2(1) > e2(end)
                                    start = e2(end);
                                    finish = e2(1);
                                    e2 = [start finish];
                                elseif e2(1) < e2(end)
                                    start = e2(1);
                                    finish = e2(end);
                                    e2 = [start finish];
                                else
                                    error('end point 2 range is same point');
                                end
                                if interactive
                                    reply = ask4input('x0',boundary,tieline);
                                    x0 = reply;
                                else
                                    x0 = [e1 e2(1)+(e2(end)-e2(1)).*rand(1) 0 0]';
                                end
                                disp(sprintf('starting search and fitting...'));
                                tic;
                                switch fit_method
                                    case 'lsqnonlin'
                                        lb = [e1-1e-8 e2(1) -Inf -Inf]';
                                        ub = [e1+1e-8 e2(end) Inf Inf]';
                                        [x,resnorm,residual,exitflag] = lsqnonlin(@parallel_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                    case 'levmar'
                                    case 'simplex'
                                        lb = [e1 e2(1) -100 -100]';
                                        ub = [e1 e2(end) 100 100]';
                                        A = [0 0 -1 -1];
                                        b = 1-(1e-8);
                                        nonlcon = [];
                                        [x,fval,exitflag] = fminsearchcon(@parallel_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                    otherwise
                                        error('invalid fit method');
                                end

                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                disp(sprintf('exitflag = %d',exitflag));
                                [fitmeasure,output] = parallel_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                plot_tielines(output);
                                display = [fitmeasure;x];
                                disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f]\n',display));
                                disp(sprintf('done!'));
                                toc;
                                return % fitting finished
                            otherwise
                                error('invalid end point 2 position');
                        end
                    case {'range','unknown'} % e1
                        switch boundary.endpt2position
                            case 'located'
                                if e1(1) > e1(end)
                                    start = e1(end);
                                    finish = e1(1);
                                    e1 = [start finish];
                                elseif e1(1) < e1(end)
                                    start = e1(1);
                                    finish = e1(end);
                                    e1 = [start finish];
                                else
                                    error('end point 1 range is same point');
                                end
                                if interactive
                                    reply = ask4input('x0',boundary,tieline);
                                    x0 = reply;
                                else
                                    x0 = [e1(1)+(e1(end)-e1(1)).*rand(1) e2 0 0]';
                                end
                                disp(sprintf('starting search and fitting...'));
                                tic;
                                switch fit_method
                                    case 'lsqnonlin'
                                        lb = [e1(1) e2-1e-8 -Inf -Inf]';
                                        ub = [e1(end) e2+1e-8 Inf Inf]';
                                        [x,resnorm,residual,exitflag] = lsqnonlin(@parallel_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                    case 'levmar'
                                    case 'simplex'
                                        lb = [e1(1) e2 -100 -100]';
                                        ub = [e1(end) e2 100 100]';
                                        A = [0 0 -1 -1];
                                        b = 1-(1e-8);
                                        nonlcon = [];
                                        [x,fval,exitflag] = fminsearchcon(@parallel_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                    otherwise
                                        error('invalid fit method');
                                end

                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                disp(sprintf('exitflag = %d',exitflag));
                                [fitmeasure,output] = parallel_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                plot_tielines(output);
                                display = [fitmeasure;x];
                                disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f]\n',display));
                                disp(sprintf('done!'));
                                toc;
                                return % fitting finished
                            case {'range','unknown'} % e2
                                if e1(1) > e1(end)
                                    start = e1(end);
                                    finish = e1(1);
                                    e1 = [start finish];
                                elseif e1(1) < e1(end)
                                    start = e1(1);
                                    finish = e1(end);
                                    e1 = [start finish];
                                else
                                    error('end point 1 range is same point');
                                end
                                if e2(1) > e2(end)
                                    start = e2(end);
                                    finish = e2(1);
                                    e2 = [start finish];
                                elseif e2(1) < e2(end)
                                    start = e2(1);
                                    finish = e2(end);
                                    e2 = [start finish];
                                else
                                    error('end point 2 range is same point');
                                end
                                if interactive
                                    reply = ask4input('x0',boundary,tieline);
                                    x0 = reply;
                                else
                                    x0 = [e1(1)+(e1(end)-e1(1)).*rand(1) e2(1)+(e2(end)-e2(1)).*rand(1) 0 0]';
                                end
                                disp(sprintf('starting search and fitting...'));
                                tic;
                                switch fit_method
                                    case 'lsqnonlin'
                                        lb = [e1(1) e2(1) -Inf -Inf]';
                                        ub = [e1(end) e2(end) Inf Inf]';
                                        [x,resnorm,residual,exitflag] = lsqnonlin(@parallel_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                    case 'levmar'
                                    case 'simplex'
                                        lb = [e1(1) e2(1) -100 -100]';
                                        ub = [e1(end) e2(end) 100 100]';
                                        A = [0 0 -1 -1];
                                        b = 1-(1e-8);
                                        nonlcon = [];
                                        [x,fval,exitflag] = fminsearchcon(@parallel_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                    otherwise
                                        error('invalid fit method');
                                end

                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                disp(sprintf('exitflag = %d',exitflag));
                                [fitmeasure,output] = parallel_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                plot_tielines(output);
                                display = [fitmeasure;x];
                                disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f]\n',display));
                                disp(sprintf('done!'));
                                toc;
                                return % fitting finished
                            otherwise
                                error('invalid end point 2 position');
                        end
                    otherwise
                        error('invalid end point 1 position');
                end
            case '2endtls'
                % if end tielines parallel, only one configuration;if not,
                % no configuration
            otherwise
                error('invalid boundary configuration');
        end
    case 'tangent'
        switch fit_method
            case 'lsqnonlin'
                options = optimset('OutputFcn',@outfun,'MaxFunEvals',6000,'MaxIter',1000,'TolFun',1e-6,'TolX',1e-6,'LargeScale','on','Display','iter','Jacobian','off');
            case 'levmar'
                options = optimset('OutputFcn',@outfun,'MaxFunEvals',6000,'MaxIter',1000,'TolFun',1e-6,'TolX',1e-6,'LargeScale','on','LevenbergMarquardt','on','Display','iter','Jacobian','off');
            case 'simplex'
                options = optimset('OutputFcn',@outfun,'MaxFunEvals',6000,'MaxIter',1000,'TolFun',1e-6,'TolX',1e-6,'Display','iter');
            otherwise
                error('invalid fit method');
        end
        switch boundary.config
            case '2critpts'
                % boundary has two search parameters b1 and b2 for the
                % position of critical point 1 and critical point 2 on the
                % boundary, respectively.
                if interactive
                    reply = ask4input('constraints',boundary,tieline);
                    boundary.critpt1position = reply.reply1;
                    boundary.critpt1parameter = reply.reply2;
                    boundary.critpt2position = reply.reply3;
                    boundary.critpt2parameter = reply.reply4;
                end
                c1 = boundary.critpt1parameter;
                c2 = boundary.critpt2parameter;
                switch boundary.critpt1position
                    case 'located' 
                        switch boundary.critpt2position
                            case 'located'
                                if interactive
                                    reply = ask4input('x0',boundary,tieline);
                                    x0 = reply;
                                else
                                    x0 = [c1 c2 0]';
                                end
                                
                                disp(sprintf('starting search and fitting...'));
                                tic;
                                switch fit_method
                                    case 'lsqnonlin'
                                        lb = [c1-1e-8 c2-1e-8 -Inf]';
                                        ub = [c1+1e-8 c2+1e-8 Inf]';
                                        [x,resnorm,residual,exitflag] = lsqnonlin(@tangent_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                    case 'levmar'
                                    case 'simplex'
                                        lb = [c1 c2 -100]';
                                        ub = [c1 c2 100]';
                                        A = [];
                                        b = [];
                                        nonlcon = [];
                                        [x,fval,exitflag] = fminsearchcon(@tangent_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                    otherwise
                                        error('invalid fit method');
                                end

                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                disp(sprintf('exitflag = %d',exitflag));
                                [fitmeasure,output] = tangent_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
%                                 plot_tielines(output);
                                display = [fitmeasure;x];
                                disp(sprintf('best fit measure = %f, with parameters [%f %f %f]\n',display));
                                disp(sprintf('done!'));
                                toc;
                                return % fitting finished
                            case {'range','unknown'} % c2
                                if c2(1) > c2(end)
                                    start = c2(end);
                                    finish = c2(1);
                                    c2 = [start finish];
                                elseif c2(1) < c2(end)
                                    start = c2(1);
                                    finish = c2(end);
                                    c2 = [start finish];
                                else
                                    error('critical point 2 range is same point');
                                end
                                if interactive
                                    reply = ask4input('x0',boundary,tieline);
                                    x0 = reply;
                                else
                                    x0 = [c1 c2(1)+(c2(end)-c2(1)).*rand(1) 0]';
                                end
                                
                                disp(sprintf('starting search and fitting...'));
                                tic;
                                switch fit_method
                                    case 'lsqnonlin'
                                        lb = [c1-1e-8 c2(1) -Inf]';
                                        ub = [c1+1e-8 c2(end) Inf]';
                                        [x,resnorm,residual,exitflag] = lsqnonlin(@tangent_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                    case 'levmar'
                                    case 'simplex'
                                        lb = [c1 c2(1) -100]';
                                        ub = [c1 c2(end) 100]';
                                        A = [];
                                        b = [];
                                        nonlcon = [];
                                        [x,fval,exitflag] = fminsearchcon(@tangent_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                    otherwise
                                        error('invalid fit method');
                                end

                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                disp(sprintf('exitflag = %d',exitflag));
                                [fitmeasure,output] = tangent_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                plot_tielines(output);
                                display = [fitmeasure;x];
                                disp(sprintf('best fit measure = %f, with parameters [%f %f %f]\n',display));
                                disp(sprintf('done!'));
                                toc;
                                return % fitting finished
                            otherwise
                                error('invalid critical point 2 position on boundary');
                        end
                    case {'range','unknown'} % c1 
                        switch boundary.critpt2position
                            case 'located'
                                if c1(1) > c1(end)
                                    start = c1(end);
                                    finish = c1(1);
                                    c1 = [start finish];
                                elseif c1(1) < c1(end)
                                    start = c1(1);
                                    finish = c1(end);
                                    c1 = [start finish];
                                else
                                    error('critical point 1 range is same point');
                                end
                                if interactive
                                    reply = ask4input('x0',boundary,tieline);
                                    x0 = reply;
                                else
                                    x0 = [c1(1)+(c1(end)-c1(1)).*rand(1) c2 0]';
                                end
                                
                                disp(sprintf('starting search and fitting...'));
                                tic;
                                switch fit_method
                                    case 'lsqnonlin'
                                        lb = [c1(1) c2-1e-8 -Inf]';
                                        ub = [c1(end) c2+1e-8 Inf]';
                                        [x,resnorm,residual,exitflag] = lsqnonlin(@tangent_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                    case 'levmar'
                                    case 'simplex'
                                        lb = [c1(1) c2 -100]';
                                        ub = [c1(end) c2 100]';
                                        A = [];
                                        b = [];
                                        nonlcon = [];
                                        [x,fval,exitflag] = fminsearchcon(@tangent_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                    otherwise
                                        error('invalid fit method');
                                end

                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                disp(sprintf('exitflag = %d',exitflag));
                                [fitmeasure,output] = tangent_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                plot_tielines(output);
                                display = [fitmeasure;x];
                                disp(sprintf('best fit measure = %f, with parameters [%f %f %f]\n',display));
                                disp(sprintf('done!'));
                                toc;
                                return % fitting finished
                            case {'range','unknown'} % for critical point 2, critical point 1 also range or unknown
                                if c1(1) > c1(end)
                                    start = c1(end);
                                    finish = c1(1);
                                    e1 = [start finish];
                                elseif c1(1) < c1(end)
                                    start = c1(1);
                                    finish = c1(end);
                                    c1 = [start finish];
                                else
                                    error('critical point 1 range is same point');
                                end
                                if c2(1) > c2(end)
                                    start = c2(end);
                                    finish = c2(1);
                                    c2 = [start finish];
                                elseif c2(1) < c2(end)
                                    start = c2(1);
                                    finish = c2(end);
                                    c2 = [start finish];
                                else
                                    error('critical point 2 range is same point');
                                end
                                if interactive
                                    reply = ask4input('x0',boundary,tieline);
                                    x0 = reply;
                                else
                                    x0 = [c1(1)+(c1(end)-c1(1))*rand(1) c2(1)+(c2(end)-c2(1))*rand(1) 0]';
                                end
                                
                                disp(sprintf('starting search and fitting...'));
                                tic;
                                switch fit_method
                                    case 'lsqnonlin'
                                        lb = [c1(1) c2(1) -Inf]';
                                        ub = [c1(end) c2(end) Inf]';
                                        [x,resnorm,residual,exitflag] = lsqnonlin(@tangent_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                    case 'levmar'
                                    case 'simplex'
                                        lb = [c1(1) c2(1) -100]';
                                        ub = [c1(end) c2(end) 100]';
                                        A = [];
                                        b = [];
                                        nonlcon = [];
                                        [x,fval,exitflag] = fminsearchcon(@tangent_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                    otherwise
                                        error('invalid fit method');
                                end

                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                disp(sprintf('exitflag = %d',exitflag));
                                [fitmeasure,output] = tangent_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                plot_tielines(output);
                                display = [fitmeasure;x];
                                disp(sprintf('best fit measure = %f, with parameters [%f %f %f]\n',display));
                                disp(sprintf('done!'));
                                toc;
                                return % fitting finished
                            otherwise
                                error('invalid critical point 2 position on boundary');
                        end 
                    otherwise
                        error('invalid critical point 1 position on boundary');
                end 
            case '1critpt/1endtl'
                if interactive
                    reply = ask4input('constraints',boundary,tieline);
                    boundary.critptposition = reply.reply1;
                    boundary.critptparameter = reply.reply2;
                    boundary.endpt1position = reply.reply3;
                    boundary.endpt1parameter = reply.reply4;
                    boundary.endpt2position = reply.reply5;
                    boundary.endpt2parameter = reply.reply6;
                end
                cp = boundary.critptparameter;
                e1 = boundary.endpt1parameter;
                e2 = boundary.endpt2parameter;
                switch boundary.critptposition
                    case 'located'
                        switch boundary.endpt1position
                            case 'located'
                                switch boundary.endpt2position
                                    case 'located'
                                        if interactive
                                            reply = ask4input('x0',boundary,tieline);
                                            x0 = reply;
                                        else
                                            x0 = [cp e1 e2 0 0]';
                                        end
                                        
                                        disp(sprintf('starting search and fitting...'));
                                        tic;
                                        switch fit_method
                                            case 'lsqnonlin'
                                                lb = [cp-1e-8 e1-1e-8 e2-1e-8 -Inf -Inf]';
                                                ub = [cp+1e-8 e1+1e-8 e2+1e-8 Inf Inf]';
                                                [x,resnorm,residual,exitflag] = lsqnonlin(@tangent_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            case 'levmar'
                                            case 'simplex'
                                                lb = [cp e1 e2 -100 -100]';
                                                ub = [cp e1 e2 100 100]';
                                                A = [0 0 0 -1 -1];
                                                b = 1-(1e-8);
                                                nonlcon = [];
                                                [x,fval,exitflag] = fminsearchcon(@tangent_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            otherwise
                                                error('invalid fit method');
                                        end

                                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                                        disp(sprintf('exitflag = %d',exitflag));
                                        [fitmeasure,output] = tangent_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
%                                         plot_tielines(output);
                                        display = [fitmeasure;x];
                                        disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f %f]\n',display));
                                        disp(sprintf('done!'));
                                        toc;
                                        return % fitting finished
                                    case {'range','unknown'} % e2, cp = 'located', e1 = 'located'
                                        if e2(1) > e2(end)
                                            start = e2(end);
                                            finish = e2(1);
                                            e2 = [start finish];
                                        elseif e2(1) < e2(end)
                                            start = e2(1);
                                            finish = e2(end);
                                            e2 = [start finish];
                                        else
                                            error('end point 2 range is same point');
                                        end
                                        if interactive
                                            reply = ask4input('x0',boundary,tieline);
                                            x0 = reply;
                                        else
                                            x0 = [cp e1 e2(1)+(e2(end)-e2(1)).*rand(1) 0 0]';
                                        end
                                        
                                        disp(sprintf('starting search and fitting...'));
                                        tic;
                                        switch fit_method
                                            case 'lsqnonlin'
                                                lb = [cp-1e-8 e1-1e-8 e2(1) -Inf -Inf]';
                                                ub = [cp+1e-8 e1+1e-8 e2(end) Inf Inf]';
                                                [x,resnorm,residual,exitflag] = lsqnonlin(@tangent_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            case 'levmar'
                                            case 'simplex'
                                                lb = [cp e1 e2(1) -100 -100]';
                                                ub = [cp e1 e2(end) 100 100]';
                                                A = [0 0 0 -1 -1];
                                                b = 1-(1e-8);
                                                nonlcon = [];
                                                [x,fval,exitflag] = fminsearchcon(@tangent_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            otherwise
                                                error('invalid fit method');
                                        end

                                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                                        disp(sprintf('exitflag = %d',exitflag));
                                        [fitmeasure,output] = tangent_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                        plot_tielines(output);
                                        display = [fitmeasure;x];
                                        disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f %f]\n',display));
                                        disp(sprintf('done!'));
                                        toc;
                                        return % fitting finished
                                    otherwise
                                        error('invalid end point 2 position');
                                end
                            case {'range','unknown'} % e1, cp = 'located'
                                switch boundary.endpt2position
                                    case 'located'
                                        if e1(1) > e1(end)
                                            start = e1(end);
                                            finish = e1(1);
                                            e1 = [start finish];
                                        elseif e1(1) < e1(end)
                                            start = e1(1);
                                            finish = e1(end);
                                            e1 = [start finish];
                                        else
                                            error('end point 1 range is same point');
                                        end
                                        if interactive
                                            reply = ask4input('x0',boundary,tieline);
                                            x0 = reply;
                                        else
                                            x0 = [cp e1(1)+(e1(end)-e1(1)).*rand(1) e2 0 0]';
                                        end
                                         
                                        disp(sprintf('starting search and fitting...'));
                                        tic;
                                        switch fit_method
                                            case 'lsqnonlin'
                                                lb = [cp-1e-8 e1(1) e2-1e-8 -Inf -Inf]';
                                                ub = [cp+1e-8 e1(end) e2+1e-8 Inf Inf]';
                                                [x,resnorm,residual,exitflag] = lsqnonlin(@tangent_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            case 'levmar'
                                            case 'simplex'
                                                lb = [cp e1(1) e2 -100 -100]';
                                                ub = [cp e1(end) e2 100 100]';
                                                A = [0 0 0 -1 -1];
                                                b = 1-(1e-8);
                                                nonlcon = [];
                                                [x,fval,exitflag] = fminsearchcon(@tangent_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            otherwise
                                                error('invalid fit method');
                                        end

                                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                                        disp(sprintf('exitflag = %d',exitflag));
                                        [fitmeasure,output] = tangent_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                        plot_tielines(output);
                                        display = [fitmeasure;x];
                                        disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f %f]\n',display));
                                        disp(sprintf('done!'));
                                        toc;
                                        return % fitting finished
                                    case {'range','unknown'} % e2, cp = 'located', e1 = 'range'/'unknown'
                                        if e1(1) > e1(end)
                                            start = e1(end);
                                            finish = e1(1);
                                            e1 = [start finish];
                                        elseif e1(1) < e1(end)
                                            start = e1(1);
                                            finish = e1(end);
                                            e1 = [start finish];
                                        else
                                            error('end point 1 range is same point');
                                        end
                                        if e2(1) > e2(end)
                                            start = e2(end);
                                            finish = e2(1);
                                            e2 = [start finish];
                                        elseif e2(1) < e2(end)
                                            start = e2(1);
                                            finish = e2(end);
                                            e2 = [start finish];
                                        else
                                            error('end point 2 range is same point');
                                        end
                                        if interactive
                                            reply = ask4input('x0',boundary,tieline);
                                            x0 = reply;
                                        else
                                            x0 = [cp e1(1)+(e1(end)-e1(1)).*rand(1) e2(1)+(e2(end)-e2(1)).*rand(1) 0 0]';
                                        end
                                        
                                        disp(sprintf('starting search and fitting...'));
                                        tic;
                                        switch fit_method
                                            case 'lsqnonlin'
                                                lb = [cp-1e-8 e1(1) e2(1) -Inf -Inf]';
                                                ub = [cp+1e-8 e1(end) e2(end) Inf Inf]';
                                                [x,resnorm,residual,exitflag] = lsqnonlin(@tangent_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            case 'levmar'
                                            case 'simplex'
                                                lb = [cp e1(1) e2(1) -100 -100]';
                                                ub = [cp e1(end) e2(end) 100 100]';
                                                A = [0 0 0 -1 -1];
                                                b = 1-(1e-8);
                                                nonlcon = [];
                                                [x,fval,exitflag] = fminsearchcon(@tangent_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            otherwise
                                                error('invalid fit method');
                                        end

                                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                                        disp(sprintf('exitflag = %d',exitflag));
                                        [fitmeasure,output] = tangent_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                        plot_tielines(output);
                                        display = [fitmeasure;x];
                                        disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f %f]\n',display));
                                        disp(sprintf('done!'));
                                        toc;
                                        return % fitting finished
                                    otherwise
                                        error('invalid end point 2 position');
                                end
                            otherwise
                                error('invalid end point 1 position');
                        end
                    case {'range','unknown'} % cp
                        if cp(1) > cp(end)
                            start = cp(end);
                            finish = cp(1);
                            cp = [start finish];
                        elseif cp(1) < cp(end)
                            start = cp(1);
                            finish = cp(end);
                            cp = [start finish];
                        else
                            error('critical point range is same point');
                        end
                        switch boundary.endpt1position
                            case 'located'
                                switch boundary.endpt2position
                                    case 'located'
                                        if interactive
                                            reply = ask4input('x0',boundary,tieline);
                                            x0 = reply;
                                        else
                                            x0 = [cp(1)+(cp(end)-cp(1)).*rand(1) e1 e2 0 0]';
                                        end
                                        
                                        disp(sprintf('starting search and fitting...'));
                                        tic;
                                        switch fit_method
                                            case 'lsqnonlin'
                                                lb = [cp(1) e1-1e-8 e2-1e-8 -Inf -Inf]';
                                                ub = [cp(end) e1+1e-8 e2+1e-8 Inf Inf]';
                                                [x,resnorm,residual,exitflag] = lsqnonlin(@tangent_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            case 'levmar'
                                            case 'simplex'
                                                lb = [cp(1) e1 e2 -100 -100]';
                                                ub = [cp(end) e1 e2 100 100]';
                                                A = [0 0 0 -1 -1];
                                                b = 1-(1e-8);
                                                nonlcon = [];
                                                [x,fval,exitflag] = fminsearchcon(@tangent_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            otherwise
                                                error('invalid fit method');
                                        end

                                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                                        disp(sprintf('exitflag = %d',exitflag));
                                        [fitmeasure,output] = tangent_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                        plot_tielines(output);
                                        display = [fitmeasure;x];
                                        disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f %f]\n',display));
                                        disp(sprintf('done!'));
                                        toc;
                                        return % fitting finished
                                    case {'range','unknown'} % e2, cp = 'range'/'unknown', e1 = 'located'
                                        if e2(1) > e2(end)
                                            start = e2(end);
                                            finish = e2(1);
                                            e2 = [start finish];
                                        elseif e2(1) < e2(end)
                                            start = e2(1);
                                            finish = e2(end);
                                            e2 = [start finish];
                                        else
                                            error('end point 2 range is same point');
                                        end
                                        if interactive
                                            reply = ask4input('x0',boundary,tieline);
                                            x0 = reply;
                                        else
                                            x0 = [cp(1)+(cp(end)-cp(1)).*rand(1) e1 e2(1)+(e2(end)-e2(1)).*rand(1) 0 0]';
                                        end
                                        
                                        disp(sprintf('starting search and fitting...'));
                                        tic;
                                        switch fit_method
                                            case 'lsqnonlin'
                                                lb = [cp(1) e1-1e-8 e2(1) -Inf -Inf]';
                                                ub = [cp(end) e1+1e-8 e2(end) Inf Inf]';
                                                [x,resnorm,residual,exitflag] = lsqnonlin(@tangent_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            case 'levmar'
                                            case 'simplex'
                                                lb = [cp(1) e1 e2(1) -100 -100]';
                                                ub = [cp(end) e1 e2(end) 100 100]';
                                                A = [0 0 0 -1 -1];
                                                b = 1-(1e-8);
                                                nonlcon = [];
                                                [x,fval,exitflag] = fminsearchcon(@tangent_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            otherwise
                                                error('invalid fit method');
                                        end

                                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                                        disp(sprintf('exitflag = %d',exitflag));
                                        [fitmeasure,output] = tangent_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                        plot_tielines(output);
                                        display = [fitmeasure;x];
                                        disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f %f]\n',display));
                                        disp(sprintf('done!'));
                                        toc;
                                        return % fitting finished
                                    otherwise
                                        error('invalid end point 2 position');
                                end
                            case {'range','unknown'} % e1, cp = 'range'/'unknown'
                                switch boundary.endpt2position
                                    case 'located'
                                        if e1(1) > e1(end)
                                            start = e1(end);
                                            finish = e1(1);
                                            e1 = [start finish];
                                        elseif e1(1) < e1(end)
                                            start = e1(1);
                                            finish = e1(end);
                                            e1 = [start finish];
                                        else
                                            error('end point 1 range is same point');
                                        end
                                        if interactive
                                            reply = ask4input('x0',boundary,tieline);
                                            x0 = reply;
                                        else
                                            x0 = [cp(1)+(cp(end)-cp(1)).*rand(1) e1(1)+(e1(end)-e1(1)).*rand(1) e2 0 0]';
                                        end
                                        
                                        disp(sprintf('starting search and fitting...'));
                                        tic;
                                        switch fit_method
                                            case 'lsqnonlin'
                                                lb = [cp(1) e1(1) e2-1e-8 -Inf -Inf]';
                                                ub = [cp(end) e1(end) e2+1e-8 Inf Inf]';
                                                [x,resnorm,residual,exitflag] = lsqnonlin(@tangent_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            case 'levmar'
                                            case 'simplex'
                                                lb = [cp(1) e1(1) e2 -100 -100]';
                                                ub = [cp(end) e1(end) e2 100 100]';
                                                A = [0 0 0 -1 -1];
                                                b = 1-(1e-8);
                                                nonlcon = [];
                                                [x,fval,exitflag] = fminsearchcon(@tangent_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            otherwise
                                                error('invalid fit method');
                                        end

                                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                                        disp(sprintf('exitflag = %d',exitflag));
                                        [fitmeasure,output] = tangent_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                        plot_tielines(output);
                                        display = [fitmeasure;x];
                                        disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f %f]\n',display));
                                        disp(sprintf('done!'));
                                        toc;
                                        return % fitting finished
                                    case {'range','unknown'} % e2, cp = 'range'/'unknown', e1 = 'range'/'unknown'
                                        if e1(1) > e1(end)
                                            start = e1(end);
                                            finish = e1(1);
                                            e1 = [start finish];
                                        elseif e1(1) < e1(end)
                                            start = e1(1);
                                            finish = e1(end);
                                            e1 = [start finish];
                                        else
                                            error('end point 1 range is same point');
                                        end
                                        if e2(1) > e2(end)
                                            start = e2(end);
                                            finish = e2(1);
                                            e2 = [start finish];
                                        elseif e2(1) < e2(end)
                                            start = e2(1);
                                            finish = e2(end);
                                            e2 = [start finish];
                                        else
                                            error('end point 2 range is same point');
                                        end
                                        if interactive
                                            reply = ask4input('x0',boundary,tieline);
                                            x0 = reply;
                                        else
                                            x0 = [cp(1)+(cp(end)-cp(1)).*rand(1) e1(1)+(e1(end)-e1(1)).*rand(1) e2(1)+(e2(end)-e2(1)).*rand(1) 0 0]';
                                        end
                                        
                                        disp(sprintf('starting search and fitting...'));
                                        tic;
                                        switch fit_method
                                            case 'lsqnonlin'
                                                lb = [cp(1) e1(1) e2(1) -Inf -Inf]';
                                                ub = [cp(end) e1(end) e2(end) Inf Inf]';
                                                [x,resnorm,residual,exitflag] = lsqnonlin(@tangent_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            case 'levmar'
                                            case 'simplex'
                                                lb = [cp(1) e1(1) e2(1) -100 -100]';
                                                ub = [cp(end) e1(end) e2(end) 100 100]';
                                                A = [0 0 0 -1 -1];
                                                b = 1-(1e-8);
                                                nonlcon = [];
                                                [x,fval,exitflag] = fminsearchcon(@tangent_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            otherwise
                                                error('invalid fit method');
                                        end

                                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                                        disp(sprintf('exitflag = %d',exitflag));
                                        [fitmeasure,output] = tangent_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                        plot_tielines(output);
                                        display = [fitmeasure;x];
                                        disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f %f]\n',display));
                                        disp(sprintf('done!'));
                                        toc;
                                        return % fitting finished
                                    otherwise
                                        error('invalid end point 2 position');
                                end
                            otherwise
                                error('invalid end point 1 position');
                        end
                    otherwise
                        error('invalid critical point position');
                end
            case '2endtls'
                % only one configuration
            otherwise
                error('invalid boundary configuration');
        end
    case 'ruled'
        switch fit_method
            case 'lsqnonlin'
                options = optimset('OutputFcn',@outfun,'MaxFunEvals',6000,'MaxIter',1000,'TolFun',1e-6,'TolX',1e-6,'LargeScale','on','Display','iter','Jacobian','off');
            case 'levmar'
                options = optimset('OutputFcn',@outfun,'MaxFunEvals',6000,'MaxIter',1000,'TolFun',1e-6,'TolX',1e-6,'LargeScale','on','LevenbergMarquardt','on','Display','iter','Jacobian','off');
            case 'simplex'
                options = optimset('OutputFcn',@outfun,'MaxFunEvals',6000,'MaxIter',1000,'TolFun',1e-6,'TolX',1e-6,'Display','iter');
            otherwise
                error('invalid fit method');
        end
        switch boundary.config
            case '2critpts'
                % boundary has two search parameters b1 and b2 for the
                % position of critical point 1 and critical point 2 on the
                % boundary, respectively.
                if interactive
                    reply = ask4input('constraints',boundary,tieline);
                    boundary.critpt1position = reply.reply1;
                    boundary.critpt1parameter = reply.reply2;
                    boundary.critpt2position = reply.reply3;
                    boundary.critpt2parameter = reply.reply4;
                end
                c1 = boundary.critpt1parameter;
                c2 = boundary.critpt2parameter;
                switch boundary.critpt1position
                    case 'located' 
                        switch boundary.critpt2position
                            case 'located'
                                if interactive
                                    reply = ask4input('x0',boundary,tieline);
                                    x0 = reply;
                                else
                                    x0 = [c1 c2 1 0]';
                                end
                                
                                disp(sprintf('starting search and fitting...'));
                                tic;
                                switch fit_method
                                    case 'lsqnonlin'
                                        lb = [c1-1e-8 c2-1e-8 1e-10 -Inf]';
                                        ub = [c1+1e-8 c2+1e-8 Inf Inf]';
                                        [x,resnorm,residual,exitflag] = lsqnonlin(@ruled_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                    case 'levmar'
                                    case 'simplex'
                                        lb = [c1 c2 1e-10 -100]';
                                        ub = [c1 c2 100 100]';
                                        A = [];
                                        b = [];
                                        nonlcon = [];
                                        [x,fval,exitflag] = fminsearchcon(@ruled_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                    otherwise
                                        error('invalid fit method');
                                end

                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                disp(sprintf('exitflag = %d',exitflag));
                                [fitmeasure,output] = ruled_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
%                                 plot_tielines(output);
                                display = [fitmeasure;x];
                                disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f]\n',display));
                                disp(sprintf('done!'));
                                toc;
                                return % fitting finished   
                            case {'range','unknown'} % c2
                                if c2(1) > c2(end)
                                    start = c2(end);
                                    finish = c2(1);
                                    c2 = [start finish];
                                elseif c2(1) < c2(end)
                                    start = c2(1);
                                    finish = c2(end);
                                    c2 = [start finish];
                                else
                                    error('critical point 2 range is same point');
                                end
                                if interactive
                                    reply = ask4input('x0',boundary,tieline);
                                    x0 = reply;
                                else
                                    x0 = [c1 c2(1)+(c2(end)-c2(1)).*rand(1) 1 0]';
                                end
                                
                                disp(sprintf('starting search and fitting...'));
                                tic;
                                switch fit_method
                                    case 'lsqnonlin'
                                        lb = [c1-1e-8 c2(1) 1e-10 -Inf]';
                                        ub = [c1+1e-8 c2(end) Inf Inf]';
                                        [x,resnorm,residual,exitflag] = lsqnonlin(@ruled_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                    case 'levmar'
                                    case 'simplex'
                                        lb = [c1 c2(1) 1e-10 -100]';
                                        ub = [c1 c2(end) 100 100]';
                                        A = [];
                                        b = [];
                                        nonlcon = [];
                                        [x,fval,exitflag] = fminsearchcon(@ruled_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                    otherwise
                                        error('invalid fit method');
                                end

                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                disp(sprintf('exitflag = %d',exitflag));
                                [fitmeasure,output] = ruled_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                plot_tielines(output);
                                display = [fitmeasure;x];
                                disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f]\n',display));
                                disp(sprintf('done!'));
                                toc;
                                return % fitting finished
                            otherwise
                                error('invalid critical point 2 position on boundary');
                        end
                    case {'range','unknown'} % c1
                        switch boundary.critpt2position
                            case 'located'
                                if c1(1) > c1(end)
                                    start = c1(end);
                                    finish = c1(1);
                                    c1 = [start finish];
                                elseif c1(1) < c1(end)
                                    start = c1(1);
                                    finish = c1(end);
                                    c1 = [start finish];
                                else
                                    error('critical point 1 range is same point');
                                end
                                if interactive
                                    reply = ask4input('x0',boundary,tieline);
                                    x0 = reply;
                                else
                                    x0 = [c1(1)+(c1(end)-c1(1)).*rand(1) c2 1 0]';
                                end
                                
                                disp(sprintf('starting search and fitting...'));
                                tic;
                                switch fit_method
                                    case 'lsqnonlin'
                                        lb = [c1(1) c2-1e-8 1e-10 -Inf]';
                                        ub = [c1(end) c2+1e-8 Inf Inf]';
                                        [x,resnorm,residual,exitflag] = lsqnonlin(@ruled_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                    case 'levmar'
                                    case 'simplex'
                                        lb = [c1(1) c2 1e-10 -100]';
                                        ub = [c1(end) c2 100 100]';
                                        A = [];
                                        b = [];
                                        nonlcon = [];
                                        [x,fval,exitflag] = fminsearchcon(@ruled_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                    otherwise
                                        error('invalid fit method');
                                end

                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                disp(sprintf('exitflag = %d',exitflag));
                                [fitmeasure,output] = ruled_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                plot_tielines(output);
                                display = [fitmeasure;x];
                                disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f]\n',display));
                                disp(sprintf('done!'));
                                toc;
                                return % fitting finished
                            case {'range','unknown'} % for critical point 2, critical point 1 also range or unknown
                                if c1(1) > c1(end)
                                    start = c1(end);
                                    finish = c1(1);
                                    e1 = [start finish];
                                elseif c1(1) < c1(end)
                                    start = c1(1);
                                    finish = c1(end);
                                    c1 = [start finish];
                                else
                                    error('critical point 1 range is same point');
                                end
                                if c2(1) > c2(end)
                                    start = c2(end);
                                    finish = c2(1);
                                    c2 = [start finish];
                                elseif c2(1) < c2(end)
                                    start = c2(1);
                                    finish = c2(end);
                                    c2 = [start finish];
                                else
                                    error('critical point 2 range is same point');
                                end
                                if interactive
                                    reply = ask4input('x0',boundary,tieline);
                                    x0 = reply;
                                else
                                    x0 = [c1(1)+(c1(end)-c1(1)).*rand(1) c2(1)+(c2(end)-c2(1)).*rand(1) 1 0]';
                                end
                                
                                disp(sprintf('starting search and fitting...'));
                                tic;
                                switch fit_method
                                    case 'lsqnonlin'
                                        lb = [c1(1) c2(1) 1e-10 -Inf]';
                                        ub = [c1(end) c2(end) Inf Inf]';
                                        [x,resnorm,residual,exitflag] = lsqnonlin(@ruled_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                    case 'levmar'
                                    case 'simplex'
                                        lb = [c1(1) c2(1) 1e-10 -100]';
                                        ub = [c1(end) c2(end) 100 100]';
                                        A = [];
                                        b = [];
                                        nonlcon = [];
                                        [x,fval,exitflag] = fminsearchcon(@ruled_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                    otherwise
                                        error('invalid fit method');
                                end

                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                disp(sprintf('exitflag = %d',exitflag));
                                [fitmeasure,output] = ruled_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                plot_tielines(output);
                                display = [fitmeasure;x];
                                disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f]\n',display));
                                disp(sprintf('done!'));
                                toc;
                                return % fitting finished
                            otherwise
                                error('invalid critical point 2 position on boundary');
                        end 
                    otherwise
                        error('invalid critical point 1 position on boundary');
                end
            case '1critpt/1endtl'
                if interactive
                    reply = ask4input('constraints',boundary,tieline);
                    boundary.critptposition = reply.reply1;
                    boundary.critptparameter = reply.reply2;
                    boundary.endpt1position = reply.reply3;
                    boundary.endpt1parameter = reply.reply4;
                    boundary.endpt2position = reply.reply5;
                    boundary.endpt2parameter = reply.reply6;
                end
                cp = boundary.critptparameter;
                e1 = boundary.endpt1parameter;
                e2 = boundary.endpt2parameter;
                switch boundary.critptposition
                    case 'located'
                        switch boundary.endpt1position
                            case 'located'
                                switch boundary.endpt2position
                                    case 'located'
                                        if interactive
                                            reply = ask4input('x0',boundary,tieline);
                                            x0 = reply;
                                        else
                                            x0 = [cp e1 e2 1 0 0]';
                                        end
                                        
                                        disp(sprintf('starting search and fitting...'));
                                        tic;
                                        switch fit_method
                                            case 'lsqnonlin'
                                                lb = [cp-1e-8 e1-1e-8 e2-1e-8 1e-10 -Inf -Inf]';
                                                ub = [cp+1e-8 e1+1e-8 e2+1e-8 Inf Inf Inf]';
                                                [x,resnorm,residual,exitflag] = lsqnonlin(@ruled_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            case 'levmar'
                                            case 'simplex'
                                                lb = [cp e1 e2 1e-10 -100 -100]';
                                                ub = [cp e1 e2 100 100 100]';
                                                A = [0 0 0 0 -1 -1];
                                                b = 1-(1e-8);
                                                nonlcon = [];
                                                [x,fval,exitflag] = fminsearchcon(@ruled_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            otherwise
                                                error('invalid fit method');
                                        end

                                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                                        disp(sprintf('exitflag = %d',exitflag));
                                        [fitmeasure,output] = ruled_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
%                                         plot_tielines(output);
                                        display = [fitmeasure;x];
                                        disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f %f %f]\n',display));
                                        disp(sprintf('done!'));
                                        toc;
                                        return % fitting finished
                                    case {'range','unknown'} % e2, cp = 'located', e1 = 'located'
                                        if e2(1) > e2(end)
                                            start = e2(end);
                                            finish = e2(1);
                                            e2 = [start finish];
                                        elseif e2(1) < e2(end)
                                            start = e2(1);
                                            finish = e2(end);
                                            e2 = [start finish];
                                        else
                                            error('end point 2 range is same point');
                                        end
                                        if interactive
                                            reply = ask4input('x0',boundary,tieline);
                                            x0 = reply;
                                        else
                                            x0 = [cp e1 e2(1)+(e2(end)-e2(1)).*rand(1) 1 0 0]';
                                        end
                                        
                                        disp(sprintf('starting search and fitting...'));
                                        tic;
                                        switch fit_method
                                            case 'lsqnonlin'
                                                lb = [cp-1e-8 e1-1e-8 e2(1) 1e-10 -Inf -Inf]';
                                                ub = [cp+1e-8 e1+1e-8 e2(end) Inf Inf Inf]';
                                                [x,resnorm,residual,exitflag] = lsqnonlin(@ruled_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            case 'levmar'
                                            case 'simplex'
                                                lb = [cp e1 e2(1) 1e-10 -100 -100]';
                                                ub = [cp e1 e2(end) 100 100 100]';
                                                A = [0 0 0 0 -1 -1];
                                                b = 1-(1e-8);
                                                nonlcon = [];
                                                [x,fval,exitflag] = fminsearchcon(@ruled_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            otherwise
                                                error('invalid fit method');
                                        end

                                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                                        disp(sprintf('exitflag = %d',exitflag));
                                        [fitmeasure,output] = ruled_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                        plot_tielines(output);
                                        display = [fitmeasure;x];
                                        disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f %f %f]\n',display));
                                        disp(sprintf('done!'));
                                        toc;
                                        return % fitting finished
                                    otherwise
                                        error('invalid end point 2 position');
                                end
                            case {'range','unknown'} % e1, cp = 'located'
                                switch boundary.endpt2position
                                    case 'located'
                                        if e1(1) > e1(end)
                                            start = e1(end);
                                            finish = e1(1);
                                            e1 = [start finish];
                                        elseif e1(1) < e1(end)
                                            start = e1(1);
                                            finish = e1(end);
                                            e1 = [start finish];
                                        else
                                            error('end point 1 range is same point');
                                        end
                                        if interactive
                                            reply = ask4input('x0',boundary,tieline);
                                            x0 = reply;
                                        else
                                            x0 = [cp e1(1)+(e1(end)-e1(1)).*rand(1) e2 1 0 0]';
                                        end
                                        
                                        disp(sprintf('starting search and fitting...'));
                                        tic;
                                        switch fit_method
                                            case 'lsqnonlin'
                                                lb = [cp-1e-8 e1(1) e2-1e-8 1e-10 -Inf -Inf]';
                                                ub = [cp+1e-8 e1(end) e2+1e-8 Inf Inf Inf]';
                                                [x,resnorm,residual,exitflag] = lsqnonlin(@ruled_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            case 'levmar'
                                            case 'simplex'
                                                lb = [cp e1(1) e2 1e-10 -100 -100]';
                                                ub = [cp e1(end) e2 100 100 100]';
                                                A = [0 0 0 0 -1 -1];
                                                b = 1-(1e-8);
                                                nonlcon = [];
                                                [x,fval,exitflag] = fminsearchcon(@ruled_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            otherwise
                                                error('invalid fit method');
                                        end

                                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                                        disp(sprintf('exitflag = %d',exitflag));
                                        [fitmeasure,output] = ruled_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                        plot_tielines(output);
                                        display = [fitmeasure;x];
                                        disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f %f %f]\n',display));
                                        disp(sprintf('done!'));
                                        toc;
                                        return % fitting finished
                                    case {'range','unknown'} % e2, cp = 'located', e1 = 'range'/'unknown'
                                        if e1(1) > e1(end)
                                            start = e1(end);
                                            finish = e1(1);
                                            e1 = [start finish];
                                        elseif e1(1) < e1(end)
                                            start = e1(1);
                                            finish = e1(end);
                                            e1 = [start finish];
                                        else
                                            error('end point 1 range is same point');
                                        end
                                        if e2(1) > e2(end)
                                            start = e2(end);
                                            finish = e2(1);
                                            e2 = [start finish];
                                        elseif e2(1) < e2(end)
                                            start = e2(1);
                                            finish = e2(end);
                                            e2 = [start finish];
                                        else
                                            error('end point 2 range is same point');
                                        end
                                        if interactive
                                            reply = ask4input('x0',boundary,tieline);
                                            x0 = reply;
                                        else
                                            x0 = [cp e1(1)+(e1(end)-e1(1)).*rand(1) e2(1)+(e2(end)-e2(1)).*rand(1) 1 0 0]';
                                        end
                                        
                                        disp(sprintf('starting search and fitting...'));
                                        tic;
                                        switch fit_method
                                            case 'lsqnonlin'
                                                lb = [cp-1e-8 e1(1) e2(1) 1e-10 -Inf -Inf]';
                                                ub = [cp+1e-8 e1(end) e2(end) Inf Inf Inf]';
                                                [x,resnorm,residual,exitflag] = lsqnonlin(@ruled_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            case 'levmar'
                                            case 'simplex'
                                                lb = [cp e1(1) e2(1) 1e-10 -100 -100]';
                                                ub = [cp e1(end) e2(end) 100 100 100]';
                                                A = [0 0 0 0 -1 -1];
                                                b = 1-(1e-8);
                                                nonlcon = [];
                                                [x,fval,exitflag] = fminsearchcon(@ruled_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            otherwise
                                                error('invalid fit method');
                                        end

                                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                                        disp(sprintf('exitflag = %d',exitflag));
                                        [fitmeasure,output] = ruled_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                        plot_tielines(output);
                                        display = [fitmeasure;x];
                                        disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f %f %f]\n',display));
                                        disp(sprintf('done!'));
                                        toc;
                                        return % fitting finished
                                    otherwise
                                        error('invalid end point 2 position');
                                end
                            otherwise
                                error('invalid end point 1 position');
                        end
                    case {'range','unknown'} % cp
                        if cp(1) > cp(end)
                            start = cp(end);
                            finish = cp(1);
                            cp = [start finish];
                        elseif cp(1) < cp(end)
                            start = cp(1);
                            finish = cp(end);
                            cp = [start finish];
                        else
                            error('critical point range is same point');
                        end
                        switch boundary.endpt1position
                            case 'located'
                                switch boundary.endpt2position
                                    case 'located'
                                        if interactive
                                            reply = ask4input('x0',boundary,tieline);
                                            x0 = reply;
                                        else
                                            x0 = [cp(1)+(cp(end)-cp(1)).*rand(1) e1 e2 1 0 0]';
                                        end
                                        
                                        disp(sprintf('starting search and fitting...'));
                                        tic;
                                        switch fit_method
                                            case 'lsqnonlin'
                                                lb = [cp(1) e1-1e-8 e2-1e-8 1e-10 -Inf -Inf]';
                                                ub = [cp(end) e1+1e-8 e2+1e-8 Inf Inf Inf]';
                                                [x,resnorm,residual,exitflag] = lsqnonlin(@ruled_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            case 'levmar'
                                            case 'simplex'
                                                lb = [cp(1) e1 e2 1e-10 -100 -100]';
                                                ub = [cp(end) e1 e2 100 100 100]';
                                                A = [0 0 0 0 -1 -1];
                                                b = 1-(1e-8);
                                                nonlcon = [];
                                                [x,fval,exitflag] = fminsearchcon(@ruled_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            otherwise
                                                error('invalid fit method');
                                        end

                                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                                        disp(sprintf('exitflag = %d',exitflag));
                                        [fitmeasure,output] = ruled_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                        plot_tielines(output);
                                        display = [fitmeasure;x];
                                        disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f %f %f]\n',display));
                                        disp(sprintf('done!'));
                                        toc;
                                        return % fitting finished
                                    case {'range','unknown'} % e2, cp = 'range'/'unknown', e1 = 'located'
                                        if e2(1) > e2(end)
                                            start = e2(end);
                                            finish = e2(1);
                                            e2 = [start finish];
                                        elseif e2(1) < e2(end)
                                            start = e2(1);
                                            finish = e2(end);
                                            e2 = [start finish];
                                        else
                                            error('end point 2 range is same point');
                                        end
                                        if interactive
                                            reply = ask4input('x0',boundary,tieline);
                                            x0 = reply;
                                        else
                                            x0 = [cp(1)+(cp(end)-cp(1)).*rand(1) e1 e2(1)+(e2(end)-e2(1)).*rand(1) 1 0 0]';
                                        end
                                        
                                        disp(sprintf('starting search and fitting...'));
                                        tic;
                                        switch fit_method
                                            case 'lsqnonlin'
                                                lb = [cp(1) e1-1e-8 e2(1) 1e-10 -Inf -Inf]';
                                                ub = [cp(end) e1+1e-8 e2(end) Inf Inf Inf]';
                                                [x,resnorm,residual,exitflag] = lsqnonlin(@ruled_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            case 'levmar'
                                            case 'simplex'
                                                lb = [cp(1) e1 e2(1) 1e-10 -100 -100]';
                                                ub = [cp(end) e1 e2(end) 100 100 100]';
                                                A = [0 0 0 0 -1 -1];
                                                b = 1-(1e-8);
                                                nonlcon = [];
                                                [x,fval,exitflag] = fminsearchcon(@ruled_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            otherwise
                                                error('invalid fit method');
                                        end

                                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                                        disp(sprintf('exitflag = %d',exitflag));
                                        [fitmeasure,output] = ruled_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                        plot_tielines(output);
                                        display = [fitmeasure;x];
                                        disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f %f %f]\n',display));
                                        disp(sprintf('done!'));
                                        toc;
                                        return % fitting finished
                                    otherwise
                                        error('invalid end point 2 position');
                                end
                            case {'range','unknown'} % e1, cp = 'range'/'unknown'
                                switch boundary.endpt2position
                                    case 'located'
                                        if e1(1) > e1(end)
                                            start = e1(end);
                                            finish = e1(1);
                                            e1 = [start finish];
                                        elseif e1(1) < e1(end)
                                            start = e1(1);
                                            finish = e1(end);
                                            e1 = [start finish];
                                        else
                                            error('end point 1 range is same point');
                                        end
                                        if interactive
                                            reply = ask4input('x0',boundary,tieline);
                                            x0 = reply;
                                        else
                                            x0 = [cp(1)+(cp(end)-cp(1)).*rand(1) e1(1)+(e1(end)-e1(1)).*rand(1) e2 1 0 0]';
                                        end
                                        
                                        disp(sprintf('starting search and fitting...'));
                                        tic;
                                        switch fit_method
                                            case 'lsqnonlin'
                                                lb = [cp(1) e1(1) e2-1e-8 1e-10 -Inf -Inf]';
                                                ub = [cp(end) e1(end) e2+1e-8 Inf Inf Inf]';
                                                [x,resnorm,residual,exitflag] = lsqnonlin(@ruled_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            case 'levmar'
                                            case 'simplex'
                                                lb = [cp(1) e1(1) e2 1e-10 -100 -100]';
                                                ub = [cp(end) e1(end) e2 100 100 100]';
                                                A = [0 0 0 0 -1 -1];
                                                b = 1-(1e-8);
                                                nonlcon = [];
                                                [x,fval,exitflag] = fminsearchcon(@ruled_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            otherwise
                                                error('invalid fit method');
                                        end

                                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                                        disp(sprintf('exitflag = %d',exitflag));
                                        [fitmeasure,output] = ruled_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                        plot_tielines(output);
                                        display = [fitmeasure;x];
                                        disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f %f %f]\n',display));
                                        disp(sprintf('done!'));
                                        toc;
                                        return % fitting finished
                                    case {'range','unknown'} % e2, cp = 'range'/'unknown', e1 = 'range'/'unknown'
                                        if e1(1) > e1(end)
                                            start = e1(end);
                                            finish = e1(1);
                                            e1 = [start finish];
                                        elseif e1(1) < e1(end)
                                            start = e1(1);
                                            finish = e1(end);
                                            e1 = [start finish];
                                        else
                                            error('end point 1 range is same point');
                                        end
                                        if e2(1) > e2(end)
                                            start = e2(end);
                                            finish = e2(1);
                                            e2 = [start finish];
                                        elseif e2(1) < e2(end)
                                            start = e2(1);
                                            finish = e2(end);
                                            e2 = [start finish];
                                        else
                                            error('end point 2 range is same point');
                                        end
                                        if interactive
                                            reply = ask4input('x0',boundary,tieline);
                                            x0 = reply;
                                        else
                                            x0 = [cp(1)+(cp(end)-cp(1)).*rand(1) e1(1)+(e1(end)-e1(1)).*rand(1) e2(1)+(e2(end)-e2(1)).*rand(1) 1 0 0]';
                                        end
                                        
                                        disp(sprintf('starting search and fitting...'));
                                        tic;
                                        switch fit_method
                                            case 'lsqnonlin'
                                                lb = [cp(1) e1(1) e2(1) 1e-10 -Inf -Inf]';
                                                ub = [cp(end) e1(end) e2(end) Inf Inf Inf]';
                                                [x,resnorm,residual,exitflag] = lsqnonlin(@ruled_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            case 'levmar'
                                            case 'simplex'
                                                lb = [cp(1) e1(1) e2(1) 1e-10 -100 -100]';
                                                ub = [cp(end) e1(end) e2(end) 100 100 100]';
                                                A = [0 0 0 0 -1 -1];
                                                b = 1-(1e-8);
                                                nonlcon = [];
                                                [x,fval,exitflag] = fminsearchcon(@ruled_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            otherwise
                                                error('invalid fit method');
                                        end

                                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                                        disp(sprintf('exitflag = %d',exitflag));
                                        [fitmeasure,output] = ruled_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
%                                         plot_tielines(output);
                                        display = [fitmeasure;x];
                                        disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f %f %f]\n',display));
                                        disp(sprintf('done!'));
                                        toc;
                                        return % fitting finished
                                    otherwise
                                        error('invalid end point 2 position');
                                end
                            otherwise
                                error('invalid end point 1 position');
                        end
                    otherwise
                        error('invalid critical point position');
                end
            case '2endtls'
            otherwise
                error('invalid boundary configuration');
        end
    otherwise
        error('invalid tieline configuration');
end
        
return