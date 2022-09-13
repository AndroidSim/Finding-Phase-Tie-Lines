function output = global_tieline_fitting(S_coex,C_coex,S_bdy,C_bdy,fit_method,varargin)
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

% if isempty(varargin)
%     fit_method = 'lsqnonlin';
% elseif ischar(varargin{1})
%     if any(strcmp(varargin{1},{'lsqnonlin';'levmar';'simplex'}))
%         fit_method = varargin{1};
%     else
%         error('fit method must be either "lsqnonlin","levmar",or "simplex"');
%     end 
% else
%     error('fit method must be either "lsqnonlin","levmar",or "simplex"');
% end

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
% area of the absorbance spectrum equals one, then convert back to the
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

disp(sprintf('initializing random number generator...\n'));
rand('state',sum(100*clock));

% declare spectra data structure that contains information about the spectra
spectra = struct('magnetic_fields',magnetic_fields,'cSbline_var',var_cSbline,'cSbline_std',std_cSbline,...
    'bSC_std',stdbS,'bSC_var',varbS,'weight_fxn_S',weight_fxn_S,'weight_fxn_fp',weight_fxn_fp);
% declare and initialize boundary, tieline, and Kp search and data structures
boundary = struct('config',0);
tieline = struct('config',0);
part_coeff = struct('config',0);

% ask user for the phase boundary configuration
reply = input(sprintf('"2critpts", "1critpt/1endtl", or "2endtls" boundary configuration?\n'),'s');
while isempty(reply) || ~any(strcmp(reply,{'2critpts','1critpt/1endtl','2endtls'}))
    reply = input(sprintf('"2critpts", "1critpt/1endtl", or "2endtls" boundary configuration?\n'),'s');
end
boundary.config = reply;
% ask user for boundary search options
% reply = input(sprintf('"continuous" or "grid" boundary search method?\n'),'s');
% while isempty(reply) || ~any(strcmp(reply,{'continuous','grid'}))
%     reply = input(sprintf('"continuous" or "grid" boundary search method?\n'),'s');
% end
% boundary.search_method = reply;
% ask user for the tieline configuration
reply = input(sprintf('"parallel", "tangent", or "ruled" tieline configuration?\n'),'s');
while isempty(reply) || ~any(strcmp(reply,{'parallel','tangent','ruled'}))
    reply = input(sprintf('"parallel", "tangent", or "ruled" tieline configuration?\n'),'s');
end
tieline.config = reply;
% ask user for the Kp configuration
reply = input(sprintf('"constant", or "function" Kp configuration?\n'),'s');
while isempty(reply) || ~any(strcmp(reply,{'constant','function'}))
    reply = input(sprintf('"constant", or "function" Kp configuration?\n'),'s');
end
part_coeff.config = reply;

switch boundary.config
    case '2critpts'
        switch tieline.config
            case 'parallel'
                % parallel tielines specify both critical points;
                % therefore, the position of the critical points will not be
                % searched.  instead the search parameter will be the slope
                % of the tielines.
                % chord length parameterization of boundary.
                disp(sprintf('\nparameterizing boundary...\n'));
                [bdy_length,b] = chord_length(C_bdy); % b = boundary parameter
                boundary.parameterization = b;
                
                reply = input(sprintf('is the slope constant? y/n\n'),'s'); 
                if reply == 'y'
                    slope_position = 'located';
                    reply = input(sprintf('what is the slope in angles? [0:180]\n'));
                    while isempty(reply) || ~isscalar(reply) || (reply < 0 | reply > 180)
                        reply = input(sprintf('what is the slope in angles? [0:180]\n'));
                    end
                    m = reply.*(1/180);  
                elseif reply == 'n'
                    reply = input(sprintf('is there a range of slope angles to search? y/n\n'),'s');
                    if reply == 'y'
                        slope_position = 'range';
                        reply = input(sprintf('what is the range of slope angles to search? [0:180] (if unknown, press enter)\n'));
                        if isempty(reply)
                            m = [0 1];
                        elseif isvector(reply) & all(reply >= 0 & reply <= 180)
                            m = reply.*(1/180);
                        else
                            error('invalid response');
                        end
                    elseif reply == 'n'
                        slope_position = 'unknown';
                        m = [0 1];
                    else
                        error('invalid response');
                    end      
                else
                    error('invalid response');
                end  
                
                switch fit_method
                    case 'lsqnonlin'
                        options = optimset('OutputFcn',@outfun,'MaxFunEvals',10000,'MaxIter',1000,'TolFun',1e-6,'TolX',1e-6,'LargeScale','on','Display','iter','Jacobian','off');
                    case 'levmar'
                        options = optimset('OutputFcn',@outfun,'MaxFunEvals',10000,'MaxIter',1000,'TolFun',1e-6,'TolX',1e-6,'LargeScale','on','LevenbergMarquardt','on','Display','iter','Jacobian','off');
                    case 'simplex'
                        options = optimset('OutputFcn',@outfun,'MaxFunEvals',10000,'MaxIter',1000,'TolFun',1e-6,'TolX',1e-6,'Display','iter');
                    otherwise
                        error('invalid fit method');
                end
                
                switch slope_position
                    case 'located'
                        switch part_coeff.config
                            case 'constant' 
                                reply = input(sprintf('what is the range of Kp to search? [1/10:10] (if unspecified, press enter)\n'));
                                if isempty(reply)
                                    kp = [0.1 10];
                                elseif isvector(reply) & all(reply >= 0.1 & reply <= 10)
                                    kp = reply;
                                else
                                    error('invalid response');
                                end
                                reply = input(sprintf('enter 2 parameter starting vector as [%4.3f %4.3f<=kp<=%4.3f]:(enter for default)\n',[m*180;kp(1);kp(end)]));
                                if isempty(reply)
                                    x0 = [m kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                elseif isvector(reply) & length(reply) == 2 & (reply(1) == m*180) & (reply(2) >= kp(1) & reply(2) <= kp(end))
                                    x0 = [reply(1).*(1/180) reply(2)]';
                                else 
                                    error('invalid response');
                                end
                                
                                disp(sprintf('starting search and fitting...'));
                                tic;
                                switch fit_method
                                    case 'lsqnonlin'
                                        lb = [m-1e-8 kp(1)]';
                                        ub = [m+1e-8 kp(end)]';
                                        [x,resnorm,residual,exitflag] = lsqnonlin(@parallel_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                    case 'levmar'
                                    case 'simplex'
                                        lb = [m kp(1)]';
                                        ub = [m kp(end)]';
                                        A = [0 -1];
                                        b = -(1e-8);
                                        nonlcon = [];
                                        [x,fval,exitflag,foutput] = fminsearchcon(@parallel_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
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
                            case 'function' % m = 'located'
                                reply = input(sprintf('enter 2 parameter starting vector as [%4.3f c]:(enter for default)\n',m*180));
                                if isempty(reply)
                                    x0 = [m 0]';
                                elseif isvector(reply) & length(reply) == 2 & (reply(1) == m*180) 
                                    x0 = [reply(1).*(1/180) reply(2)]';
                                else 
                                    error('invalid response');
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
                                plot_tielines(output);
                                display = [fitmeasure;x];
                                disp(sprintf('best fit measure = %f, with parameters [%f %f]\n',display));
                                disp(sprintf('done!'));
                                toc;
                                return % fitting finished
                            otherwise
                                error('invalid Kp configuration');
                        end
                    case {'range','unknown'} % m
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
                        switch part_coeff.config
                            case 'constant' % m = 'range','unknown'
                                reply = input(sprintf('what is the range of Kp to search? [1/10:10] (if unspecified, press enter)\n'));
                                if isempty(reply)
                                    kp = [0.1 10];
                                elseif isvector(reply) & all(reply >= 0.1 & reply <= 10)
                                    kp = reply;
                                else
                                    error('invalid response');
                                end
                                reply = input(sprintf('enter 2 parameter starting vector as [%4.3f<=m<=%4.3f %4.3f<=kp<=%4.3f]:(enter for default)\n',[m(1)*180;m(end)*180;kp(1);kp(end)]));
                                if isempty(reply)
                                    x0 = [m(1)+(m(end)-m(1)).*rand(1) kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                elseif isvector(reply) & length(reply) == 2 & (reply(1) >= m(1)*180 & reply(1) <= m(end)*180) & (reply(2) >= kp(1) & reply(2) <= kp(end))
                                    x0 = [reply(1).*(1/180) reply(2)]';
                                else 
                                    error('invalid response');
                                end
                                
                                disp(sprintf('starting search and fitting...'));
                                tic;
                                switch fit_method
                                    case 'lsqnonlin'
                                        lb = [m(1) kp(1)]';
                                        ub = [m(end) kp(end)]';
                                        [x,resnorm,residual,exitflag] = lsqnonlin(@parallel_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                    case 'levmar'
                                    case 'simplex'
                                        lb = [m(1) kp(1)]';
                                        ub = [m(end) kp(end)]';
                                        A = [0 -1];
                                        b = -(1e-8);
                                        nonlcon = [];
                                        [x,fval,exitflag,foutput] = fminsearchcon(@parallel_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
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
                            case 'function' % m = 'range','unknown'
                                reply = input(sprintf('enter 2 parameter starting vector as [%4.3f<=m<=%4.3f c]:(enter for default)\n',[m(1)*180;m(end)*180]));
                                if isempty(reply)
                                    x0 = [m(1)+(m(end)-m(1)).*rand(1) 0]';
                                elseif isvector(reply) & length(reply) == 2 & (reply(1) >= m(1)*180 & reply(1) <= m(end)*180)
                                    x0 = [reply(1).*(1/180) reply(2)]';
                                else 
                                    error('invalid response');
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
                                error('invalid Kp configuration');
                        end
                    otherwise
                        error('invalid end point 2 position');
                end
            case 'tangent' % 2critpts
                % boundary has two search parameters b1 and b2 for the
                % position of critical point 1 and critical point 2 on the
                % boundary, respectively.
                % loop over number of critical points and get their
                % coordinates.
                % chord length parameterization of boundary.
                disp(sprintf('\nparameterizing boundary...\n'));
                [bdy_length,b] = chord_length(C_bdy); % b = boundary parameter
                boundary.parameterization = b;
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
                
                c1 = boundary.critpt1parameter;
                c2 = boundary.critpt2parameter;
                
                switch fit_method
                    case 'lsqnonlin'
                        options = optimset('OutputFcn',@outfun,'MaxFunEvals',10000,'MaxIter',1000,'TolFun',1e-6,'TolX',1e-6,'LargeScale','on','Display','iter','Jacobian','off');
                    case 'levmar'
                        options = optimset('OutputFcn',@outfun,'MaxFunEvals',10000,'MaxIter',1000,'TolFun',1e-6,'TolX',1e-6,'LargeScale','on','LevenbergMarquardt','on','Display','iter','Jacobian','off');
                    case 'simplex'
                        options = optimset('OutputFcn',@outfun,'MaxFunEvals',10000,'MaxIter',1000,'TolFun',1e-6,'TolX',1e-6,'Display','iter');
                    otherwise
                        error('invalid fit method');
                end
                
                switch boundary.critpt1position
                    case 'located' 
                        switch boundary.critpt2position
                            case 'located'
                                switch part_coeff.config
                                    case 'constant' % c1 = 'located', c2 = 'located'
                                        reply = input(sprintf('what is the range of Kp to search? [1/10:10] (if unspecified, press enter)\n'));
                                        if isempty(reply)
                                            kp = [0.1 10];
                                        elseif isvector(reply) & all(reply >= 0.1 & reply <= 10)
                                            kp = reply;
                                        else
                                            error('invalid response');
                                        end
                                        reply = input(sprintf('enter 3 parameter starting vector as [%4.3f %4.3f %4.3f<=kp<=%4.3f]:(enter for default)\n',[c1;c2;kp(1);kp(end)]));
                                        if isempty(reply)
                                            x0 = [c1 c2 kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                        elseif isvector(reply) & length(reply) == 3 & (reply(1) == c1) & (reply(2) == c2) & (reply(3) >= kp(1) & reply(3) <= kp(end))
                                            x0 = reply';
                                        else 
                                            error('invalid response');
                                        end
                                        
                                        disp(sprintf('starting search and fitting...'));
                                        tic;
                                        switch fit_method
                                            case 'lsqnonlin'
                                                lb = [c1-1e-8 c2-1e-8 kp(1)]';
                                                ub = [c1+1e-8 c2+1e-8 kp(end)]';
                                                [x,resnorm,residual,exitflag] = lsqnonlin(@tangent_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            case 'levmar'
                                            case 'simplex'
                                                lb = [c1 c2 kp(1)]';
                                                ub = [c1 c2 kp(end)]';
                                                A = [0 0 -1];
                                                b = -(1e-8);
                                                nonlcon = [];
                                                [x,fval,exitflag,foutput] = fminsearchcon(@tangent_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
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
                                    case 'function' % c1 = 'located', c2 = 'located'
                                        reply = input(sprintf('enter 3 parameter starting vector as [%4.3f %4.3f c]:(enter for default)\n',[c1;c2]));
                                        if isempty(reply)
                                            x0 = [c1 c2 0]';
                                        elseif isvector(reply) & length(reply) == 3 & (reply(1) == c1) & (reply(2) == c2)
                                            x0 = reply';
                                        else 
                                            error('invalid response');
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
                                        plot_tielines(output);
                                        display = [fitmeasure;x];
                                        disp(sprintf('best fit measure = %f, with parameters [%f %f %f]\n',display));
                                        disp(sprintf('done!'));
                                        toc;
                                        return % fitting finished
                                    otherwise
                                        error('invalid Kp configuration');
                                end
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
                                switch part_coeff.config
                                    case 'constant' % c1 = 'located', c2 = 'range','unknown'
                                        reply = input(sprintf('what is the range of Kp to search? [1/10:10] (if unspecified, press enter)\n'));
                                        if isempty(reply)
                                            kp = [0.1 10];
                                        elseif isvector(reply) & all(reply >= 0.1 & reply <= 10)
                                            kp = reply;
                                        else
                                            error('invalid response');
                                        end
                                        reply = input(sprintf('enter 3 parameter starting vector as [%4.3f %4.3f<=c2<=%4.3f %4.3f<=kp<=%4.3f]:(enter for default)\n',[c1;c2(1);c2(end);kp(1);kp(end)]));
                                        if isempty(reply)
                                            x0 = [c1 c2(1)+(c2(end)-c2(1)).*rand(1) kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                        elseif isvector(reply) & length(reply) == 3 & (reply(1) == c1) & (reply(2) >= c2(1) & reply(2) <= c2(end)) & (reply(3) >= kp(1) & reply(3) <= kp(end))
                                            x0 = reply';
                                        else 
                                            error('invalid response');
                                        end
                                        
                                        disp(sprintf('starting search and fitting...'));
                                        tic;
                                        switch fit_method
                                            case 'lsqnonlin'
                                                lb = [c1-1e-8 c2(1) kp(1)]';
                                                ub = [c1+1e-8 c2(end) kp(end)]';
                                                [x,resnorm,residual,exitflag] = lsqnonlin(@tangent_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            case 'levmar'
                                            case 'simplex'
                                                lb = [c1 c2(1) kp(1)]';
                                                ub = [c1 c2(end) kp(end)]';
                                                A = [0 0 -1];
                                                b = -(1e-8);
                                                nonlcon = [];
                                                [x,fval,exitflag,foutput] = fminsearchcon(@tangent_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
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
                                    case 'function' % c1 = 'located', c2 = 'range','unknown'
                                        reply = input(sprintf('enter 3 parameter starting vector as [%4.3f %4.3f<=c2<=%4.3f c]:(enter for default)\n',[c1;c2(1);c2(end)]));
                                        if isempty(reply)
                                            x0 = [c1 c2(1)+(c2(end)-c2(1)).*rand(1) 0]';
                                        elseif isvector(reply) & length(reply) == 3 & (reply(1) == c1) & (reply(2) >= c2(1) & reply(2) <= c2(end))
                                            x0 = reply';
                                        else 
                                            error('invalid response');
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
                                        error('invalid Kp configuration');
                                end
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
                                switch part_coeff.config
                                    case 'constant' % c1 = 'range','unknown', c2 = 'located'
                                        reply = input(sprintf('what is the range of Kp to search? [1/10:10] (if unspecified, press enter)\n'));
                                        if isempty(reply)
                                            kp = [0.1 10];
                                        elseif isvector(reply) & all(reply >= 0.1 & reply <= 10)
                                            kp = reply;
                                        else
                                            error('invalid response');
                                        end
                                        reply = input(sprintf('enter 3 parameter starting vector as [%4.3f<=c1<=%4.3f %4.3f %4.3f<=kp<=%4.3f]:(enter for default)\n',[c1(1);c1(end);c2;kp(1);kp(end)]));
                                        if isempty(reply)
                                            x0 = [c1(1)+(c1(end)-c1(1)).*rand(1) c2 kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                        elseif isvector(reply) & length(reply) == 3 & (reply(1) >= c1(1) & reply(1) <= c1(end)) & (reply(2) == c2) & (reply(3) >= kp(1) & reply(3) <= kp(end))
                                            x0 = reply';
                                        else 
                                            error('invalid response');
                                        end
                                        
                                        disp(sprintf('starting search and fitting...'));
                                        tic;
                                        switch fit_method
                                            case 'lsqnonlin'
                                                lb = [c1(1) c2-1e-8 kp(1)]';
                                                ub = [c1(end) c2+1e-8 kp(end)]';
                                                [x,resnorm,residual,exitflag] = lsqnonlin(@tangent_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            case 'levmar'
                                            case 'simplex'
                                                lb = [c1(1) c2 kp(1)]';
                                                ub = [c1(end) c2 kp(end)]';
                                                A = [0 0 -1];
                                                b = -(1e-8);
                                                nonlcon = [];
                                                [x,fval,exitflag,foutput] = fminsearchcon(@tangent_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
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
                                    case 'function' % c1 = 'range','unknown', c2 = 'located'
                                        reply = input(sprintf('enter 3 parameter starting vector as [%4.3f<=c1<=%4.3f %4.3f c]:(enter for default)\n',[c1(1);c1(end);c2]));
                                        if isempty(reply)
                                            x0 = [c1(1)+(c1(end)-c1(1)).*rand(1) c2 0]';
                                        elseif isvector(reply) & length(reply) == 3 & (reply(1) >= c1(1) & reply(1) <= c1(end)) & (reply(2) == c2)
                                            x0 = reply';
                                        else 
                                            error('invalid response');
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
                                    otherwise
                                        error('invalid Kp configuration');
                                end 
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
                                switch part_coeff.config
                                    case 'constant' % c1 = 'range','unknown', % c2 = 'range','unknown'
                                        reply = input(sprintf('what is the range of Kp to search? [1/10:10] (if unspecified, press enter)\n'));
                                        if isempty(reply)
                                            kp = [0.1 10];
                                        elseif isvector(reply) & all(reply >= 0.1 & reply <= 10)
                                            kp = reply;
                                        else
                                            error('invalid response');
                                        end
                                        reply = input(sprintf('enter 3 parameter starting vector as [%4.3f<=c1<=%4.3f %4.3f<=c2<=%4.3f %4.3f<=kp<=%4.3f]:(enter for default)\n',[c1(1);c1(end);c2(1);c2(end);kp(1);kp(end)]));
                                        if isempty(reply)
                                            x0 = [c1(1)+(c1(end)-c1(1))*rand(1) c2(1)+(c2(end)-c2(1))*rand(1) kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                        elseif isvector(reply) & length(reply) == 3 & (reply(1) >= c1(1) & reply(1) <= c1(end)) & (reply(2) >= c2(1) & reply(2) <= c2(end)) & (reply(3) >= kp(1) & reply(3) <= kp(end))
                                            x0 = reply';
                                        else 
                                            error('invalid response');
                                        end
                                        
                                        disp(sprintf('starting search and fitting...'));
                                        tic;
                                        switch fit_method
                                            case 'lsqnonlin'
                                                lb = [c1(1) c2(1) kp(1)]';
                                                ub = [c1(end) c2(end) kp(end)]';
                                                [x,resnorm,residual,exitflag] = lsqnonlin(@tangent_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            case 'levmar'
                                            case 'simplex'
                                                lb = [c1(1) c2(1) kp(1)]';
                                                ub = [c1(end) c2(end) kp(end)]';
                                                A = [0 0 -1];
                                                b = -(1e-8);
                                                nonlcon = [];
                                                [x,fval,exitflag,foutput] = fminsearchcon(@tangent_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
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
                                    case 'function' % c1 = 'range','unknown', % c2 = 'range','unknown'
                                        reply = input(sprintf('enter 3 parameter starting vector as [%4.3f<=c1<=%4.3f %4.3f<=c2<=%4.3f c]:(enter for default)\n',[c1(1);c1(end);c2(1);c2(end)]));
                                        if isempty(reply)
                                            x0 = [c1(1)+(c1(end)-c1(1))*rand(1) c2(1)+(c2(end)-c2(1))*rand(1) 0]';
                                        elseif isvector(reply) & length(reply) == 3 & (reply(1) >= c1(1) & reply(1) <= c1(end)) & (reply(2) >= c2(1) & reply(2) <= c2(end))
                                            x0 = reply';
                                        else 
                                            error('invalid response');
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
                                        error('invalid Kp configuration');
                                end
                            otherwise
                                error('invalid critical point 2 position on boundary');
                        end 
                    otherwise
                        error('invalid critical point 1 position on boundary');
                end
            case 'ruled'
                % boundary has two search parameters b1 and b2 for the
                % position of critical point 1 and critical point 2 on the
                % boundary, respectively.
                % loop over number of critical points and get their
                % coordinates.
                % chord length parameterization of boundary.
                disp(sprintf('\nparameterizing boundary...\n'));
                [bdy_length,b] = chord_length(C_bdy); % b = boundary parameter
                boundary.parameterization = b;
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
                
                c1 = boundary.critpt1parameter;
                c2 = boundary.critpt2parameter;
                
                switch fit_method
                    case 'lsqnonlin'
                        options = optimset('OutputFcn',@outfun,'MaxFunEvals',10000,'MaxIter',1000,'TolFun',1e-6,'TolX',1e-6,'LargeScale','on','Display','iter','Jacobian','off');
                    case 'levmar'
                        options = optimset('OutputFcn',@outfun,'MaxFunEvals',10000,'MaxIter',1000,'TolFun',1e-6,'TolX',1e-6,'LargeScale','on','LevenbergMarquardt','on','Display','iter','Jacobian','off');
                    case 'simplex'
                        options = optimset('OutputFcn',@outfun,'MaxFunEvals',10000,'MaxIter',1000,'TolFun',1e-6,'TolX',1e-6,'Display','iter');
                    otherwise
                        error('invalid fit method');
                end
                
                switch boundary.critpt1position
                    case 'located' 
                        switch boundary.critpt2position
                            case 'located'
                                switch part_coeff.config
                                    case 'constant' % c1 = 'located', c2 = 'located'
                                        reply = input(sprintf('what is the range of Kp to search? [1/10:10] (if unspecified, press enter)\n'));
                                        if isempty(reply)
                                            kp = [0.1 10];
                                        elseif isvector(reply) & all(reply >= 0.1 & reply <= 10)
                                            kp = reply;
                                        else
                                            error('invalid response');
                                        end
                                        reply = input(sprintf('enter 4 parameter starting vector as [%4.3f %4.3f a>0 %4.3f<=kp<=%4.3f]:(enter for default)\n',[c1;c2;kp(1);kp(end)]));
                                        if isempty(reply)
                                            x0 = [c1 c2 1 kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                        elseif isvector(reply) & length(reply) == 4 & (reply(1) == c1) & (reply(2) == c2) & (reply(3) > 0) & (reply(4) >= kp(1) & reply(4) <= kp(end))
                                            x0 = reply';
                                        else 
                                            error('invalid response');
                                        end
                                        
                                        disp(sprintf('starting search and fitting...'));
                                        tic;
                                        switch fit_method
                                            case 'lsqnonlin'
                                                lb = [c1-1e-8 c2-1e-8 1e-10 kp(1)]';
                                                ub = [c1+1e-8 c2+1e-8 Inf kp(end)]';
                                                [x,resnorm,residual,exitflag] = lsqnonlin(@ruled_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            case 'levmar'
                                            case 'simplex'
                                                lb = [c1 c2 1e-10 kp(1)]';
                                                ub = [c1 c2 100 kp(end)]';
                                                A = [0 0 0 -1];
                                                b = -(1e-8);
                                                nonlcon = [];
                                                [x,fval,exitflag,foutput] = fminsearchcon(@ruled_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
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
                                    case 'function' % c1 = 'located', c2 = 'located'
                                        reply = input(sprintf('enter 4 parameter starting vector as [%4.3f %4.3f a>0 c]:(enter for default)\n',[c1;c2]));
                                        if isempty(reply)
                                            x0 = [c1 c2 1 0]';
                                        elseif isvector(reply) & length(reply) == 4 & (reply(1) == c1) & (reply(2) == c2) & (reply(3) > 0)
                                            x0 = reply';
                                        else 
                                            error('invalid response');
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
                                        plot_tielines(output);
                                        display = [fitmeasure;x];
                                        disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f]\n',display));
                                        disp(sprintf('done!'));
                                        toc;
                                        return % fitting finished
                                    otherwise
                                        error('invalid Kp configuration');
                                end      
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
                                switch part_coeff.config
                                    case 'constant' % c1 = 'located', c2 = 'range'/'unknown'
                                        reply = input(sprintf('what is the range of Kp to search? [1/10:10] (if unspecified, press enter)\n'));
                                        if isempty(reply)
                                            kp = [0.1 10];
                                        elseif isvector(reply) & all(reply >= 0.1 & reply <= 10)
                                            kp = reply;
                                        else
                                            error('invalid response');
                                        end
                                        reply = input(sprintf('enter 4 parameter starting vector as [%4.3f %4.3f<=c2<=%4.3f a>0 %4.3f<=kp<=%4.3f]:(enter for default)\n',[c1;c2(1);c2(end);kp(1);kp(end)]));
                                        if isempty(reply)
                                            x0 = [c1 c2(1)+(c2(end)-c2(1)).*rand(1) 1 kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                        elseif isvector(reply) & length(reply) == 4 & (reply(1) == c1) & (reply(2) >= c2(1) & reply(2) <= c2(end)) & (reply(3) > 0) & (reply(4) >= kp(1) & reply(4) <= kp(end))
                                            x0 = reply';
                                        else 
                                            error('invalid response');
                                        end
                                        
                                        disp(sprintf('starting search and fitting...'));
                                        tic;
                                        switch fit_method
                                            case 'lsqnonlin'
                                                lb = [c1-1e-8 c2(1) 1e-10 kp(1)]';
                                                ub = [c1+1e-8 c2(end) Inf kp(end)]';
                                                [x,resnorm,residual,exitflag] = lsqnonlin(@ruled_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            case 'levmar'
                                            case 'simplex'
                                                lb = [c1 c2(1) 1e-10 kp(1)]';
                                                ub = [c1 c2(end) 100 kp(end)]';
                                                A = [0 0 0 -1];
                                                b = -(1e-8);
                                                nonlcon = [];
                                                [x,fval,exitflag,foutput] = fminsearchcon(@ruled_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
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
                                    case 'function' % c1 = 'located', c2 = 'range'/'unknown'
                                        reply = input(sprintf('enter 4 parameter starting vector as [%4.3f %4.3f<=c2<=%4.3f a>0 c]:(enter for default)\n',[c1;c2(1);c2(end)]));
                                        if isempty(reply)
                                            x0 = [c1 c2(1)+(c2(end)-c2(1)).*rand(1) 1 0]';
                                        elseif isvector(reply) & length(reply) == 4 & (reply(1) == c1) & (reply(2) >= c2(1) & reply(2) <= c2(end)) & (reply(3) > 0)
                                            x0 = reply';
                                        else 
                                            error('invalid response');
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
                                        error('invalid Kp configuration');
                                end
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
                                switch part_coeff.config
                                    case 'constant' % c1 = 'range','unknown', c2 = 'located'
                                        reply = input(sprintf('what is the range of Kp to search? [1/10:10] (if unspecified, press enter)\n'));
                                        if isempty(reply)
                                            kp = [0.1 10];
                                        elseif isvector(reply) & all(reply >= 0.1 & reply <= 10)
                                            kp = reply;
                                        else
                                            error('invalid response');
                                        end
                                        reply = input(sprintf('enter 4 parameter starting vector as [%4.3f<=c1<=%4.3f %4.3f a>0 %4.3f<=kp<=%4.3f]:(enter for default)\n',[c1(1);c1(end);c2;kp(1);kp(end)]));
                                        if isempty(reply)
                                            x0 = [c1(1)+(c1(end)-c1(1)).*rand(1) c2 1 kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                        elseif isvector(reply) & length(reply) == 4 & (reply(1) >= c1(1) & reply(1) <= c1(end)) & (reply(2) == c2) & (reply(3) > 0) & (reply(4) >= kp(1) & reply(4) <= kp(end))
                                            x0 = reply';
                                        else 
                                            error('invalid response');
                                        end
                                        
                                        disp(sprintf('starting search and fitting...'));
                                        tic;
                                        switch fit_method
                                            case 'lsqnonlin'
                                                lb = [c1(1) c2-1e-8 1e-10 kp(1)]';
                                                ub = [c1(end) c2+1e-8 Inf kp(end)]';
                                                [x,resnorm,residual,exitflag] = lsqnonlin(@ruled_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            case 'levmar'
                                            case 'simplex'
                                                lb = [c1(1) c2 1e-10 kp(1)]';
                                                ub = [c1(end) c2 100 kp(end)]';
                                                A = [0 0 0 -1];
                                                b = -(1e-8);
                                                nonlcon = [];
                                                [x,fval,exitflag,foutput] = fminsearchcon(@ruled_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
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
                                    case 'function' % c1 = 'range','unknown', c2 = 'located'
                                        reply = input(sprintf('enter 4 parameter starting vector as [%4.3f<=c1<=%4.3f %4.3f a>0 c]:(enter for default)\n',[c1(1);c1(end);c2]));
                                        if isempty(reply)
                                            x0 = [c1(1)+(c1(end)-c1(1)).*rand(1) c2 1 0]';
                                        elseif isvector(reply) & length(reply) == 4 & (reply(1) >= c1(1) & reply(1) <= c1(end)) & (reply(2) == c2) & (reply(3) > 0)
                                            x0 = reply';
                                        else 
                                            error('invalid response');
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
                                    otherwise
                                        error('invalid Kp configuration');
                                end
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
                                switch part_coeff.config
                                    case 'constant' % c1 = 'range','unknown', % c2 = 'range','unknown'
                                        reply = input(sprintf('what is the range of Kp to search? [1/10:10] (if unspecified, press enter)\n'));
                                        if isempty(reply)
                                            kp = [0.1 10];
                                        elseif isvector(reply) & all(reply >= 0.1 & reply <= 10)
                                            kp = reply;
                                        else
                                            error('invalid response');
                                        end
                                        reply = input(sprintf('enter 4 parameter starting vector as [%4.3f<=c1<=%4.3f %4.3f<=c2<=%4.3f a>0 %4.3f<=kp<=%4.3f]:(enter for default)\n',[c1(1);c1(end);c2(1);c2(end);kp(1);kp(end)]));
                                        if isempty(reply)
                                            x0 = [c1(1)+(c1(end)-c1(1)).*rand(1) c2(1)+(c2(end)-c2(1)).*rand(1) 1 kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                        elseif isvector(reply) & length(reply) == 4 & (reply(1) >= c1(1) & reply(1) <= c1(end)) & (reply(2) >= c2(1) & reply(2) <= c2(end)) & (reply(3) > 0) & (reply(4) >= kp(1) & reply(4) <= kp(end))
                                            x0 = reply';
                                        else 
                                            error('invalid response');
                                        end
                                        
                                        disp(sprintf('starting search and fitting...'));
                                        tic;
                                        switch fit_method
                                            case 'lsqnonlin'
                                                lb = [c1(1) c2(1) 1e-10 kp(1)]';
                                                ub = [c1(end) c2(end) Inf kp(end)]';
                                                [x,resnorm,residual,exitflag] = lsqnonlin(@ruled_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            case 'levmar'
                                            case 'simplex'
                                                lb = [c1(1) c2(1) 1e-10 kp(1)]';
                                                ub = [c1(end) c2(end) 100 kp(end)]';
                                                A = [0 0 0 -1];
                                                b = -(1e-8);
                                                nonlcon = [];
                                                [x,fval,exitflag,foutput] = fminsearchcon(@ruled_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
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
                                    case 'function' % c1 = 'range','unknown', % c2 = 'range','unknown'
                                        reply = input(sprintf('enter 4 parameter starting vector as [%4.3f<=c1<=%4.3f %4.3f<=c2<=%4.3f a>0 c]:(enter for default)\n',[c1(1);c1(end);c2(1);c2(end)]));
                                        if isempty(reply)
                                            x0 = [c1(1)+(c1(end)-c1(1)).*rand(1) c2(1)+(c2(end)-c2(1)).*rand(1) 1 0]';
                                        elseif isvector(reply) & length(reply) == 4 & (reply(1) >= c1(1) & reply(1) <= c1(end)) & (reply(2) >= c2(1) & reply(2) <= c2(end)) & (reply(3) > 0)
                                            x0 = reply';
                                        else 
                                            error('invalid response');
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
                                        error('invalid Kp configuration');
                                end
                            otherwise
                                error('invalid critical point 2 position on boundary');
                        end 
                    otherwise
                        error('invalid critical point 1 position on boundary');
                end
            otherwise
                error('invalid tieline configuration');
        end
    case '1critpt/1endtl'
        % initial chord length parameterization of boundary.
        disp(sprintf('\nparameterizing boundary...\n'));
        [bdy_length,b] = chord_length(C_bdy); % b = boundary parameter
        boundary.parameterization = b;
        
        disp(sprintf('with one critical point, there must be one end tieline with two end points:\n'));
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
        
        cp = boundary.critptparameter;
        e1 = boundary.endpt1parameter;
        e2 = boundary.endpt2parameter;
        
        switch fit_method
            case 'lsqnonlin'
                options = optimset('OutputFcn',@outfun,'MaxFunEvals',10000,'MaxIter',10000,'TolFun',1e-6,'TolX',1e-6,'LargeScale','on','Display','iter','Jacobian','off');
            case 'levmar'
                options = optimset('OutputFcn',@outfun,'MaxFunEvals',10000,'MaxIter',10000,'TolFun',1e-6,'TolX',1e-6,'LargeScale','on','LevenbergMarquardt','on','Display','iter','Jacobian','off');
            case 'simplex'
                options = optimset('OutputFcn',@outfun,'MaxFunEvals',10000,'MaxIter',10000,'TolFun',1e-6,'TolX',1e-6,'Display','iter');
            otherwise
                error('invalid fit method');
        end
        
        switch tieline.config
            case 'parallel'
                switch boundary.endpt1position
                    case 'located'
                        switch boundary.endpt2position
                            case 'located'
                                switch part_coeff.config
                                    case 'constant' % e1 = 'located', e2 = 'located'
                                        reply = input(sprintf('what is the range of Kp to search? [1/10:10] (if unspecified, press enter)\n'));
                                        if isempty(reply)
                                            kp = [0.1 10];
                                        elseif isvector(reply) & all(reply >= 0.1 & reply <= 10)
                                            kp = reply;
                                        else
                                            error('invalid response');
                                        end
                                        reply = input(sprintf('enter 3 parameter starting vector as [%4.3f %4.3f %4.3f<=kp<=%4.3f]:(enter for default)\n',[e1;e2;kp(1);kp(end)]));
                                        if isempty(reply)
                                            x0 = [e1 e2 kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                        elseif isvector(reply) & length(reply) == 3 & (reply(1) == e1) & (reply(2) == e2) & (reply(3) >= kp(1) & reply(3) <= kp(end))
                                            x0 = reply';
                                        else 
                                            error('invalid response');
                                        end
                                        
                                        disp(sprintf('starting search and fitting...'));
                                        switch fit_method
                                            case 'lsqnonlin'
                                                lb = [e1-1e-8 e2-1e-8 kp(1)]';
                                                ub = [e1+1e-8 e2+1e-8 kp(end)]';
                                                [x,resnorm,residual,exitflag] = lsqnonlin(@parallel_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            case 'levmar'
                                            case 'simplex'
                                                lb = [e1 e2 kp(1)]';
                                                ub = [e1 e2 kp(end)]';
                                                A = [0 0 -1];
                                                b = -(1e-8);
                                                nonlcon = [];
                                                [x,fval,exitflag,foutput] = fminsearchcon(@parallel_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            otherwise
                                                error('invalid fit method');
                                        end

                                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                                        disp(sprintf('exitflag = %d',exitflag));
                                        [fitmeasure,output] = parallel_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                        plot_tielines(output);
                                        display = [fitmeasure;x];
                                        disp(sprintf('best fit measure = %f, with parameters [%f %f %f]\n',display));
                                        disp(sprintf('done!'));
                                        return % fitting finished
                                    case 'function' % e1 = 'located', e2 = 'located'
                                        reply = input(sprintf('enter 4 parameter starting vector as [%4.3f %4.3f c1 c2]:(enter for default)\n',[e1;e2]));
                                        if isempty(reply)
                                            x0 = [e1 e2 0 0]';
                                        elseif isvector(reply) & length(reply) == 4 & (reply(1) == e1) & (reply(2) == e2)
                                            x0 = reply';
                                        else 
                                            error('invalid response');
                                        end 
                                        
                                        disp(sprintf('starting search and fitting...'));
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
                                        plot_tielines(output);
                                        display = [fitmeasure;x];
                                        disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f]\n',display));
                                        disp(sprintf('done!'));
                                        return % fitting finished
                                    otherwise
                                        error('invalid Kp configuration');
                                end
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
                                switch part_coeff.config
                                    case 'constant' % e1 = 'located', e2 = 'range','unknown'
                                        reply = input(sprintf('what is the range of Kp to search? [1/10:10] (if unspecified, press enter)\n'));
                                        if isempty(reply)
                                            kp = [0.1 10];
                                        elseif isvector(reply) & all(reply >= 0.1 & reply <= 10)
                                            kp = reply;
                                        else
                                            error('invalid response');
                                        end
                                        reply = input(sprintf('enter 3 parameter starting vector as [%4.3f %4.3f<=e2<=%4.3f %4.3f<=kp<=%4.3f]:(enter for default)\n',[e1;e2(1);e2(end);kp(1);kp(end)]));
                                        if isempty(reply)
                                            x0 = [e1 e2(1)+(e2(end)-e2(1)).*rand(1) kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                        elseif isvector(reply) & length(reply) == 3 & (reply(1) == e1) & (reply(2) >= e2(1) & reply(2) <= e2(end)) & (reply(3) >= kp(1) & reply(3) <= kp(end))
                                            x0 = reply';
                                        else 
                                            error('invalid response');
                                        end
                                        
                                        disp(sprintf('starting search and fitting...'));
                                        switch fit_method
                                            case 'lsqnonlin'
                                                lb = [e1-1e-8 e2(1) kp(1)]';
                                                ub = [e1+1e-8 e2(end) kp(end)]';
                                                [x,resnorm,residual,exitflag] = lsqnonlin(@parallel_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            case 'levmar'
                                            case 'simplex'
                                                lb = [e1 e2(1) kp(1)]';
                                                ub = [e1 e2(end) kp(end)]';
                                                A = [0 0 -1];
                                                b = -(1e-8);
                                                nonlcon = [];
                                                [x,fval,exitflag,foutput] = fminsearchcon(@parallel_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            otherwise
                                                error('invalid fit method');
                                        end

                                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                                        disp(sprintf('exitflag = %d',exitflag));
                                        [fitmeasure,output] = parallel_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                        plot_tielines(output);
                                        display = [fitmeasure;x];
                                        disp(sprintf('best fit measure = %f, with parameters [%f %f %f]\n',display));
                                        disp(sprintf('done!'));
                                        return % fitting finished
                                    case 'function' % e1 = 'located', e2 = 'range','unknown'
                                        reply = input(sprintf('enter 4 parameter starting vector as [%4.3f %4.3f<=e2<=%4.3f c1 c2]:(enter for default)\n',[e1;e2(1);e2(end)]));
                                        if isempty(reply)
                                            x0 = [e1 e2(1)+(e2(end)-e2(1)).*rand(1) 0 0]';
                                        elseif isvector(reply) & length(reply) == 4 & (reply(1) == e1) & (reply(2) >= e2(1) & reply(2) <= e2(end))
                                            x0 = reply';
                                        else 
                                            error('invalid response');
                                        end 
                                        
                                        disp(sprintf('starting search and fitting...'));
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
                                        return % fitting finished
                                    otherwise
                                        error('invalid Kp configuration');
                                end
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
                                switch part_coeff.config
                                    case 'constant' % e1 = 'range','unknown', e2 = 'located'
                                        reply = input(sprintf('what is the range of Kp to search? [1/10:10] (if unspecified, press enter)\n'));
                                        if isempty(reply)
                                            kp = [0.1 10];
                                        elseif isvector(reply) & all(reply >= 0.1 & reply <= 10)
                                            kp = reply;
                                        else
                                            error('invalid response');
                                        end
                                        reply = input(sprintf('enter 3 parameter starting vector as [%4.3f<=e1<=%4.3f %4.3f %4.3f<=kp<=%4.3f]:(enter for default)\n',[e1(1);e1(end);e2;kp(1);kp(end)]));
                                        if isempty(reply)
                                            x0 = [e1(1)+(e1(end)-e1(1)).*rand(1) e2 kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                        elseif isvector(reply) & length(reply) == 3 & (reply(1) >= e1(1) & reply(1) <= e1(end)) & (reply(2) == e2) & (reply(3) >= kp(1) & reply(3) <= kp(end))
                                            x0 = reply';
                                        else 
                                            error('invalid response');
                                        end
                                        
                                        disp(sprintf('starting search and fitting...'));
                                        switch fit_method
                                            case 'lsqnonlin'
                                                lb = [e1(1) e2-1e-8 kp(1)]';
                                                ub = [e1(end) e2+1e-8 kp(end)]';
                                                [x,resnorm,residual,exitflag] = lsqnonlin(@parallel_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            case 'levmar'
                                            case 'simplex'
                                                lb = [e1(1) e2 kp(1)]';
                                                ub = [e1(end) e2 kp(end)]';
                                                A = [0 0 -1];
                                                b = -(1e-8);
                                                nonlcon = [];
                                                [x,fval,exitflag,foutput] = fminsearchcon(@parallel_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            otherwise
                                                error('invalid fit method');
                                        end

                                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                                        disp(sprintf('exitflag = %d',exitflag));
                                        [fitmeasure,output] = parallel_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                        plot_tielines(output);
                                        display = [fitmeasure;x];
                                        disp(sprintf('best fit measure = %f, with parameters [%f %f %f]\n',display));
                                        disp(sprintf('done!'));
                                        return % fitting finished
                                    case 'function' % e1 = 'range','unknown', e2 = 'located'
                                        reply = input(sprintf('enter 4 parameter starting vector as [%4.3f<=e1<=%4.3f %4.3f c1 c2]:(enter for default)\n',[e1(1);e1(end);e2]));
                                        if isempty(reply)
                                            x0 = [e1(1)+(e1(end)-e1(1)).*rand(1) e2 0 0]';
                                        elseif isvector(reply) & length(reply) == 4 & (reply(1) >= e1(1) & reply(1) <= e1(end)) & (reply(2) == e2)
                                            x0 = reply';
                                        else 
                                            error('invalid response');
                                        end 
                                        
                                        disp(sprintf('starting search and fitting...'));
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
                                        return % fitting finished
                                    otherwise
                                        error('invalid Kp configuration');
                                end
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
                                switch part_coeff.config
                                    case 'constant' % e1 = 'range','unknown', % e2 = 'range','unknown'
                                        reply = input(sprintf('what is the range of Kp to search? [1/10:10] (if unspecified, press enter)\n'));
                                        if isempty(reply)
                                            kp = [0.1 10];
                                        elseif isvector(reply) & all(reply >= 0.1 & reply <= 10)
                                            kp = reply;
                                        else
                                            error('invalid response');
                                        end
                                        reply = input(sprintf('enter 3 parameter starting vector as [%4.3f<=e1<=%4.3f %4.3f<=e2<=%4.3f %4.3f<=kp<=%4.3f]:(enter for default)\n',[e1(1);e1(end);e2(1);e2(end);kp(1);kp(end)]));
                                        if isempty(reply)
                                            x0 = [e1(1)+(e1(end)-e1(1)).*rand(1) e2(1)+(e2(end)-e2(1)).*rand(1) kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                        elseif isvector(reply) & length(reply) == 3 & (reply(1) >= e1(1) & reply(1) <= e1(end)) & (reply(2) >= e2(1) & reply(2) <= e2(end)) & (reply(3) >= kp(1) & reply(3) <= kp(end))
                                            x0 = reply';
                                        else 
                                            error('invalid response');
                                        end
                                        
                                        disp(sprintf('starting search and fitting...'));
                                        switch fit_method
                                            case 'lsqnonlin'
                                                lb = [e1(1) e2(1) kp(1)]';
                                                ub = [e1(end) e2(end) kp(end)]';
                                                [x,resnorm,residual,exitflag] = lsqnonlin(@parallel_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            case 'levmar'
                                            case 'simplex'
                                                lb = [e1(1) e2(1) kp(1)]';
                                                ub = [e1(end) e2(end) kp(end)]';
                                                A = [0 0 -1];
                                                b = -(1e-8);
                                                nonlcon = [];
                                                [x,fval,exitflag,foutput] = fminsearchcon(@parallel_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                            otherwise
                                                error('invalid fit method');
                                        end

                                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                                        disp(sprintf('exitflag = %d',exitflag));
                                        [fitmeasure,output] = parallel_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                        plot_tielines(output);
                                        display = [fitmeasure;x];
                                        disp(sprintf('best fit measure = %f, with parameters [%f %f %f]\n',display));
                                        disp(sprintf('done!'));
                                        return % fitting finished
                                    case 'function' % e1 = 'range','unknown', % e2 = 'range','unknown'
                                        reply = input(sprintf('enter 4 parameter starting vector as [%4.3f<=e1<=%4.3f %4.3f<=e2<=%4.3f c1 c2]:(enter for default)\n',[e1(1);e1(end);e2(1);e2(end)]));
                                        if isempty(reply)
                                            x0 = [e1(1)+(e1(end)-e1(1)).*rand(1) e2(1)+(e2(end)-e2(1)).*rand(1) 0 0]';
                                        elseif isvector(reply) & length(reply) == 4 & (reply(1) >= e1(1) & reply(1) <= e1(end)) & (reply(2) >= e2(1) & reply(2) <= e2(end))
                                            x0 = reply';
                                        else 
                                            error('invalid response');
                                        end 
                                        
                                        disp(sprintf('starting search and fitting...'));
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
                                        return % fitting finished
                                    otherwise
                                        error('invalid Kp configuration');
                                end
                            otherwise
                                error('invalid end point 2 position');
                        end
                    otherwise
                        error('invalid end point 1 position');
                end
            case 'tangent'
                switch boundary.critptposition
                    case 'located'
                        switch boundary.endpt1position
                            case 'located'
                                switch boundary.endpt2position
                                    case 'located'
                                        switch part_coeff.config
                                            case 'constant' % cp = 'located', e1 = 'located', e2 = 'located'
                                                reply = input(sprintf('what is the range of Kp to search? [1/10:10] (if unspecified, press enter)\n'));
                                                if isempty(reply)
                                                    kp = [0.1 10];
                                                elseif isvector(reply) & all(reply >= 0.1 & reply <= 10)
                                                    kp = reply;
                                                else
                                                    error('invalid response');
                                                end
                                                reply = input(sprintf('enter 4 parameter starting vector as [%4.3f %4.3f %4.3f %4.3f<=kp<=%4.3f]:(enter for default)\n',[cp;e1;e2;kp(1);kp(end)]));
                                                if isempty(reply)
                                                    x0 = [cp e1 e2 1 kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                                elseif isvector(reply) & length(reply) == 4 & (reply(1) == cp) & (reply(2) == e1) & (reply(3) == e2) & (reply(4) >= kp(1) & reply(4) <= kp(end))
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end
                                                
                                                disp(sprintf('starting search and fitting...'));
                                                switch fit_method
                                                    case 'lsqnonlin'
                                                        lb = [cp-1e-8 e1-1e-8 e2-1e-8 kp(1)]';
                                                        ub = [cp+1e-8 e1+1e-8 e2+1e-8 kp(end)]';
                                                        [x,resnorm,residual,exitflag] = lsqnonlin(@tangent_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                    case 'levmar'
                                                    case 'simplex'
                                                        lb = [cp e1 e2 kp(1)]';
                                                        ub = [cp e1 e2 kp(end)]';
                                                        A = [0 0 0 -1];
                                                        b = -(1e-8);
                                                        nonlcon = [];
                                                        [x,fval,exitflag,foutput] = fminsearchcon(@tangent_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                    otherwise
                                                        error('invalid fit method');
                                                end
                                                
                                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                                disp(sprintf('exitflag = %d',exitflag));
                                                [fitmeasure,output] = tangent_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                plot_tielines(output);
                                                display = [fitmeasure;x];
                                                disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f]\n',display));
                                                disp(sprintf('done!'));
                                                return % fitting finished
                                            case 'function' % cp = 'located', e1 = 'located', e2 = 'located'
                                                reply = input(sprintf('enter 5 parameter starting vector as [%4.3f %4.3f %4.3f c1 c2]:(enter for default)\n',[cp;e1;e2]));
                                                if isempty(reply)
                                                    x0 = [cp e1 e2 0 0]';
                                                elseif isvector(reply) & length(reply) == 5 & (reply(1) == cp) & (reply(2) == e1) & (reply(3) == e2) 
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end 
                                                
                                                disp(sprintf('starting search and fitting...'));
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
                                                plot_tielines(output);
                                                display = [fitmeasure;x];
                                                disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f %f]\n',display));
                                                disp(sprintf('done!'));
                                                return % fitting finished
                                            otherwise
                                                error('invalid Kp configuration');
                                        end
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
                                        switch part_coeff.config
                                            case 'constant' % cp = 'located', e1 = 'located', e2 = 'range'/'unknown'
                                                reply = input(sprintf('what is the range of Kp to search? [1/10:10] (if unspecified, press enter)\n'));
                                                if isempty(reply)
                                                    kp = [0.1 10];
                                                elseif isvector(reply) & all(reply >= 0.1 & reply <= 10)
                                                    kp = reply;
                                                else
                                                    error('invalid response');
                                                end
                                                reply = input(sprintf('enter 4 parameter starting vector as [%4.3f %4.3f %4.3f<=e2<=%4.3f %4.3f<=kp<=%4.3f]:(enter for default)\n',[cp;e1;e2(1);e2(end);kp(1);kp(end)]));
                                                if isempty(reply)
                                                    x0 = [cp e1 e2(1)+(e2(end)-e2(1)).*rand(1) kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                                elseif isvector(reply) & length(reply) == 4 & (reply(1) == cp) & (reply(2) == e1) & (reply(3) >= e2(1) & reply(3) <= e2(end)) & (reply(4) >= kp(1) & reply(4) <= kp(end))
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end
                                                
                                                disp(sprintf('starting search and fitting...'));
                                                switch fit_method
                                                    case 'lsqnonlin'
                                                        lb = [cp-1e-8 e1-1e-8 e2(1) kp(1)]';
                                                        ub = [cp+1e-8 e1+1e-8 e2(end) kp(end)]';
                                                        [x,resnorm,residual,exitflag] = lsqnonlin(@tangent_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                    case 'levmar'
                                                    case 'simplex'
                                                        lb = [cp e1 e2(1) kp(1)]';
                                                        ub = [cp e1 e2(end) kp(end)]';
                                                        A = [0 0 0 -1];
                                                        b = -(1e-8);
                                                        nonlcon = [];
                                                        [x,fval,exitflag,foutput] = fminsearchcon(@tangent_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                    otherwise
                                                        error('invalid fit method');
                                                end
                                                
                                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                                disp(sprintf('exitflag = %d',exitflag));
                                                [fitmeasure,output] = tangent_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                plot_tielines(output);
                                                display = [fitmeasure;x];
                                                disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f]\n',display));
                                                disp(sprintf('done!'));
                                                return % fitting finished
                                            case 'function' % cp = 'located', e1 = 'located', e2 = 'range'/'unknown'
                                                reply = input(sprintf('enter 5 parameter starting vector as [%4.3f %4.3f %4.3f<=e2<=%4.3f c1 c2]:(enter for default)\n',[cp;e1;e2(1);e2(end)]));
                                                if isempty(reply)
                                                    x0 = [cp e1 e2(1)+(e2(end)-e2(1)).*rand(1) 0 0]';
                                                elseif isvector(reply) & length(reply) == 5 & (reply(1) == cp) & (reply(2) == e1) & (reply(3) >= e2(1) & reply(3) <= e2(end))
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end 
                                                
                                                disp(sprintf('starting search and fitting...'));
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
                                                return % fitting finished
                                            otherwise
                                                error('invalid Kp configuration');
                                        end
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
                                        switch part_coeff.config
                                            case 'constant' % cp = 'located', e2 = 'located', e1 = 'range'/'unknown'
                                                reply = input(sprintf('what is the range of Kp to search? [1/10:10] (if unspecified, press enter)\n'));
                                                if isempty(reply)
                                                    kp = [0.1 10];
                                                elseif isvector(reply) & all(reply >= 0.1 & reply <= 10)
                                                    kp = reply;
                                                else
                                                    error('invalid response');
                                                end
                                                reply = input(sprintf('enter 4 parameter starting vector as [%4.3f %4.3f<=e1<=%4.3f %4.3f %4.3f<=kp<=%4.3f]:(enter for default)\n',[cp;e1(1);e1(end);e2;kp(1);kp(end)]));
                                                if isempty(reply)
                                                    x0 = [cp e1(1)+(e1(end)-e1(1)).*rand(1) e2 kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                                elseif isvector(reply) & length(reply) == 5 & (reply(1) == cp) & (reply(2) >= e1(1) & reply(2) <= e1(end)) & (reply(3) == e2) & (reply(4) >= kp(1) & reply(4) <= kp(end))
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end
                                                
                                                disp(sprintf('starting search and fitting...'));
                                                switch fit_method
                                                    case 'lsqnonlin'
                                                        lb = [cp-1e-8 e1(1) e2-1e-8 kp(1)]';
                                                        ub = [cp+1e-8 e1(end) e2+1e-8 kp(end)]';
                                                        [x,resnorm,residual,exitflag] = lsqnonlin(@tangent_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                    case 'levmar'
                                                    case 'simplex'
                                                        lb = [cp e1(1) e2 kp(1)]';
                                                        ub = [cp e1(end) e2 kp(end)]';
                                                        A = [0 0 0 -1];
                                                        b = -(1e-8);
                                                        nonlcon = [];
                                                        [x,fval,exitflag,foutput] = fminsearchcon(@tangent_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                    otherwise
                                                        error('invalid fit method');
                                                end
                                                
                                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                                disp(sprintf('exitflag = %d',exitflag));
                                                [fitmeasure,output] = tangent_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                plot_tielines(output);
                                                display = [fitmeasure;x];
                                                disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f]\n',display));
                                                disp(sprintf('done!'));
                                                return % fitting finished
                                            case 'function' % cp = 'located', e2 = 'located', e1 = 'range'/'unknown'
                                                reply = input(sprintf('enter 5 parameter starting vector as [%4.3f %4.3f<=e1<=%4.3f %4.3f c1 c2]:(enter for default)\n',[cp;e1(1);e1(end);e2]));
                                                if isempty(reply)
                                                    x0 = [cp e1(1)+(e1(end)-e1(1)).*rand(1) e2 0 0]';
                                                elseif isvector(reply) & length(reply) == 5 & (reply(1) == cp) & (reply(2) >= e1(1) & reply(2) <= e1(end)) & (reply(3) == e2) 
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end 
                                                
                                                disp(sprintf('starting search and fitting...'));
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
                                                return % fitting finished
                                            otherwise
                                                error('invalid Kp configuration');
                                        end
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
                                        switch part_coeff.config
                                            case 'constant' % cp = 'located', e2 = 'range'/'unknown', e1 = 'range'/'unknown'
                                                reply = input(sprintf('what is the range of Kp to search? [1/10:10] (if unspecified, press enter)\n'));
                                                if isempty(reply)
                                                    kp = [0.1 10];
                                                elseif isvector(reply) & all(reply >= 0.1 & reply <= 10)
                                                    kp = reply;
                                                else
                                                    error('invalid response');
                                                end
                                                reply = input(sprintf('enter 4 parameter starting vector as [%4.3f %4.3f<=e1<=%4.3f %4.3f<=e2<=%4.3f %4.3f<=kp<=%4.3f]:(enter for default)\n',[cp;e1(1);e1(end);e2(1);e2(end);kp(1);kp(end)]));
                                                if isempty(reply)
                                                    x0 = [cp e1(1)+(e1(end)-e1(1)).*rand(1) e2(1)+(e2(end)-e2(1)).*rand(1) kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                                elseif isvector(reply) & length(reply) == 4 & (reply(1) == cp) & (reply(2) >= e1(1) & reply(2) <= e1(end)) & (reply(3) >= e2(1) & reply(3) <= e2(end)) & (reply(4) >= kp(1) & reply(4) <= kp(end))
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end
                                                
                                                disp(sprintf('starting search and fitting...'));
                                                switch fit_method
                                                    case 'lsqnonlin'
                                                        lb = [cp-1e-8 e1(1) e2(1) kp(1)]';
                                                        ub = [cp+1e-8 e1(end) e2(end) kp(end)]';
                                                        [x,resnorm,residual,exitflag] = lsqnonlin(@tangent_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                    case 'levmar'
                                                    case 'simplex'
                                                        lb = [cp e1(1) e2(1) kp(1)]';
                                                        ub = [cp e1(end) e2(end) kp(end)]';
                                                        A = [0 0 0 -1];
                                                        b = -(1e-8);
                                                        nonlcon = [];
                                                        [x,fval,exitflag,foutput] = fminsearchcon(@tangent_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                    otherwise
                                                        error('invalid fit method');
                                                end
                                                
                                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                                disp(sprintf('exitflag = %d',exitflag));
                                                [fitmeasure,output] = tangent_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                plot_tielines(output);
                                                display = [fitmeasure;x];
                                                disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f]\n',display));
                                                disp(sprintf('done!'));
                                                return % fitting finished
                                            case 'function' % cp = 'located', e2 = 'range'/'unknown', e1 = 'range'/'unknown'
                                                reply = input(sprintf('enter 5 parameter starting vector as [%4.3f %4.3f<=e1<=%4.3f %4.3f<=e2<=%4.3f c1 c2]:(enter for default)\n',[cp;e1(1);e1(end);e2(1);e2(end)]));
                                                if isempty(reply)
                                                    x0 = [cp e1(1)+(e1(end)-e1(1)).*rand(1) e2 0 0]';
                                                elseif isvector(reply) & length(reply) == 5 & (reply(1) == cp) & (reply(2) >= e1(1) & reply(2) <= e1(end)) & (reply(3) >= e2(1) & reply(3) <= e2(end))
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end 
                                                
                                                disp(sprintf('starting search and fitting...'));
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
                                                return % fitting finished
                                            otherwise
                                                error('invalid Kp configuration');
                                        end
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
                                        switch part_coeff.config
                                            case 'constant' % cp = 'range'/'unknown', e1 = 'located', e2 = 'located'
                                                reply = input(sprintf('what is the range of Kp to search? [1/10:10] (if unspecified, press enter)\n'));
                                                if isempty(reply)
                                                    kp = [0.1 10];
                                                elseif isvector(reply) & all(reply >= 0.1 & reply <= 10)
                                                    kp = reply;
                                                else
                                                    error('invalid response');
                                                end
                                                reply = input(sprintf('enter 4 parameter starting vector as [%4.3f<=cp<=%4.3f %4.3f %4.3f %4.3f<=kp<=%4.3f]:(enter for default)\n',[cp(1);cp(end);e1;e2;kp(1);kp(end)]));
                                                if isempty(reply)
                                                    x0 = [cp(1)+(cp(end)-cp(1)).*rand(1) e1 e2 kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                                elseif isvector(reply) & length(reply) == 4 & (reply(1) >= cp(1) & reply(1) <= cp(end)) & (reply(2) == e1) & (reply(3) == e2) & (reply(4) >= kp(1) & reply(4) <= kp(end))
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end
                                                
                                                disp(sprintf('starting search and fitting...'));
                                                switch fit_method
                                                    case 'lsqnonlin'
                                                        lb = [cp(1) e1-1e-8 e2-1e-8 kp(1)]';
                                                        ub = [cp(end) e1+1e-8 e2+1e-8 kp(end)]';
                                                        [x,resnorm,residual,exitflag] = lsqnonlin(@tangent_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                    case 'levmar'
                                                    case 'simplex'
                                                        lb = [cp(1) e1 e2 kp(1)]';
                                                        ub = [cp(end) e1 e2 kp(end)]';
                                                        A = [0 0 0 -1];
                                                        b = -(1e-8);
                                                        nonlcon = [];
                                                        [x,fval,exitflag,foutput] = fminsearchcon(@tangent_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                    otherwise
                                                        error('invalid fit method');
                                                end
                                                
                                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                                disp(sprintf('exitflag = %d',exitflag));
                                                [fitmeasure,output] = tangent_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                plot_tielines(output);
                                                display = [fitmeasure;x];
                                                disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f]\n',display));
                                                disp(sprintf('done!'));
                                                return % fitting finished
                                            case 'function' % cp = 'range'/'unknown', e1 = 'located', e2 = 'located'
                                                reply = input(sprintf('enter 5 parameter starting vector as [%4.3f<=cp<=%4.3f %4.3f %4.3f c1 c2]:(enter for default)\n',[cp(1);cp(end);e1;e2]));
                                                if isempty(reply)
                                                    x0 = [cp(1)+(cp(end)-cp(1)).*rand(1) e1 e2 0 0]';
                                                elseif isvector(reply) & length(reply) == 5 & (reply(1) >= cp(1) & reply(1) <= cp(end)) & (reply(2) == e1) & (reply(3) == e2)
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end 
                                                
                                                disp(sprintf('starting search and fitting...'));
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
                                                return % fitting finished
                                            otherwise
                                                error('invalid Kp configuration');
                                        end
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
                                        switch part_coeff.config
                                            case 'constant' % cp = 'range'/'unknown', e1 = 'located', e2 = 'range'/'unknown'
                                                reply = input(sprintf('what is the range of Kp to search? [1/10:10] (if unspecified, press enter)\n'));
                                                if isempty(reply)
                                                    kp = [0.1 10];
                                                elseif isvector(reply) & all(reply >= 0.1 & reply <= 10)
                                                    kp = reply;
                                                else
                                                    error('invalid response');
                                                end
                                                reply = input(sprintf('enter 4 parameter starting vector as [%4.3f<=cp<=%4.3f %4.3f %4.3f<=e2<=%4.3f %4.3f<=kp<=%4.3f]:(enter for default)\n',[cp(1);cp(end);e1;e2(1);e2(end);kp(1);kp(end)]));
                                                if isempty(reply)
                                                    x0 = [cp(1)+(cp(end)-cp(1)).*rand(1) e1 e2(1)+(e2(end)-e2(1)).*rand(1) kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                                elseif isvector(reply) & length(reply) == 4 & (reply(1) >= cp(1) & reply(1) <= cp(end)) & (reply(2) == e1) & (reply(3) >= e2(1) & reply(3) <= e2(end)) & (reply(4) >= kp(1) & reply(4) <= kp(end))
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end
                                                
                                                disp(sprintf('starting search and fitting...'));
                                                switch fit_method
                                                    case 'lsqnonlin'
                                                        lb = [cp(1) e1-1e-8 e2(1) kp(1)]';
                                                        ub = [cp(end) e1+1e-8 e2(end) kp(end)]';
                                                        [x,resnorm,residual,exitflag] = lsqnonlin(@tangent_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                    case 'levmar'
                                                    case 'simplex'
                                                        lb = [cp(1) e1 e2(1) kp(1)]';
                                                        ub = [cp(end) e1 e2(end) kp(end)]';
                                                        A = [0 0 0 -1];
                                                        b = -(1e-8);
                                                        nonlcon = [];
                                                        [x,fval,exitflag,foutput] = fminsearchcon(@tangent_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                    otherwise
                                                        error('invalid fit method');
                                                end
                                                
                                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                                disp(sprintf('exitflag = %d',exitflag));
                                                [fitmeasure,output] = tangent_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                plot_tielines(output);
                                                display = [fitmeasure;x];
                                                disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f]\n',display));
                                                disp(sprintf('done!'));
                                                return % fitting finished
                                            case 'function' % cp = 'range'/'unknown', e1 = 'located', e2 = 'range'/'unknown'
                                                reply = input(sprintf('enter 5 parameter starting vector as [%4.3f<=cp<=%4.3f %4.3f %4.3f<=e2<=%4.3f a>0 c1 c2]:(enter for default)\n',[cp(1);cp(end);e1;e2(1);e2(end)]));
                                                if isempty(reply)
                                                    x0 = [cp(1)+(cp(end)-cp(1)).*rand(1) e1 e2(1)+(e2(end)-e2(1)).*rand(1) 0 0]';
                                                elseif isvector(reply) & length(reply) == 5 & (reply(1) >= cp(1) & reply(1) <= cp(end)) & (reply(2) == e1) & (reply(3) >= e2(1) & reply(3) <= e2(end))
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end 
                                                
                                                disp(sprintf('starting search and fitting...'));
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
                                                return % fitting finished
                                            otherwise
                                                error('invalid Kp configuration');
                                        end
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
                                        switch part_coeff.config
                                            case 'constant' % cp = 'range'/'unknown', e2 = 'located', e1 = 'range'/'unknown'
                                                reply = input(sprintf('what is the range of Kp to search? [1/10:10] (if unspecified, press enter)\n'));
                                                if isempty(reply)
                                                    kp = [0.1 10];
                                                elseif isvector(reply) & all(reply >= 0.1 & reply <= 10)
                                                    kp = reply;
                                                else
                                                    error('invalid response');
                                                end
                                                reply = input(sprintf('enter 4 parameter starting vector as [%4.3f<=cp<=%4.3f %4.3f<=e1<=%4.3f %4.3f %4.3f<=kp<=%4.3f]:(enter for default)\n',[cp(1);cp(end);e1(1);e1(end);e2;kp(1);kp(end)]));
                                                if isempty(reply)
                                                    x0 = [cp(1)+(cp(end)-cp(1)).*rand(1) e1(1)+(e1(end)-e1(1)).*rand(1) e2 kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                                elseif isvector(reply) & length(reply) == 4 & (reply(1) >= cp(1) & reply(1) <= cp(end)) & (reply(2) >= e1(1) & reply(2) <= e1(end)) & (reply(3) == e2) & (reply(4) >= kp(1) & reply(4) <= kp(end))
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end
                                                
                                                disp(sprintf('starting search and fitting...'));
                                                switch fit_method
                                                    case 'lsqnonlin'
                                                        lb = [cp(1) e1(1) e2-1e-8 kp(1)]';
                                                        ub = [cp(end) e1(end) e2+1e-8 kp(end)]';
                                                        [x,resnorm,residual,exitflag] = lsqnonlin(@tangent_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                    case 'levmar'
                                                    case 'simplex'
                                                        lb = [cp(1) e1(1) e2 kp(1)]';
                                                        ub = [cp(end) e1(end) e2 kp(end)]';
                                                        A = [0 0 0 -1];
                                                        b = -(1e-8);
                                                        nonlcon = [];
                                                        [x,fval,exitflag,foutput] = fminsearchcon(@tangent_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                    otherwise
                                                        error('invalid fit method');
                                                end
                                                
                                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                                disp(sprintf('exitflag = %d',exitflag));
                                                [fitmeasure,output] = tangent_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                plot_tielines(output);
                                                display = [fitmeasure;x];
                                                disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f]\n',display));
                                                disp(sprintf('done!'));
                                                return % fitting finished
                                            case 'function' % cp = 'range'/'unknown', e2 = 'located', e1 = 'range'/'unknown'
                                                reply = input(sprintf('enter 5 parameter starting vector as [%4.3f<=cp<=%4.3f %4.3f<=e1<=%4.3f %4.3f c1 c2]:(enter for default)\n',[cp(1);cp(end);e1(1);e1(end);e2]));
                                                if isempty(reply)
                                                    x0 = [cp(1)+(cp(end)-cp(1)).*rand(1) e1(1)+(e1(end)-e1(1)).*rand(1) e2 0 0]';
                                                elseif isvector(reply) & length(reply) == 5 & (reply(1) >= cp(1) & reply(1) <= cp(end)) & (reply(2) >= e1(1) & reply(2) <= e1(end)) & (reply(3) == e2) 
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end 
                                                
                                                disp(sprintf('starting search and fitting...'));
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
                                                return % fitting finished
                                            otherwise
                                                error('invalid Kp configuration');
                                        end
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
                                        switch part_coeff.config
                                            case 'constant' % cp = 'range'/'unknown', e2 = 'range'/'unknown', e1 = 'range'/'unknown'
                                                reply = input(sprintf('what is the range of Kp to search? [1/10:10] (if unspecified, press enter)\n'));
                                                if isempty(reply)
                                                    kp = [0.1 10];
                                                elseif isvector(reply) & all(reply >= 0.1 & reply <= 10)
                                                    kp = reply;
                                                else
                                                    error('invalid response');
                                                end
                                                reply = input(sprintf('enter 4 parameter starting vector as [%4.3f<=cp<=%4.3f %4.3f<=e1<=%4.3f %4.3f<=e2<=%4.3f %4.3f<=kp<=%4.3f]:(enter for default)\n',[cp(1);cp(end);e1(1);e1(end);e2(1);e2(end);kp(1);kp(end)]));
                                                if isempty(reply)
                                                    x0 = [cp(1)+(cp(end)-cp(1)).*rand(1) e1(1)+(e1(end)-e1(1)).*rand(1) e2(1)+(e2(end)-e2(1)).*rand(1) kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                                elseif isvector(reply) & length(reply) == 4 & (reply(1) >= cp(1) & reply(1) <= cp(end)) & (reply(2) >= e1(1) & reply(2) <= e1(end)) & (reply(3) >= e2(1) & reply(3) <= e2(end)) & (reply(4) >= kp(1) & reply(4) <= kp(end))
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end
                                                
                                                disp(sprintf('starting search and fitting...'));
                                                switch fit_method
                                                    case 'lsqnonlin'
                                                        lb = [cp(1) e1(1) e2(1) kp(1)]';
                                                        ub = [cp(end) e1(end) e2(end) kp(end)]';
                                                        [x,resnorm,residual,exitflag] = lsqnonlin(@tangent_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                    case 'levmar'
                                                    case 'simplex'
                                                        lb = [cp(1) e1(1) e2(1) kp(1)]';
                                                        ub = [cp(end) e1(end) e2(end) kp(end)]';
                                                        A = [0 0 0 -1];
                                                        b = -(1e-8);
                                                        nonlcon = [];
                                                        [x,fval,exitflag,foutput] = fminsearchcon(@tangent_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                    otherwise
                                                        error('invalid fit method');
                                                end
                                                
                                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                                disp(sprintf('exitflag = %d',exitflag));
                                                [fitmeasure,output] = tangent_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                plot_tielines(output);
                                                display = [fitmeasure;x];
                                                disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f]\n',display));
                                                disp(sprintf('done!'));
                                                return % fitting finished
                                            case 'function' % cp = 'range'/'unknown', e2 = 'range'/'unknown', e1 = 'range'/'unknown'
                                                reply = input(sprintf('enter 5 parameter starting vector as [%4.3f<=cp<=%4.3f %4.3f<=e1<=%4.3f %4.3f<=e2<=%4.3f c1 c2]:(enter for default)\n',[cp(1);cp(end);e1(1);e1(end);e2(1);e2(end)]));
                                                if isempty(reply)
                                                    x0 = [cp(1)+(cp(end)-cp(1)).*rand(1) e1(1)+(e1(end)-e1(1)).*rand(1) e2(1)+(e2(end)-e2(1)).*rand(1) 0 0]';
                                                elseif isvector(reply) & length(reply) == 5 & (reply(1) >= cp(1) & reply(1) <= cp(end)) & (reply(2) >= e1(1) & reply(2) <= e1(end)) & (reply(3) >= e2(1) & reply(3) <= e2(end)) 
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end 
                                                
                                                disp(sprintf('starting search and fitting...'));
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
                                                return % fitting finished
                                            otherwise
                                                error('invalid Kp configuration');
                                        end
                                    otherwise
                                        error('invalid end point 2 position');
                                end
                            otherwise
                                error('invalid end point 1 position');
                        end
                    otherwise
                        error('invalid critical point position');
                end
            case 'ruled'
                switch boundary.critptposition
                    case 'located'
                        switch boundary.endpt1position
                            case 'located'
                                switch boundary.endpt2position
                                    case 'located'
                                        switch part_coeff.config
                                            case 'constant' % cp = 'located', e1 = 'located', e2 = 'located'
                                                reply = input(sprintf('what is the range of Kp to search? [1/10:10] (if unspecified, press enter)\n'));
                                                if isempty(reply)
                                                    kp = [0.1 10];
                                                elseif isvector(reply) & all(reply >= 0.1 & reply <= 10)
                                                    kp = reply;
                                                else
                                                    error('invalid response');
                                                end
                                                reply = input(sprintf('enter 5 parameter starting vector as [%4.3f %4.3f %4.3f a>0 %4.3f<=kp<=%4.3f]:(enter for default)\n',[cp;e1;e2;kp(1);kp(end)]));
                                                if isempty(reply)
                                                    x0 = [cp e1 e2 1 kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                                elseif isvector(reply) & length(reply) == 5 & (reply(1) == cp) & (reply(2) == e1) & (reply(3) == e2) & (reply(4) > 0) & (reply(5) >= kp(1) & reply(5) <= kp(end))
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end
                                                
                                                disp(sprintf('starting search and fitting...'));
                                                switch fit_method
                                                    case 'lsqnonlin'
                                                        lb = [cp-1e-8 e1-1e-8 e2-1e-8 1e-10 kp(1)]';
                                                        ub = [cp+1e-8 e1+1e-8 e2+1e-8 Inf kp(end)]';
                                                        [x,resnorm,residual,exitflag] = lsqnonlin(@ruled_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                    case 'levmar'
                                                    case 'simplex'
                                                        lb = [cp e1 e2 1e-10 kp(1)]';
                                                        ub = [cp e1 e2 100 kp(end)]';
                                                        A = [0 0 0 0 -1];
                                                        b = -(1e-8);
                                                        nonlcon = [];
                                                        [x,fval,exitflag,foutput] = fminsearchcon(@ruled_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                    otherwise
                                                        error('invalid fit method');
                                                end
                                                
                                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                                disp(sprintf('exitflag = %d',exitflag));
                                                [fitmeasure,output] = ruled_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                plot_tielines(output);
                                                display = [fitmeasure;x];
                                                disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f %f]\n',display));
                                                disp(sprintf('done!'));
                                                return % fitting finished
                                            case 'function' % cp = 'located', e1 = 'located', e2 = 'located'
                                                reply = input(sprintf('enter 6 parameter starting vector as [%4.3f %4.3f %4.3f a>0 c1 c2]:(enter for default)\n',[cp;e1;e2]));
                                                if isempty(reply)
                                                    x0 = [cp e1 e2 1 0 0]';
                                                elseif isvector(reply) & length(reply) == 6 & (reply(1) == cp) & (reply(2) == e1) & (reply(3) == e2) & (reply(4) > 0)
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
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
                                                plot_tielines(output);
                                                display = [fitmeasure;x];
                                                disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f %f %f]\n',display));
                                                disp(sprintf('done!'));
                                                toc;
                                                return % fitting finished
                                            otherwise
                                                error('invalid Kp configuration');
                                        end
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
                                        switch part_coeff.config
                                            case 'constant' % cp = 'located', e1 = 'located', e2 = 'range'/'unknown'
                                                reply = input(sprintf('what is the range of Kp to search? [1/10:10] (if unspecified, press enter)\n'));
                                                if isempty(reply)
                                                    kp = [0.1 10];
                                                elseif isvector(reply) & all(reply >= 0.1 & reply <= 10)
                                                    kp = reply;
                                                else
                                                    error('invalid response');
                                                end
                                                reply = input(sprintf('enter 5 parameter starting vector as [%4.3f %4.3f %4.3f<=e2<=%4.3f a>0 %4.3f<=kp<=%4.3f]:(enter for default)\n',[cp;e1;e2(1);e2(end);kp(1);kp(end)]));
                                                if isempty(reply)
                                                    x0 = [cp e1 e2(1)+(e2(end)-e2(1)).*rand(1) 1 kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                                elseif isvector(reply) & length(reply) == 5 & (reply(1) == cp) & (reply(2) == e1) & (reply(3) >= e2(1) & reply(3) <= e2(end)) & (reply(4) > 0) & (reply(5) >= kp(1) & reply(5) <= kp(end))
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end
                                                
                                                disp(sprintf('starting search and fitting...'));
                                                switch fit_method
                                                    case 'lsqnonlin'
                                                        lb = [cp-1e-8 e1-1e-8 e2(1) 1e-10 kp(1)]';
                                                        ub = [cp+1e-8 e1+1e-8 e2(end) Inf kp(end)]';
                                                        [x,resnorm,residual,exitflag] = lsqnonlin(@ruled_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                    case 'levmar'
                                                    case 'simplex'
                                                        lb = [cp e1 e2(1) 1e-10 kp(1)]';
                                                        ub = [cp e1 e2(end) 100 kp(end)]';
                                                        A = [0 0 0 0 -1];
                                                        b = -(1e-8);
                                                        nonlcon = [];
                                                        [x,fval,exitflag,foutput] = fminsearchcon(@ruled_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                    otherwise
                                                        error('invalid fit method');
                                                end
                                                
                                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                                disp(sprintf('exitflag = %d',exitflag));
                                                [fitmeasure,output] = ruled_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                plot_tielines(output);
                                                display = [fitmeasure;x];
                                                disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f %f]\n',display));
                                                disp(sprintf('done!'));
                                                return % fitting finished
                                            case 'function' % cp = 'located', e1 = 'located', e2 = 'range'/'unknown'
                                                reply = input(sprintf('enter 6 parameter starting vector as [%4.3f %4.3f %4.3f<=e2<=%4.3f a>0 c1 c2]:(enter for default)\n',[cp;e1;e2(1);e2(end)]));
                                                if isempty(reply)
                                                    x0 = [cp e1 e2(1)+(e2(end)-e2(1)).*rand(1) 1 0 0]';
                                                elseif isvector(reply) & length(reply) == 6 & (reply(1) == cp) & (reply(2) == e1) & (reply(3) >= e2(1) & reply(3) <= e2(end)) & (reply(4) > 0)
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end
                                                
                                                disp(sprintf('starting search and fitting...'));
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
                                                return % fitting finished
                                            otherwise
                                                error('invalid Kp configuration');
                                        end
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
                                        switch part_coeff.config
                                            case 'constant' % cp = 'located', e2 = 'located', e1 = 'range'/'unknown'
                                                reply = input(sprintf('what is the range of Kp to search? [1/10:10] (if unspecified, press enter)\n'));
                                                if isempty(reply)
                                                    kp = [0.1 10];
                                                elseif isvector(reply) & all(reply >= 0.1 & reply <= 10)
                                                    kp = reply;
                                                else
                                                    error('invalid response');
                                                end
                                                reply = input(sprintf('enter 5 parameter starting vector as [%4.3f %4.3f<=e1<=%4.3f %4.3f a>0 %4.3f<=kp<=%4.3f]:(enter for default)\n',[cp;e1(1);e1(end);e2;kp(1);kp(end)]));
                                                if isempty(reply)
                                                    x0 = [cp e1(1)+(e1(end)-e1(1)).*rand(1) e2 1 kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                                elseif isvector(reply) & length(reply) == 5 & (reply(1) == cp) & (reply(2) >= e1(1) & reply(2) <= e1(end)) & (reply(3) == e2) & (reply(4) > 0) & (reply(5) >= kp(1) & reply(5) <= kp(end))
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end
                                                
                                                disp(sprintf('starting search and fitting...'));
                                                switch fit_method
                                                    case 'lsqnonlin'
                                                        lb = [cp-1e-8 e1(1) e2-1e-8 1e-10 kp(1)]';
                                                        ub = [cp+1e-8 e1(end) e2+1e-8 Inf kp(end)]';
                                                        [x,resnorm,residual,exitflag] = lsqnonlin(@ruled_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                    case 'levmar'
                                                    case 'simplex'
                                                        lb = [cp e1(1) e2 1e-10 kp(1)]';
                                                        ub = [cp e1(end) e2 100 kp(end)]';
                                                        A = [0 0 0 0 -1];
                                                        b = -(1e-8);
                                                        nonlcon = [];
                                                        [x,fval,exitflag,foutput] = fminsearchcon(@ruled_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                    otherwise
                                                        error('invalid fit method');
                                                end
                                                
                                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                                disp(sprintf('exitflag = %d',exitflag));
                                                [fitmeasure,output] = ruled_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                plot_tielines(output);
                                                display = [fitmeasure;x];
                                                disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f %f]\n',display));
                                                disp(sprintf('done!'));
                                                return % fitting finished
                                            case 'function' % cp = 'located', e2 = 'located', e1 = 'range'/'unknown'
                                                reply = input(sprintf('enter 6 parameter starting vector as [%4.3f %4.3f<=e1<=%4.3f %4.3f a>0 c1 c2]:(enter for default)\n',[cp;e1(1);e1(end);e2]));
                                                if isempty(reply)
                                                    x0 = [cp e1(1)+(e1(end)-e1(1)).*rand(1) e2 1 0 0]';
                                                elseif isvector(reply) & length(reply) == 6 & (reply(1) == cp) & (reply(2) >= e1(1) & reply(2) <= e1(end)) & (reply(3) == e2) & (reply(4) > 0)
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end 
                                                
                                                disp(sprintf('starting search and fitting...'));
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
                                                return % fitting finished
                                            otherwise
                                                error('invalid Kp configuration');
                                        end
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
                                        switch part_coeff.config
                                            case 'constant' % cp = 'located', e2 = 'range'/'unknown', e1 = 'range'/'unknown'
                                                reply = input(sprintf('what is the range of Kp to search? [1/10:10] (if unspecified, press enter)\n'));
                                                if isempty(reply)
                                                    kp = [0.1 10];
                                                elseif isvector(reply) & all(reply >= 0.1 & reply <= 10)
                                                    kp = reply;
                                                else
                                                    error('invalid response');
                                                end
                                                reply = input(sprintf('enter 5 parameter starting vector as [%4.3f %4.3f<=e1<=%4.3f %4.3f<=e2<=%4.3f a>0 %4.3f<=kp<=%4.3f]:(enter for default)\n',[cp;e1(1);e1(end);e2(1);e2(end);kp(1);kp(end)]));
                                                if isempty(reply)
                                                    x0 = [cp e1(1)+(e1(end)-e1(1)).*rand(1) e2(1)+(e2(end)-e2(1)).*rand(1) 1 kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                                elseif isvector(reply) & length(reply) == 5 & (reply(1) == cp) & (reply(2) >= e1(1) & reply(2) <= e1(end)) & (reply(3) >= e2(1) & reply(3) <= e2(end)) & (reply(4) > 0) & (reply(5) >= kp(1) & reply(5) <= kp(end))
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end
                                                
                                                disp(sprintf('starting search and fitting...'));
                                                switch fit_method
                                                    case 'lsqnonlin'
                                                        lb = [cp-1e-8 e1(1) e2(1) 1e-10 kp(1)]';
                                                        ub = [cp+1e-8 e1(end) e2(end) Inf kp(end)]';
                                                        [x,resnorm,residual,exitflag] = lsqnonlin(@ruled_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                    case 'levmar'
                                                    case 'simplex'
                                                        lb = [cp e1(1) e2(1) 1e-10 kp(1)]';
                                                        ub = [cp e1(end) e2(end) 100 kp(end)]';
                                                        A = [0 0 0 0 -1];
                                                        b = -(1e-8);
                                                        nonlcon = [];
                                                        [x,fval,exitflag,foutput] = fminsearchcon(@ruled_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                    otherwise
                                                        error('invalid fit method');
                                                end
                                                
                                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                                disp(sprintf('exitflag = %d',exitflag));
                                                [fitmeasure,output] = ruled_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                plot_tielines(output);
                                                display = [fitmeasure;x];
                                                disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f %f]\n',display));
                                                disp(sprintf('done!'));
                                                return % fitting finished
                                            case 'function' % cp = 'located', e2 = 'range'/'unknown', e1 = 'range'/'unknown'
                                                reply = input(sprintf('enter 6 parameter starting vector as [%4.3f %4.3f<=e1<=%4.3f %4.3f<=e2<=%4.3f a>0 c1 c2]:(enter for default)\n',[cp;e1(1);e1(end);e2(1);e2(end)]));
                                                if isempty(reply)
                                                    x0 = [cp e1(1)+(e1(end)-e1(1)).*rand(1) e2(1)+(e2(end)-e2(1)).*rand(1) 1 0 0]';
                                                elseif isvector(reply) & length(reply) == 6 & (reply(1) == cp) & (reply(2) >= e1(1) & reply(2) <= e1(end)) & (reply(3) >= e2(1) & reply(3) <= e2(end)) & (reply(4) > 0)
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end 
                                                
                                                disp(sprintf('starting search and fitting...'));
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
                                                return % fitting finished
                                            otherwise
                                                error('invalid Kp configuration');
                                        end
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
                                        switch part_coeff.config
                                            case 'constant' % cp = 'range'/'unknown', e1 = 'located', e2 = 'located'
                                                reply = input(sprintf('what is the range of Kp to search? [1/10:10] (if unspecified, press enter)\n'));
                                                if isempty(reply)
                                                    kp = [0.1 10];
                                                elseif isvector(reply) & all(reply >= 0.1 & reply <= 10)
                                                    kp = reply;
                                                else
                                                    error('invalid response');
                                                end
                                                reply = input(sprintf('enter 5 parameter starting vector as [%4.3f<=cp<=%4.3f %4.3f %4.3f a>0 %4.3f<=kp<=%4.3f]:(enter for default)\n',[cp(1);cp(end);e1;e2;kp(1);kp(end)]));
                                                if isempty(reply)
                                                    x0 = [cp(1)+(cp(end)-cp(1)).*rand(1) e1 e2 1 kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                                elseif isvector(reply) & length(reply) == 5 & (reply(1) >= cp(1) & reply(1) <= cp(end)) & (reply(2) == e1) & (reply(3) == e2) & (reply(4) > 0) & (reply(5) >= kp(1) & reply(5) <= kp(end))
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end
                                                
                                                disp(sprintf('starting search and fitting...'));
                                                switch fit_method
                                                    case 'lsqnonlin'
                                                        lb = [cp(1) e1-1e-8 e2-1e-8 1e-10 kp(1)]';
                                                        ub = [cp(end) e1+1e-8 e2+1e-8 Inf kp(end)]';
                                                        [x,resnorm,residual,exitflag] = lsqnonlin(@ruled_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                    case 'levmar'
                                                    case 'simplex'
                                                        lb = [cp(1) e1 e2 1e-10 kp(1)]';
                                                        ub = [cp(end) e1 e2 100 kp(end)]';
                                                        A = [0 0 0 0 -1];
                                                        b = -(1e-8);
                                                        nonlcon = [];
                                                        [x,fval,exitflag,foutput] = fminsearchcon(@ruled_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                    otherwise
                                                        error('invalid fit method');
                                                end
                                                
                                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                                disp(sprintf('exitflag = %d',exitflag));
                                                [fitmeasure,output] = ruled_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                plot_tielines(output);
                                                display = [fitmeasure;x];
                                                disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f %f]\n',display));
                                                disp(sprintf('done!'));
                                                return % fitting finished
                                            case 'function' % cp = 'range'/'unknown', e1 = 'located', e2 = 'located'
                                                reply = input(sprintf('enter 6 parameter starting vector as [%4.3f<=cp<=%4.3f %4.3f %4.3f a>0 c1 c2]:(enter for default)\n',[cp(1);cp(end);e1;e2]));
                                                if isempty(reply)
                                                    x0 = [cp(1)+(cp(end)-cp(1)).*rand(1) e1 e2 1 0 0]';
                                                elseif isvector(reply) & length(reply) == 6 & (reply(1) >= cp(1) & reply(1) <= cp(end)) & (reply(2) == e1) & (reply(3) == e2) & (reply(4) > 0)
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end 
                                                
                                                disp(sprintf('starting search and fitting...'));
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
                                                return % fitting finished
                                            otherwise
                                                error('invalid Kp configuration');
                                        end
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
                                        switch part_coeff.config
                                            case 'constant' % cp = 'range'/'unknown', e1 = 'located', e2 = 'range'/'unknown'
                                                reply = input(sprintf('what is the range of Kp to search? [1/10:10] (if unspecified, press enter)\n'));
                                                if isempty(reply)
                                                    kp = [0.1 10];
                                                elseif isvector(reply) & all(reply >= 0.1 & reply <= 10)
                                                    kp = reply;
                                                else
                                                    error('invalid response');
                                                end
                                                reply = input(sprintf('enter 5 parameter starting vector as [%4.3f<=cp<=%4.3f %4.3f %4.3f<=e2<=%4.3f a>0 %4.3f<=kp<=%4.3f]:(enter for default)\n',[cp(1);cp(end);e1;e2(1);e2(end);kp(1);kp(end)]));
                                                if isempty(reply)
                                                    x0 = [cp(1)+(cp(end)-cp(1)).*rand(1) e1 e2(1)+(e2(end)-e2(1)).*rand(1) 1 kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                                elseif isvector(reply) & length(reply) == 5 & (reply(1) >= cp(1) & reply(1) <= cp(end)) & (reply(2) == e1) & (reply(3) >= e2(1) & reply(3) <= e2(end)) & (reply(4) > 0) & (reply(5) >= kp(1) & reply(5) <= kp(end))
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end
                                                
                                                disp(sprintf('starting search and fitting...'));
                                                switch fit_method
                                                    case 'lsqnonlin'
                                                        lb = [cp(1) e1-1e-8 e2(1) 1e-10 kp(1)]';
                                                        ub = [cp(end) e1+1e-8 e2(end) Inf kp(end)]';
                                                        [x,resnorm,residual,exitflag] = lsqnonlin(@ruled_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                    case 'levmar'
                                                    case 'simplex'
                                                        lb = [cp(1) e1 e2(1) 1e-10 kp(1)]';
                                                        ub = [cp(end) e1 e2(end) 100 kp(end)]';
                                                        A = [0 0 0 0 -1];
                                                        b = -(1e-8);
                                                        nonlcon = [];
                                                        [x,fval,exitflag,foutput] = fminsearchcon(@ruled_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                    otherwise
                                                        error('invalid fit method');
                                                end
                                                
                                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                                disp(sprintf('exitflag = %d',exitflag));
                                                [fitmeasure,output] = ruled_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                plot_tielines(output);
                                                display = [fitmeasure;x];
                                                disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f %f]\n',display));
                                                disp(sprintf('done!'));
                                                return % fitting finished
                                            case 'function' % cp = 'range'/'unknown', e1 = 'located', e2 = 'range'/'unknown'
                                                reply = input(sprintf('enter 6 parameter starting vector as [%4.3f<=cp<=%4.3f %4.3f %4.3f<=e2<=%4.3f a>0 c1 c2]:(enter for default)\n',[cp(1);cp(end);e1;e2(1);e2(end)]));
                                                if isempty(reply)
                                                    x0 = [cp(1)+(cp(end)-cp(1)).*rand(1) e1 e2(1)+(e2(end)-e2(1)).*rand(1) 1 0 0]';
                                                elseif isvector(reply) & length(reply) == 6 & (reply(1) >= cp(1) & reply(1) <= cp(end)) & (reply(2) == e1) & (reply(3) >= e2(1) & reply(3) <= e2(end)) & (reply(4) > 0)
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end 
                                                
                                                disp(sprintf('starting search and fitting...'));
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
                                                return % fitting finished
                                            otherwise
                                                error('invalid Kp configuration');
                                        end
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
                                        switch part_coeff.config
                                            case 'constant' % cp = 'range'/'unknown', e2 = 'located', e1 = 'range'/'unknown'
                                                reply = input(sprintf('what is the range of Kp to search? [1/10:10] (if unspecified, press enter)\n'));
                                                if isempty(reply)
                                                    kp = [0.1 10];
                                                elseif isvector(reply) & all(reply >= 0.1 & reply <= 10)
                                                    kp = reply;
                                                else
                                                    error('invalid response');
                                                end
                                                reply = input(sprintf('enter 5 parameter starting vector as [%4.3f<=cp<=%4.3f %4.3f<=e1<=%4.3f %4.3f a>0 %4.3f<=kp<=%4.3f]:(enter for default)\n',[cp(1);cp(end);e1(1);e1(end);e2;kp(1);kp(end)]));
                                                if isempty(reply)
                                                    x0 = [cp(1)+(cp(end)-cp(1)).*rand(1) e1(1)+(e1(end)-e1(1)).*rand(1) e2 1 kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                                elseif isvector(reply) & length(reply) == 5 & (reply(1) >= cp(1) & reply(1) <= cp(end)) & (reply(2) >= e1(1) & reply(2) <= e1(end)) & (reply(3) == e2) & (reply(4) > 0) & (reply(5) >= kp(1) & reply(5) <= kp(end))
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end
                                                
                                                disp(sprintf('starting search and fitting...'));
                                                switch fit_method
                                                    case 'lsqnonlin'
                                                        lb = [cp(1) e1(1) e2-1e-8 1e-10 kp(1)]';
                                                        ub = [cp(end) e1(end) e2+1e-8 Inf kp(end)]';
                                                        [x,resnorm,residual,exitflag] = lsqnonlin(@ruled_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                    case 'levmar'
                                                    case 'simplex'
                                                        lb = [cp(1) e1(1) e2 1e-10 kp(1)]';
                                                        ub = [cp(end) e1(end) e2 100 kp(end)]';
                                                        A = [0 0 0 0 -1];
                                                        b = -(1e-8);
                                                        nonlcon = [];
                                                        [x,fval,exitflag,foutput] = fminsearchcon(@ruled_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                    otherwise
                                                        error('invalid fit method');
                                                end
                                                
                                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                                disp(sprintf('exitflag = %d',exitflag));
                                                [fitmeasure,output] = ruled_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                plot_tielines(output);
                                                display = [fitmeasure;x];
                                                disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f %f]\n',display));
                                                disp(sprintf('done!'));
                                                return % fitting finished
                                            case 'function' % cp = 'range'/'unknown', e2 = 'located', e1 = 'range'/'unknown'
                                                reply = input(sprintf('enter 6 parameter starting vector as [%4.3f<=cp<=%4.3f %4.3f<=e1<=%4.3f %4.3f a>0 c1 c2]:(enter for default)\n',[cp(1);cp(end);e1(1);e1(end);e2]));
                                                if isempty(reply)
                                                    x0 = [cp(1)+(cp(end)-cp(1)).*rand(1) e1(1)+(e1(end)-e1(1)).*rand(1) e2 1 0 0]';
                                                elseif isvector(reply) & length(reply) == 6 & (reply(1) >= cp(1) & reply(1) <= cp(end)) & (reply(2) >= e1(1) & reply(2) <= e1(end)) & (reply(3) == e2) & (reply(4) > 0)
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end 
                                                
                                                disp(sprintf('starting search and fitting...'));
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
                                                return % fitting finished
                                            otherwise
                                                error('invalid Kp configuration');
                                        end
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
                                        switch part_coeff.config
                                            case 'constant' % cp = 'range'/'unknown', e2 = 'range'/'unknown', e1 = 'range'/'unknown'
                                                reply = input(sprintf('what is the range of Kp to search? [1/10:10] (if unspecified, press enter)\n'));
                                                if isempty(reply)
                                                    kp = [0.1 10];
                                                elseif isvector(reply) & all(reply >= 0.1 & reply <= 10)
                                                    kp = reply;
                                                else
                                                    error('invalid response');
                                                end
                                                reply = input(sprintf('enter 5 parameter starting vector as [%4.3f<=cp<=%4.3f %4.3f<=e1<=%4.3f %4.3f<=e2<=%4.3f a>0 %4.3f<=kp<=%4.3f]:(enter for default)\n',[cp(1);cp(end);e1(1);e1(end);e2(1);e2(end);kp(1);kp(end)]));
                                                if isempty(reply)
                                                    x0 = [cp(1)+(cp(end)-cp(1)).*rand(1) e1(1)+(e1(end)-e1(1)).*rand(1) e2(1)+(e2(end)-e2(1)).*rand(1) 1 kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                                elseif isvector(reply) & length(reply) == 5 & (reply(1) >= cp(1) & reply(1) <= cp(end)) & (reply(2) >= e1(1) & reply(2) <= e1(end)) & (reply(3) >= e2(1) & reply(3) <= e2(end)) & (reply(4) > 0) & (reply(5) >= kp(1) & reply(5) <= kp(end))
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end
                                                
                                                disp(sprintf('starting search and fitting...'));
                                                switch fit_method
                                                    case 'lsqnonlin'
                                                        lb = [cp(1) e1(1) e2(1) 1e-10 kp(1)]';
                                                        ub = [cp(end) e1(end) e2(end) Inf kp(end)]';
                                                        [x,resnorm,residual,exitflag] = lsqnonlin(@ruled_tieline_fit_fxn,x0,lb,ub,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                    case 'levmar'
                                                    case 'simplex'
                                                        lb = [cp(1) e1(1) e2(1) 1e-10 kp(1)]';
                                                        ub = [cp(end) e1(end) e2(end) 100 kp(end)]';
                                                        A = [0 0 0 0 -1];
                                                        b = -(1e-8);
                                                        nonlcon = [];
                                                        [x,fval,exitflag,foutput] = fminsearchcon(@ruled_tieline_fit_fxn,x0,lb,ub,A,b,nonlcon,options,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                    otherwise
                                                        error('invalid fit method');
                                                end
                                                
                                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                                disp(sprintf('exitflag = %d',exitflag));
                                                [fitmeasure,output] = ruled_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method);
                                                plot_tielines(output);
                                                display = [fitmeasure;x];
                                                disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f %f]\n',display));
                                                disp(sprintf('done!'));
                                                return % fitting finished
                                            case 'function' % cp = 'range'/'unknown', e2 = 'range'/'unknown', e1 = 'range'/'unknown'
                                                reply = input(sprintf('enter 6 parameter starting vector as [%4.3f<=cp<=%4.3f %4.3f<=e1<=%4.3f %4.3f<=e2<=%4.3f a>0 c1 c2]:(enter for default)\n',[cp(1);cp(end);e1(1);e1(end);e2(1);e2(end)]));
                                                if isempty(reply)
                                                    x0 = [cp(1)+(cp(end)-cp(1)).*rand(1) e1(1)+(e1(end)-e1(1)).*rand(1) e2(1)+(e2(end)-e2(1)).*rand(1) 1 0 0]';
                                                elseif isvector(reply) & length(reply) == 6 & (reply(1) >= cp(1) & reply(1) <= cp(end)) & (reply(2) >= e1(1) & reply(2) <= e1(end)) & (reply(3) >= e2(1) & reply(3) <= e2(end)) & (reply(4) > 0)
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
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
                                                plot_tielines(output);
                                                display = [fitmeasure;x];
                                                disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f %f %f]\n',display));
                                                disp(sprintf('done!'));
                                                toc;
                                                return % fitting finished
                                            otherwise
                                                error('invalid Kp configuration');
                                        end
                                    otherwise
                                        error('invalid end point 2 position');
                                end
                            otherwise
                                error('invalid end point 1 position');
                        end
                    otherwise
                        error('invalid critical point position');
                end
            otherwise
                error('invalid tieline configuration');
        end
    case '2endtls'
        switch tieline.config
            case 'parallel'
                % if end tielines parallel, only one configuration;if not,
                % no configuration
            case 'tangent'
                % only one configuration
            case 'ruled'
                
            otherwise
                error('invalid tieline configuration');
        end
    otherwise
        error('invalid boundary configuration');
end
        
return