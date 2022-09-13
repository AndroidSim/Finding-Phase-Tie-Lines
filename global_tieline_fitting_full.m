function [fit_output,tieline_output,spectra_output] = global_tieline_fitting(S_coex,C_coex,S_bdy,C_bdy,nc,np)
% [fit_out,tieline_out] = global_tieline_fitting(S_coex,C_coex,S_bdy,C_bdy,np,nc,varargin)
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

if ~isscalar(nc)
    error('number of chemical components must be a scalar');
end

if ~isscalar(np)
    error('number of phases must be a scalar');
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
% area of the absorbance spectrum equals one, then convert back to the
% derivative spectrum.
disp(sprintf('normalizing...\n'));
allS = normalize_spectra(allS,'d');
S_coex = allS(:,1:ncScol);
S_bdy = allS(:,ncScol+1:ncScol+nbScol);
clear allS;
magnetic_fields = S_coex(:,1); % or S_bdy, magnetic field values the same for all spectra because of alignment
% reduced all spectra to a vector of intensity values (magnetic fields not needed for fitting)
S_coex = S_coex(:,2:2:end);
S_bdy = S_bdy(:,2:2:end);

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
weight_fxn = ones(size(stdbS));
% weight_fxn = stdbS./sum(stdbS);
% w = zeros(size(weight_fxn));
% w(weight_fxn >= 10*mean(std_bSbline)) = 1;
% w(weight_fxn >= 100*mean(std_bSbline)) = 2;
% w(weight_fxn >= 1000*mean(std_bSbline)) = 3;
% w(weight_fxn >= 10000*mean(std_bSbline)) = 4;
% coeff_var = stdbS./abs(mean(ISb')');
% absSb = deriv2abs(Sb);
% IabsSb = absSb(:,2:2:end);
% stdabsbS = std(IabsSb')';
% coeff_var_absSb = stdabsbS./abs(mean(IabsSb')');

% linearly transform ternary(i.e. homogeneous barycentric coordinates) to 
% cartesian coordinates using 1st (Sm) and 3rd (chol) mole fractions
C_bdy = tern2cart(C_bdy,1);
C_coex = tern2cart(C_coex,1);

% determine sizes of inputs before fitting
[ncC,ncd] = size(C_coex);
[ncB,ncS] = size(S_coex);
[nbC,nbd] = size(C_bdy);
[nbB,nbS] = size(S_bdy);

% reshape coexistence spectra matrix for fitting
weighting = repmat(sqrt(weight_fxn),1,ncS);
std_cS = repmat(std_cSbline*sqrt(ncB)*sqrt(ncS),ncB,1);
S_coex = (S_coex./(std_cS)).*weighting;
data = reshape(S_coex,ncB*ncS,1);

% declare spectra data structure that contains information about the spectra
spectra = struct('magnetic_fields',magnetic_fields,'cSbline_var',var_cSbline,'cSbline_std',std_cSbline,'weight_fxn',weight_fxn);
% declare and initialize boundary, tieline, and Kp search and data structures
boundary = struct('config',0,'search_method',0);
tieline = struct('config',0);
part_coeff = struct('config',0);

% ask user for the phase boundary configuration
reply = input(sprintf('"2critpts", "1critpt/1endtl", or "2endtls" boundary configuration?\n'),'s');
while isempty(reply) || ~any(strcmp(reply,{'2critpts','1critpt/1endtl','2endtls'}))
    reply = input(sprintf('"2critpts", "1critpt/1endtl", or "2endtls" boundary configuration?\n'),'s');
end
boundary.config = reply;
% ask user for boundary search options
reply = input(sprintf('"continuous" or "grid" boundary search method?\n'),'s');
while isempty(reply) || ~any(strcmp(reply,{'continuous','grid'}))
    reply = input(sprintf('"continuous" or "grid" boundary search method?\n'),'s');
end
boundary.search_method = reply;
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
                switch boundary.search_method % same as tieline.search_method
                    case 'continuous'
                        reply = input(sprintf('what is the range of slope angles to search? [0:180] (if unspecified, press enter)\n'));
                        disp(sprintf('parameterizing tieline slope...\n'));
                        if isempty(reply)
                            % default is 50 points for grid search
                            % parameter also thought of as angles from 0 to 180
                            x = [0 1];
                        elseif isvector(reply) & all(reply >= 0 & reply <= 180)
                            x = reply.*(1/180);
                        else
                            error('invalid response');
                        end
                        options = optimset('MaxFunEvals',500,'TolFun',1e-8,'TolX',1e-8,'MaxIter',500,'LargeScale','on','LevenbergMarquardt','on','Display','iter','Jacobian','off');
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
                                reply = input(sprintf('enter 2 parameter starting vector as [%4.3f<=slope<=%4.3f %4.2f<=kp<=%4.2f]:(enter for default)\n',[x(1)*180;x(end)*180;kp(1);kp(end)]));
                                if isempty(reply)
                                    x0 = [x(1)+(x(end)-x(1)).*rand(1) kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                elseif isvector(reply) & length(reply) == 2 & all(reply(1).*(1/180) >= x(1) & reply(1).*(1/180) <= x(end)) & all(reply(2) >= kp(1) & reply(2) <= kp(end))
                                    x0 = [reply(1).*(1/180) reply(2)]';
                                else 
                                    error('invalid response');
                                end
                                lb = [x(1) kp(1)]';
                                ub = [x(end) kp(end)]';
                                disp(sprintf('starting search and fitting...'));
                                [x,resnorm] = lsqcurvefit(@parallel_tieline_fit_fxn,x0,C_coex,data,lb,ub,options,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                            case 'function'
                                reply = input(sprintf('enter 2 parameter starting vector as [%4.3f<=slope<=%4.3f c]:(enter for default)\n',[x(1)*180;x(end)*180]));
                                if isempty(reply)
                                    x0 = [x(1)+(x(end)-x(1)).*rand(1) 0]';% for 2 critical points, Kp function has 1 parameter
                                elseif isvector(reply) & length(reply) == 2 & all(reply(1).*(1/180) >= x(1) & reply(1).*(1/180) <= x(end))
                                    x0 = [reply(1).*(1/180) reply(2)]';
                                else 
                                    error('invalid response');
                                end 
                                lb = [x(1) -Inf]';
                                ub = [x(end) Inf]';
                                disp(sprintf('starting search and fitting...'));
                                [x,resnorm] = lsqcurvefit(@parallel_tieline_fit_fxn,x0,C_coex,data,lb,ub,options,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                            otherwise
                                error('invalid Kp configuration');
                        end
                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                        [expcS,output] = parallel_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                        plot_tielines(output.aC,output.bC,C_bdy);
                        display = [resnorm;x];
                        disp(sprintf('best chisquare = %f, with parameters [%f %f]\n',display));
                        % pack up output structures
                        fit_output = struct('Chisquare',output.ChiSquare,'search_params',x','condn',output.condn);
                        spectra_output = struct('magnetic_field',spectra.magnetic_fields,'aS',output.aS,'bS',output.bS,...
                                    'obs_cS',S_coex,'exp_cS',output.exp_cS);
                        tieline_output = struct('aC',output.aC,'bC',output.bC,'fpafpb',output.fpafpb,'fafb',output.fafb,'Kp',output.Kp);
                        disp(sprintf('done!'));
                        return % fitting finished
                    case 'grid'
                        reply = input(sprintf('what is the range of slope angles to search? [0:180] (if unspecified, press enter)\n'));
                        disp(sprintf('parameterizing tieline slope...\n'));
                        if isempty(reply)
                            % default is 30 points for grid search
                            % parameter also thought of as angles from 0 to 180
                            s = linspace(0,1,30);
                        elseif isvector(reply) & any(reply >= 0 & reply <= 180)
                            s = reply.*(1/180);
                        else
                            error('invalid response');
                        end
                        % ask user whether to search Kp continuously or in a grid
                        reply = input(sprintf('slope is gridded, "continuous" or "grid" Kp search?\n'),'s');
                        while isempty(reply) || ~any(strcmp(reply,{'continuous','grid'}))
                            reply = input(sprintf('slope is gridded, "continuous" or "grid" Kp search?\n'),'s');
                        end
                        switch reply
                            case 'continuous'
                                options = optimset('MaxFunEvals',500,'TolFun',1e-8,'TolX',1e-8,'MaxIter',500,'LargeScale','on','LevenbergMarquardt','on','Display','off','Jacobian','off');
                                data = reshape(S_coex,ncB*ncS,1);
                                switch part_coeff.config
                                    case 'constant'
                                        reply = input(sprintf('what is the range of Kp to search? [1/10:10] (if unspecified, press enter)\n'));
                                        if isempty(reply)
                                            kp = [0.1 10];
                                        elseif isvector(reply) & any(reply >= 0.1 & reply <= 10)
                                            kp = reply;
                                        else
                                            error('invalid response');
                                        end
                                        disp(sprintf('starting search and fitting...'));
                                        % loop over slopes
                                        nslopes = length(s);
                                        for i = 1:nslopes
                                            x0 = [s(i) kp(1)+(kp(end)-kp(1)).*0.5]';% kp(1)+(kp(end)-kp(1)).*rand(1)
                                            lb = [s(i)-0.001 kp(1)]';
                                            ub = [s(i)+0.001 kp(end)]';
                                            [x,resnorm] = lsqcurvefit(@parallel_tieline_fit_fxn,x0,C_coex,data,lb,ub,options,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                            [expcS,output] = parallel_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                            % output =
                                            % struct('aC',alphaC,'bC',betaC,'aS',alphaS,'bS',betaS,
                                            % 'exp_cS',expcS,'fpafpb',fpafpb,'fafb',faf
                                            % b,'Kp',Kp,'ChiSquare',ChiSquare,'condn',condn);
                                            Chi2(i) = output.ChiSquare;
                                            % pack up variables for output
                                            fit_output(i) = struct('Chisquare',output.ChiSquare,'search_params',x','condn',output.condn);
                                            spectra_output(i) = struct('magnetic_field',spectra.magnetic_fields,'aS',output.aS,'bS',output.bS,...
                                                'obs_cS',S_coex,'exp_cS',output.exp_cS);
                                            tieline_output(i) = struct('aC',output.aC,'bC',output.bC,'fpafpb',output.fpafpb,'fafb',output.fafb,'Kp',output.Kp);
                                            pack;
                                            disp(sprintf('%f%% done',(i/nslopes)*100));
                                        end
                                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                                        % analyze results of search and fitting   
                                        i = find(min(Chi2(:)) == Chi2);
                                        if isvector(i) 
                                            for ii = 1:length(i)
                                                best_chisquare(ii,1) = Chi2(i(ii));
                                                best_indices(ii,1) = i(ii);
                                                best_params(ii,:) = fit_output(i(ii)).search_params;
                                                % plot boundary and alpha and beta phase compositions
                                                % and plot tielines.
                                                plot_tielines(tieline_output(i(ii)).aC,tieline_output(i(ii)).bC,C_bdy);
                                            end 
                                        else
                                            best_chisquare = Chi2(i);
                                            best_indices = i;
                                            best_params = fit_output(i).search_params;
                                            % plot boundary and alpha and beta phase compositions
                                            plot_tielines(tieline_output(i).aC,tieline_output(i).bC,C_bdy);
                                        end
                                        display = [best_chisquare';best_indices';best_params'];
                                        disp(sprintf('best chisquare = %f, with index %d and parameters [%f %f]\n',display));
                %                         disp(sprintf('best slope = %f => angle = %f\n',[best_params;slope2angle(best_params)]));
                                        disp(sprintf('done!'));
                                        return % fitting finished
                                    case 'function'
                                        % loop over slopes
                                        disp(sprintf('starting search and fitting...'));
                                        nslopes = length(s);
                                        for i = 1:nslopes
                                            x0 = [s(i) 0]';
                                            lb = [s(i)-0.001 -Inf]';
                                            ub = [s(i)+0.001 Inf]';
                                            [x,resnorm] = lsqcurvefit(@parallel_tieline_fit_fxn,x0,C_coex,data,lb,ub,options,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                            [expcS,output] = parallel_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                            % output =
                                            % struct('aC',alphaC,'bC',betaC,'aS',alphaS,'bS',betaS,
                                            % 'exp_cS',expcS,'fpafpb',fpafpb,'fafb',faf
                                            % b,'Kp',Kp,'ChiSquare',ChiSquare,'condn',condn);
                                            Chi2(i) = output.ChiSquare;
                                            % pack up variables for output
                                            fit_output(i) = struct('Chisquare',output.ChiSquare,'search_params',x','condn',output.condn);
                                            spectra_output(i) = struct('magnetic_field',spectra.magnetic_fields,'aS',output.aS,'bS',output.bS,...
                                                'obs_cS',S_coex,'exp_cS',output.exp_cS);
                                            tieline_output(i) = struct('aC',output.aC,'bC',output.bC,'fpafpb',output.fpafpb,'fafb',output.fafb,'Kp',output.Kp);
                                            pack;
                                            disp(sprintf('%f%% done',(i/nslopes)*100));
                                        end
                                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                                        % analyze results of search and fitting   
                                        i = find(min(Chi2(:)) == Chi2);
                                        if isvector(i) 
                                            for ii = 1:length(i)
                                                best_chisquare(ii,1) = Chi2(i(ii));
                                                best_indices(ii,1) = i(ii);
                                                best_params(ii,:) = fit_output(i(ii)).search_params;
                                                % plot boundary and alpha and beta phase compositions
                                                % and plot tielines.
                                                plot_tielines(tieline_output(i(ii)).aC,tieline_output(i(ii)).bC,C_bdy);
                                            end 
                                        else
                                            best_chisquare = Chi2(i);
                                            best_indices = i;
                                            best_params = fit_output(i).search_params;
                                            % plot boundary and alpha and beta phase compositions
                                            plot_tielines(tieline_output(i).aC,tieline_output(i).bC,C_bdy);
                                        end
                                        display = [best_chisquare';best_indices';best_params'];
                                        disp(sprintf('best chisquare = %f, with index %d and parameters [%f %f]\n',display));
                %                         disp(sprintf('best slope = %f => angle = %f\n',[best_params;slope2angle(best_params)]));
                                        disp(sprintf('done!'));
                                        return % fitting finished
                                    otherwise
                                        error('invalid Kp configuration');
                                end     
                            case 'grid'
                                reply = input(sprintf('what is the range of Kp to search? [1/10:10] (if unspecified, press enter)\n'));
                                if isempty(reply)
                                    % default is 50 points for grid search
                                    kp = linspace(0.1,10,21);
                                elseif isvector(reply)
                                    kp = reply;
                                else
                                    error('invalid response');
                                end
                                disp(sprintf('starting search and fitting...'));
                                % loop over slopes
                                nslopes = length(s);
                                part_coeff.config = 'constant';
                                nkp = length(kp);
                                for i = 1:nslopes
                                    for k = 1:nkp
                                        [expcS,output] = parallel_tieline_fit_fxn([s(i);kp(k)],C_coex,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                        % output =
                                        % struct('aC',alphaC,'bC',betaC,'aS',alphaS,'bS',betaS,
                                        % 'exp_cS',expcS,'fpafpb',fpafpb,'fafb',faf
                                        % b,'Kp',Kp,'ChiSquare',ChiSquare,'condn',condn);
                                        Chi2(i,k) = output.ChiSquare;
                                        % pack up variables for output
                                        fit_output(i,k) = struct('Chisquare',output.ChiSquare,'search_params',[s(i) kp(k)],'condn',output.condn);
                                        spectra_output(i,k) = struct('magnetic_field',spectra.magnetic_fields,'aS',output.aS,'bS',output.bS,...
                                            'obs_cS',S_coex,'exp_cS',output.exp_cS);
                                        tieline_output(i,k) = struct('aC',output.aC,'bC',output.bC,'fpafpb',output.fpafpb,'fafb',output.fafb,'Kp',output.Kp);
                                    end
                                    pack;
                                    disp(sprintf('%f%% done',(i/nslopes)*100));
                                end
                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                % analyze results of search and fitting   
                                [i,k] = find(min(Chi2(:)) == Chi2);
                                if isvector(i) 
                                    for ii = 1:length(i)
                                        best_chisquare(ii,1) = Chi2(i(ii),k(ii));
                                        best_indices(ii,:) = [i(ii) k(ii)];
                                        best_params(ii,:) = fit_output(i(ii),k(ii)).search_params;
                                        % plot boundary and alpha and beta phase compositions
                                        % and plot tielines.
                                        plot_tielines(tieline_output(i(ii),k(ii)).aC,tieline_output(i(ii),k(ii)).bC,C_bdy);
                                    end 
                                else
                                    best_chisquare = Chi2(i,k);
                                    best_indices = [i k];
                                    best_params = fit_output(i,k).search_params;
                                    % plot boundary and alpha and beta phase compositions
                                    plot_tielines(tieline_output(i,k).aC,tieline_output(i,k).bC,C_bdy);
                                end
                                display = [best_chisquare';best_indices';best_params'];
                                disp(sprintf('best chisquare = %f, with indices [%d %d] and parameters [%f %f]\n',display));
        %                         disp(sprintf('best slope = %f => angle = %f\n',[best_params;slope2angle(best_params)]));
                                disp(sprintf('done!'));
                                return % fitting finished
                            otherwise
                                error('invalid response');
                        end
                    otherwise
                        error('invalid boundary search method');
                end
            case 'tangent'
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
                % analyze boundary structure to organize search
                switch boundary.critpt1position
                    case 'located' 
                        switch boundary.critpt2position
                            case 'located'
                                if boundary.critpt1parameter == boundary.critpt2parameter
                                    error('both critical points cannot be the same point');
                                end
                                c1 = boundary.critpt1parameter;
                                c2 = boundary.critpt2parameter;
                                % ask user whether to search Kp
                                % continuously or in a grid
                                reply = input(sprintf('both critical points are located, "continuous" or "grid" Kp search?\n'),'s');
                                while isempty(reply) || ~any(strcmp(reply,{'continuous','grid'}))
                                    reply = input(sprintf('both critical points are located, "continuous" or "grid" Kp search?\n'),'s');
                                end
                                switch reply
                                    case 'continuous'
                                        % perform search and fit
                                        disp(sprintf('starting search and fitting...'));
                                        options = optimset('MaxFunEvals',500,'TolFun',1e-8,'TolX',1e-8,'MaxIter',500,'LargeScale','on','LevenbergMarquardt','on','Display','iter','Jacobian','off');
                                        data = reshape(S_coex,ncB*ncS,1);
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
                                                reply = input(sprintf('enter 3 parameter starting vector as [%4.3f<=cpt1<=%4.3f %4.3f<=cpt2<=%4.3f %4.2f<=kp<=%4.2f]:(enter for default)\n',[c1(1);c1(end);c2(1);c2(end);kp(1);kp(end)]));
                                                if isempty(reply)
                                                    x0 = [c1(1)+(c1(end)-c1(1))*rand(1) c2(1)+(c2(end)-c2(1))*rand(1) kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                                elseif isvector(reply)
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end
                                                lb = [c1-0.001 c2-0.001 kp(1)]';
                                                ub = [c1+0.001 c2+0.001 kp(end)]';
                                                [x,resnorm] = lsqcurvefit(@tangent_tieline_fit_fxn,x0,C_coex,data,lb,ub,options,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                            case 'function'
                                                reply = input(sprintf('enter 3 parameter starting vector as [%4.3f<=cpt1<=%4.3f %4.3f<=cpt2<=%4.3f c]:(enter for default)\n',[c1(1);c1(end);c2(1);c2(end)]));
                                                if isempty(reply)
                                                    x0 = [c1(1)+(c1(end)-c1(1))*rand(1) c2(1)+(c2(end)-c2(1))*rand(1) 0]';% for 2 critical points, Kp function has 1 parameter 
                                                elseif isvector(reply)
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end 
                                                lb = [c1-0.001 c2-0.001 -Inf]';
                                                ub = [c1+0.001 c2+0.001 Inf]';
                                                [x,resnorm] = lsqcurvefit(@tangent_tieline_fit_fxn,x0,C_coex,data,lb,ub,options,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                            otherwise
                                                error('invalid Kp configuration');
                                        end
                                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                                        [expcS,output] = tangent_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                        plot_tielines(output.aC,output.bC,C_bdy);
                                        display = [resnorm;x];
                                        disp(sprintf('best chisquare = %f, with parameters [%f %f %f]\n',display));
                                        % pack up output structures
                                        fit_output = struct('Chisquare',output.ChiSquare,'search_params',x','condn',output.condn);
                                        spectra_output = struct('magnetic_field',spectra.magnetic_fields,'aS',output.aS,'bS',output.bS,...
                                                    'obs_cS',S_coex,'exp_cS',output.exp_cS);
                                        tieline_output = struct('aC',output.aC,'bC',output.bC,'fpafpb',output.fpafpb,'fafb',output.fafb,'Kp',output.Kp);
                                        disp(sprintf('done!'));
                                        return % fitting finished
                                    case 'grid'
                                        reply = input(sprintf('what is the range of Kp to search? [1/10:10] (if unspecified, press enter)\n'));
                                        if isempty(reply)
                                            % default is 50 points for grid search
                                            kp = linspace(0.1,10,21);
                                        elseif isvector(reply)
                                            kp = reply;
                                        else
                                            error('invalid response');
                                        end
                                        disp(sprintf('starting search and fitting...'));
                                        % loop over Kps
                                        part_coeff.config = 'constant'; % grid with fxn not practical
                                        nkp = length(kp);
                                        for k = 1:nkp
                                            [expcS,output] = tangent_tieline_fit_fxn([c1;c2;kp(k)],C_coex,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                            % output =
                                            % struct('aC',alphaC,'bC',betaC,'aS',alphaS,'bS',betaS,
                                            % 'exp_cS',expcS,'fpafpb',fpafpb,'fafb',faf
                                            % b,'Kp',Kp,'ChiSquare',ChiSquare,'condn',condn);
                                            Chi2(k) = output.ChiSquare;
                                            % pack up variables for output
                                            fit_output(k) = struct('Chisquare',output.ChiSquare,'search_params',[c1 c2 kp(k)],'condn',output.condn);
                                            spectra_output(k) = struct('magnetic_field',spectra.magnetic_fields,'aS',output.aS,'bS',output.bS,...
                                                'obs_cS',S_coex,'exp_cS',output.exp_cS);
                                            tieline_output(k) = struct('aC',output.aC,'bC',output.bC,'fpafpb',output.fpafpb,'fafb',output.fafb,'Kp',output.Kp);
                                            pack;
                                            disp(sprintf('%f%% done',(k/nkp)*100));
                                        end
                                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                                        % analyze results of search and fitting   
                                        k = find(min(Chi2) == Chi2);
                                        if isvector(k) 
                                            for i = 1:length(k)
                                                best_chisquare(i,1) = Chi2(k(i));
                                                best_indices(i,1) = k(i);
                                                best_params(i,:) = fit_output(k(i)).search_params;
                                                % plot boundary and alpha and beta phase compositions
                                                % and plot tielines.
                                                plot_tielines(tieline_output(k(i)).aC,tieline_output(k(i)).bC,C_bdy);
                                            end 
                                        else
                                            best_chisquare = Chi2(k);
                                            best_indices = k;
                                            best_params = fit_output(k).search_params;
                                            % plot boundary and alpha and beta phase compositions
                                            plot_tielines(tieline_output(k).aC,tieline_output(k).bC,C_bdy);
                                        end
                                        display = [best_chisquare';best_indices';best_params'];
                                        disp(sprintf('best chisquare = %f, with index %d and parameters [%f %f %f]\n',display));
                                        disp(sprintf('done!'));
                                        return % fitting finished
                                    otherwise
                                        error('invalid response');
                                end         
                            case {'range','unknown'}
                                if length(boundary.critpt2parameter) == 2 % picked from figure as 2 points or unknown
                                    if boundary.critpt2parameter(1) > boundary.critpt2parameter(2)
                                        start = boundary.critpt2parameter(2);
                                        finish = boundary.critpt2parameter(1);
                                    elseif boundary.critpt2parameter(1) < boundary.critpt2parameter(2)
                                        start = boundary.critpt2parameter(1);
                                        finish = boundary.critpt2parameter(2);
                                    else
                                        error('critical point 2 range is same point');
                                    end
                                    boundary.critpt2parameter = linspace(start,finish,50*(finish-start));
                                end
                                c1 = boundary.critpt1parameter;
                                c2 = boundary.critpt2parameter;
                                % perform search and fit
                                switch boundary.search_method
                                    case 'continuous'
                                        disp(sprintf('starting search and fitting...'));
                                        options = optimset('MaxFunEvals',500,'TolFun',1e-8,'TolX',1e-8,'MaxIter',500,'LargeScale','on','LevenbergMarquardt','on','Display','iter','Jacobian','off');
                                        data = reshape(S_coex,ncB*ncS,1);
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
                                                reply = input(sprintf('enter 3 parameter starting vector as [%4.3f<=cpt1<=%4.3f %4.3f<=cpt2<=%4.3f %4.2f<=kp<=%4.2f]:(enter for default)\n',[c1(1);c1(end);c2(1);c2(end);kp(1);kp(end)]));
                                                if isempty(reply)
                                                    x0 = [c1(1)+(c1(end)-c1(1))*rand(1) c2(1)+(c2(end)-c2(1))*rand(1) kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                                elseif isvector(reply)
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end
                                                lb = [c1-0.001 c2(1) kp(1)]';
                                                ub = [c1+0.001 c2(end) kp(end)]';
                                                [x,resnorm] = lsqcurvefit(@tangent_tieline_fit_fxn,x0,C_coex,data,lb,ub,options,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                            case 'function'
                                                reply = input(sprintf('enter 3 parameter starting vector as [%4.3f<=cpt1<=%4.3f %4.3f<=cpt2<=%4.3f c]:(enter for default)\n',[c1(1);c1(end);c2(1);c2(end)]));
                                                if isempty(reply)
                                                    x0 = [c1(1)+(c1(end)-c1(1))*rand(1) c2(1)+(c2(end)-c2(1))*rand(1) 0]';% for 2 critical points, Kp function has 1 parameter 
                                                elseif isvector(reply)
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end
                                                lb = [c1-0.001 c2(1) -Inf]';
                                                ub = [c1+0.001 c2(end) Inf]';
                                                [x,resnorm] = lsqcurvefit(@tangent_tieline_fit_fxn,x0,C_coex,data,lb,ub,options,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                            otherwise
                                                error('invalid Kp configuration');
                                        end
                                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                                        [expcS,output] = tangent_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                        plot_tielines(output.aC,output.bC,C_bdy);
                                        display = [resnorm;x];
                                        disp(sprintf('best chisquare = %f, with parameters [%f %f %f]\n',display));
                                        % pack up output structures
                                        fit_output = struct('Chisquare',output.ChiSquare,'search_params',x','condn',output.condn);
                                        spectra_output = struct('magnetic_field',spectra.magnetic_fields,'aS',output.aS,'bS',output.bS,...
                                                    'obs_cS',S_coex,'exp_cS',output.exp_cS);
                                        tieline_output = struct('aC',output.aC,'bC',output.bC,'fpafpb',output.fpafpb,'fafb',output.fafb,'Kp',output.Kp);
                                        disp(sprintf('done!'));
                                        return % fitting finished
                                    case 'grid'
                                        switch part_coeff.config
                                            case 'constant'
                                                % ask user whether to search Kp
                                                % continuously or in a grid
                                                reply = input(sprintf('since one critical point is gridded, "continuous" or "grid" Kp search?\n'),'s');
                                                while isempty(reply) || ~any(strcmp(reply,{'continuous','grid'}))
                                                    reply = input(sprintf('since one critical point is gridded, "continuous" or "grid" Kp search?\n'),'s');
                                                end
                                                switch reply
                                                    case 'continuous'
                                                        reply = input(sprintf('what is the range of Kp to search? [1/10:10] (if unspecified, press enter)\n'));
                                                        if isempty(reply)
                                                            kp = [0.1 10];
                                                        elseif isvector(reply) & all(reply >= 0.1 & reply <= 10)
                                                            kp = reply;
                                                        else
                                                            error('invalid response');
                                                        end
                                                        disp(sprintf('starting search and fitting...'));
                                                        options = optimset('MaxFunEvals',500,'TolFun',1e-8,'TolX',1e-8,'MaxIter',500,'LargeScale','on','LevenbergMarquardt','on','Display','off','Jacobian','off');
                                                        data = reshape(S_coex,ncB*ncS,1);
                                                        npositions = length(c2);
                                                        for i = 1:npositions
                                                            x0 = [c1 c2(i) kp(1)+(kp(end)-kp(1)).*0.5]';% kp(1)+(kp(end)-kp(1)).*rand(1)
                                                            lb = [c1-0.001 c2(i)-0.001 kp(1)]';
                                                            ub = [c1+0.001 c2(i)+0.001 kp(end)]';
                                                            [x,resnorm] = lsqcurvefit(@tangent_tieline_fit_fxn,x0,C_coex,data,lb,ub,options,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                                            [expcS,output] = tangent_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                                            % output =
                                                            % struct('aC',alphaC,'bC',betaC,'aS',alphaS,'bS',betaS,
                                                            % 'exp_cS',expcS,'fpafpb',fpafpb,'fafb',faf
                                                            % b,'Kp',Kp,'ChiSquare',ChiSquare,'condn',condn);
                                                            Chi2(i) = output.ChiSquare;
                                                            % pack up variables for output
                                                            fit_output(i) = struct('Chisquare',output.ChiSquare,'search_params',x','condn',output.condn);
                                                            spectra_output(i) = struct('magnetic_field',spectra.magnetic_fields,'aS',output.aS,'bS',output.bS,...
                                                                'obs_cS',S_coex,'exp_cS',output.exp_cS);
                                                            tieline_output(i) = struct('aC',output.aC,'bC',output.bC,'fpafpb',output.fpafpb,'fafb',output.fafb,'Kp',output.Kp);
                                                            pack;
                                                            disp(sprintf('%f%% done',(i/npositions)*100));
                                                        end
                                                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                                                        % analyze results of search and fitting   
                                                        i = find(min(Chi2) == Chi2);
                                                        if isvector(i) 
                                                            for ii = 1:length(i)
                                                                best_chisquare(ii,1) = Chi2(i(ii));
                                                                best_indices(ii) = i(ii);
                                                                best_params(ii,:) = fit_output(i(ii)).search_params;
                                                                % plot boundary and alpha and beta phase compositions
                                                                % and plot tielines.
                                                                plot_tielines(tieline_output(i(ii)).aC,tieline_output(i(ii)).bC,C_bdy);
                                                            end 
                                                        else
                                                            best_chisquare = Chi2(i);
                                                            best_indices = i;
                                                            best_params = fit_output(i).search_params;
                                                            % plot boundary and alpha and beta phase compositions
                                                            plot_tielines(tieline_output(i).aC,tieline_output(i).bC,C_bdy);
                                                        end
                                                        display = [best_chisquare';best_indices';best_params'];
                                                        disp(sprintf('best chisquare = %f, with index %d and parameters [%f %f %f]\n',display));
                                %                         disp(sprintf('best slope = %f => angle = %f\n',[best_params;slope2angle(best_params)]));
                                                        disp(sprintf('done!'));
                                                        return % fitting finished
                                                    case 'grid'
                                                        reply = input(sprintf('what is the range of Kp to search? [1/10:10] (if unspecified, press enter)\n'));
                                                        if isempty(reply)
                                                            % default is 21 points for grid search
                                                            kp = linspace(0.1,10,21);
                                                        elseif isvector(reply)
                                                            kp = reply;
                                                        else
                                                            error('invalid response');
                                                        end
                                                        disp(sprintf('starting search and fitting...'));
                                                        npositions = length(c2);
                                                        nkp = length(kp);
                                                        ngridpts = nkp*npositions;
                                                        for i = 1:npositions
                                                            for k = 1:nkp
                                                                [expcS,output] = tangent_tieline_fit_fxn([c1;c2(i);kp(k)],C_coex,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                                                % output =
                                                                % struct('aC',alphaC,'bC',betaC,'aS',alphaS,'bS',betaS,
                                                                % 'exp_cS',expcS,'fpafpb',fpafpb,'fafb',faf
                                                                % b,'Kp',Kp,'ChiSquare',ChiSquare,'condn',condn);
                                                                Chi2(i,k) = output.ChiSquare;
                                                                % pack up variables for output
                                                                fit_output(i,k) = struct('Chisquare',output.ChiSquare,'search_params',[c1 c2(i) kp(k)],'condn',output.condn);
                                                                spectra_output(i,k) = struct('magnetic_field',spectra.magnetic_fields,'aS',output.aS,'bS',output.bS,...
                                                                    'obs_cS',S_coex,'exp_cS',output.exp_cS);
                                                                tieline_output(i,k) = struct('aC',output.aC,'bC',output.bC,'fpafpb',output.fpafpb,'fafb',output.fafb,'Kp',output.Kp);
                                                            end
                                                            pack;
                                                            disp(sprintf('%f%% done',(i/npositions)*100));
                                                        end
                                                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                                                        % analyze results of search and fitting   
                                                        [i,k] = find(min(Chi2(:)) == Chi2);
                                                        if isvector(i) 
                                                            for ii = 1:length(i)
                                                                best_chisquare(ii,1) = Chi2(i(ii),k(ii));
                                                                best_indices(ii,:) = [i(ii) k(ii)];
                                                                best_params(ii,:) = fit_output(i(ii),k(ii)).search_params;
                                                                % plot boundary and alpha and beta phase compositions
                                                                % and plot tielines.
                                                                plot_tielines(tieline_output(i(ii),k(ii)).aC,tieline_output(i(ii),k(ii)).bC,C_bdy);
                                                            end 
                                                        else
                                                            best_chisquare = Chi2(i,k);
                                                            best_indices = [i k];
                                                            best_params = fit_output(i,k).search_params;
                                                            % plot boundary and alpha and beta phase compositions
                                                            plot_tielines(tieline_output(i,k).aC,tieline_output(i,k).bC,C_bdy);
                                                        end
                                                        display = [best_chisquare';best_indices';best_params'];
                                                        disp(sprintf('best chisquare = %f, with indices [%d %d] and parameters [%f %f %f]\n',display));
                                %                         disp(sprintf('best slope = %f => angle = %f\n',[best_params;slope2angle(best_params)]));
                                                        disp(sprintf('done!'));
                                                        return % fitting finished
                                                    otherwise
                                                        error('invalid response');
                                                end
                                            case 'function'
                                                disp(sprintf('starting search and fitting...'));
                                                options = optimset('MaxFunEvals',500,'TolFun',1e-8,'TolX',1e-8,'MaxIter',500,'LargeScale','on','LevenbergMarquardt','on','Display','off','Jacobian','off');
                                                data = reshape(S_coex,ncB*ncS,1);
                                                npositions = length(c2);
                                                for i = 1:npositions
                                                    x0 = [c1 c2(i) 0]';
                                                    lb = [c1-0.001 c2(i)-0.001 -Inf]';
                                                    ub = [c1+0.001 c2(i)+0.001 Inf]';
                                                    [x,resnorm] = lsqcurvefit(@tangent_tieline_fit_fxn,x0,C_coex,data,lb,ub,options,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                                    [expcS,output] = tangent_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                                    % output =
                                                    % struct('aC',alphaC,'bC',betaC,'aS',alphaS,'bS',betaS,
                                                    % 'exp_cS',expcS,'fpafpb',fpafpb,'fafb',faf
                                                    % b,'Kp',Kp,'ChiSquare',ChiSquare,'condn',condn);
                                                    Chi2(i) = output.ChiSquare;
                                                    % pack up variables for output
                                                    fit_output(i) = struct('Chisquare',output.ChiSquare,'search_params',x','condn',output.condn);
                                                    spectra_output(i) = struct('magnetic_field',spectra.magnetic_fields,'aS',output.aS,'bS',output.bS,...
                                                        'obs_cS',S_coex,'exp_cS',output.exp_cS);
                                                    tieline_output(i) = struct('aC',output.aC,'bC',output.bC,'fpafpb',output.fpafpb,'fafb',output.fafb,'Kp',output.Kp);
                                                    pack;
                                                    disp(sprintf('%f%% done',(i/npositions)*100));
                                                end
                                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                                % analyze results of search and fitting   
                                                i = find(min(Chi2) == Chi2);
                                                if isvector(i) 
                                                    for ii = 1:length(i)
                                                        best_chisquare(ii,1) = Chi2(i(ii));
                                                        best_indices(ii) = i(ii);
                                                        best_params(ii,:) = fit_output(i(ii)).search_params;
                                                        % plot boundary and alpha and beta phase compositions
                                                        % and plot tielines.
                                                        plot_tielines(tieline_output(i(ii),k(ii)).aC,tieline_output(i(ii),k(ii)).bC,C_bdy);
                                                    end 
                                                else
                                                    best_chisquare = Chi2(i);
                                                    best_indices = i;
                                                    best_params = fit_output(i).search_params;
                                                    % plot boundary and alpha and beta phase compositions
                                                    plot_tielines(tieline_output(i,k).aC,tieline_output(i,k).bC,C_bdy);
                                                end
                                                display = [best_chisquare';best_indices';best_params'];
                                                disp(sprintf('best chisquare = %f, with indices %d and parameters [%f %f %f]\n',display));
                        %                         disp(sprintf('best slope = %f => angle = %f\n',[best_params;slope2angle(best_params)]));
                                                disp(sprintf('done!'));
                                                return % fitting finished  
                                            otherwise
                                                error('invalid Kp configuration');
                                        end
                                    otherwise
                                        error('invalid boundary search method');
                                end
                            otherwise
                                error('invalid critical point 2 position on boundary');
                        end
                    case {'range','unknown'} % for critical point 1
                        if length(boundary.critpt1parameter) == 2 % picked from figure as 2 points or unknown
                            if boundary.critpt1parameter > boundary.critpt1parameter(2)
                                start = boundary.critpt1parameter(2);
                                finish = boundary.critpt1parameter(1);
                            elseif boundary.critpt1parameter(1) < boundary.critpt1parameter(2)
                                start = boundary.critpt1parameter(1);
                                finish = boundary.critpt1parameter(2);
                            else
                                error('critical point 1 range is same point');
                            end
                            boundary.critpt1parameter = linspace(start,finish,50*(finish-start));
                        end
                        switch boundary.critpt2position
                            case 'located'
%                                 if any(boundary.critpt1parameter == boundary.critpt2parameter)
%                                     error('both critical points cannot be the same point');
%                                 end
                                % perform search and fit
                                c1 = boundary.critpt1parameter;
                                c2 = boundary.critpt2parameter;
                                switch boundary.search_method
                                    case 'continuous'
                                        disp(sprintf('starting search and fitting...'));
                                        options = optimset('MaxFunEvals',500,'TolFun',1e-8,'TolX',1e-8,'MaxIter',500,'LargeScale','on','LevenbergMarquardt','on','Display','iter','Jacobian','off');
                                        data = reshape(S_coex,ncB*ncS,1);
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
                                                reply = input(sprintf('enter 3 parameter starting vector as [%4.3f<=cpt1<=%4.3f %4.3f<=cpt2<=%4.3f %4.2f<=kp<=%4.2f]:(enter for default)\n',[c1(1);c1(end);c2(1);c2(end);kp(1);kp(end)]));
                                                if isempty(reply)
                                                    x0 = [c1(1)+(c1(end)-c1(1))*rand(1) c2(1)+(c2(end)-c2(1))*rand(1) kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                                elseif isvector(reply)
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end
                                                lb = [c1(1) c2-0.001 kp(1)]';
                                                ub = [c1(end) c2+0.001 kp(end)]';
                                                [x,resnorm] = lsqcurvefit(@tangent_tieline_fit_fxn,x0,C_coex,data,lb,ub,options,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                            case 'function'
                                                reply = input(sprintf('enter 3 parameter starting vector as [%4.3f<=cpt1<=%4.3f %4.3f<=cpt2<=%4.3f c]:(enter for default)\n',[c1(1);c1(end);c2(1);c2(end)]));
                                                if isempty(reply)
                                                    x0 = [c1(1)+(c1(end)-c1(1))*rand(1) c2(1)+(c2(end)-c2(1))*rand(1) 0]';% for 2 critical points, Kp function has 1 parameter 
                                                elseif isvector(reply)
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end
                                                lb = [c1(1) c2-0.001 -Inf]';
                                                ub = [c1(end) c2+0.001 Inf]';
                                                [x,resnorm] = lsqcurvefit(@tangent_tieline_fit_fxn,x0,C_coex,data,lb,ub,options,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                            otherwise
                                                error('invalid Kp configuration');
                                        end
                                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                                        [expcS,output] = tangent_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                        plot_tielines(output.aC,output.bC,C_bdy);
                                        display = [resnorm;x];
                                        disp(sprintf('best chisquare = %f, with parameters [%f %f %f]\n',display));
                                        % pack up output structures
                                        fit_output = struct('Chisquare',output.ChiSquare,'search_params',x','condn',output.condn);
                                        spectra_output = struct('magnetic_field',spectra.magnetic_fields,'aS',output.aS,'bS',output.bS,...
                                                    'obs_cS',S_coex,'exp_cS',output.exp_cS);
                                        tieline_output = struct('aC',output.aC,'bC',output.bC,'fpafpb',output.fpafpb,'fafb',output.fafb,'Kp',output.Kp);
                                        disp(sprintf('done!'));
                                        return % fitting finished
                                    case 'grid'
                                        switch part_coeff.config
                                            case 'constant'
                                                % ask user whether to search Kp
                                                % continuously or in a grid
                                                reply = input(sprintf('one critical point is gridded, "continuous" or "grid" Kp search?\n'),'s');
                                                while isempty(reply) || ~any(strcmp(reply,{'continuous','grid'}))
                                                    reply = input(sprintf('one critical point is gridded, "continuous" or "grid" Kp search?\n'),'s');
                                                end
                                                switch reply
                                                    case 'continuous'
                                                        reply = input(sprintf('what is the range of Kp to search? [1/10:10] (if unspecified, press enter)\n'));
                                                        if isempty(reply)
                                                            kp = [0.1 10];
                                                        elseif isvector(reply) & all(reply >= 0.1 & reply <= 10)
                                                            kp = reply;
                                                        else
                                                            error('invalid response');
                                                        end
                                                        disp(sprintf('starting search and fitting...'));
                                                        options = optimset('MaxFunEvals',500,'TolFun',1e-8,'TolX',1e-8,'MaxIter',500,'LargeScale','on','LevenbergMarquardt','on','Display','off','Jacobian','off');
                                                        data = reshape(S_coex,ncB*ncS,1);
                                                        npositions = length(c1);
                                                        for i = 1:npositions
                                                            x0 = [c1(i) c2 kp(1)+(kp(end)-kp(1)).*0.5]';% kp(1)+(kp(end)-kp(1)).*rand(1)
                                                            lb = [c1(i)-0.001 c2-0.001 kp(1)]';
                                                            ub = [c1(i)+0.001 c2+0.001 kp(end)]';
                                                            [x,resnorm] = lsqcurvefit(@tangent_tieline_fit_fxn,x0,C_coex,data,lb,ub,options,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                                            [expcS,output] = tangent_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                                            % output =
                                                            % struct('aC',alphaC,'bC',betaC,'aS',alphaS,'bS',betaS,
                                                            % 'exp_cS',expcS,'fpafpb',fpafpb,'fafb',faf
                                                            % b,'Kp',Kp,'ChiSquare',ChiSquare,'condn',condn);
                                                            Chi2(i) = output.ChiSquare;
                                                            % pack up variables for output
                                                            fit_output(i) = struct('Chisquare',output.ChiSquare,'search_params',x','condn',output.condn);
                                                            spectra_output(i) = struct('magnetic_field',spectra.magnetic_fields,'aS',output.aS,'bS',output.bS,...
                                                                'obs_cS',S_coex,'exp_cS',output.exp_cS);
                                                            tieline_output(i) = struct('aC',output.aC,'bC',output.bC,'fpafpb',output.fpafpb,'fafb',output.fafb,'Kp',output.Kp);
                                                            pack;
                                                            disp(sprintf('%f%% done',(i/npositions)*100));
                                                        end
                                                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                                                        % analyze results of search and fitting   
                                                        i = find(min(Chi2) == Chi2);
                                                        if isvector(i) 
                                                            for ii = 1:length(i)
                                                                best_chisquare(ii,1) = Chi2(i(ii));
                                                                best_indices(ii) = i(ii);
                                                                best_params(ii,:) = fit_output(i(ii)).search_params;
                                                                % plot boundary and alpha and beta phase compositions
                                                                % and plot tielines.
                                                                plot_tielines(tieline_output(i(ii)).aC,tieline_output(i(ii)).bC,C_bdy);
                                                            end 
                                                        else
                                                            best_chisquare = Chi2(i);
                                                            best_indices = i;
                                                            best_params = fit_output(i).search_params;
                                                            % plot boundary and alpha and beta phase compositions
                                                            plot_tielines(tieline_output(i).aC,tieline_output(i).bC,C_bdy);
                                                        end
                                                        display = [best_chisquare';best_indices';best_params'];
                                                        disp(sprintf('best chisquare = %f, with index %d and parameters [%f %f %f]\n',display));
                                %                         disp(sprintf('best slope = %f => angle = %f\n',[best_params;slope2angle(best_params)]));
                                                        disp(sprintf('done!'));
                                                        return % fitting finished
                                                    case 'grid'
                                                        reply = input(sprintf('what is the range of Kp to search? [1/10:10] (if unspecified, press enter)\n'));
                                                        if isempty(reply)
                                                            % default is 21 points for grid search
                                                            kp = linspace(0.1,10,21);
                                                        elseif isvector(reply)
                                                            kp = reply;
                                                        else
                                                            error('invalid response');
                                                        end
                                                        disp(sprintf('starting search and fitting...'));
                                                        npositions = length(c1);
                                                        nkp = length(kp);
                                                        ngridpts = nkp*npositions;
                                                        for i = 1:npositions
                                                            for k = 1:nkp
                                                                [expcS,output] = tangent_tieline_fit_fxn([c1(i);c2;kp(k)],C_coex,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                                                % output =
                                                                % struct('aC',alphaC,'bC',betaC,'aS',alphaS,'bS',betaS,
                                                                % 'exp_cS',expcS,'fpafpb',fpafpb,'fafb',faf
                                                                % b,'Kp',Kp,'ChiSquare',ChiSquare,'condn',condn);
                                                                Chi2(i,k) = output.ChiSquare;
                                                                % pack up variables for output
                                                                fit_output(i,k) = struct('Chisquare',output.ChiSquare,'search_params',[c1(i) c2 kp(k)],'condn',output.condn);
                                                                spectra_output(i,k) = struct('magnetic_field',spectra.magnetic_fields,'aS',output.aS,'bS',output.bS,...
                                                                    'obs_cS',S_coex,'exp_cS',output.exp_cS);
                                                                tieline_output(i,k) = struct('aC',output.aC,'bC',output.bC,'fpafpb',output.fpafpb,'fafb',output.fafb,'Kp',output.Kp);
                                                            end
                                                            pack;
                                                            disp(sprintf('%f%% done',(i/npositions)*100));
                                                        end
                                                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                                                        % analyze results of search and fitting   
                                                        [i,k] = find(min(Chi2(:)) == Chi2);
                                                        if isvector(i) 
                                                            for ii = 1:length(i)
                                                                best_chisquare(ii,1) = Chi2(i(ii),k(ii));
                                                                best_indices(ii,:) = [i(ii) k(ii)];
                                                                best_params(ii,:) = fit_output(i(ii),k(ii)).search_params;
                                                                % plot boundary and alpha and beta phase compositions
                                                                % and plot tielines.
                                                                plot_tielines(tieline_output(i(ii),k(ii)).aC,tieline_output(i(ii),k(ii)).bC,C_bdy);
                                                            end 
                                                        else
                                                            best_chisquare = Chi2(i,k);
                                                            best_indices = [i k];
                                                            best_params = fit_output(i,k).search_params;
                                                            % plot boundary and alpha and beta phase compositions
                                                            plot_tielines(tieline_output(i,k).aC,tieline_output(i,k).bC,C_bdy);
                                                        end
                                                        display = [best_chisquare';best_indices';best_params'];
                                                        disp(sprintf('best chisquare = %f, with indices [%d %d] and parameters [%f %f %f]\n',display));
                                %                         disp(sprintf('best slope = %f => angle = %f\n',[best_params;slope2angle(best_params)]));
                                                        disp(sprintf('done!'));
                                                        return % fitting finished
                                                    otherwise
                                                        error('invalid response');
                                                end
                                            case 'function'
                                                disp(sprintf('starting search and fitting...'));
                                                options = optimset('MaxFunEvals',500,'TolFun',1e-8,'TolX',1e-8,'MaxIter',500,'LargeScale','on','LevenbergMarquardt','on','Display','off','Jacobian','off');
                                                data = reshape(S_coex,ncB*ncS,1);
                                                npositions = length(c2);
                                                for i = 1:npositions
                                                    x0 = [c1(i) c2 0]';
                                                    lb = [c1(i)-0.001 c2-0.001 -Inf]';
                                                    ub = [c1(i)+0.001 c2+0.001 Inf]';
                                                    [x,resnorm] = lsqcurvefit(@tangent_tieline_fit_fxn,x0,C_coex,data,lb,ub,options,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                                    [expcS,output] = tangent_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                                    % output =
                                                    % struct('aC',alphaC,'bC',betaC,'aS',alphaS,'bS',betaS,
                                                    % 'exp_cS',expcS,'fpafpb',fpafpb,'fafb',faf
                                                    % b,'Kp',Kp,'ChiSquare',ChiSquare,'condn',condn);
                                                    Chi2(i) = output.ChiSquare;
                                                    % pack up variables for output
                                                    fit_output(i) = struct('Chisquare',output.ChiSquare,'search_params',x','condn',output.condn);
                                                    spectra_output(i) = struct('magnetic_field',spectra.magnetic_fields,'aS',output.aS,'bS',output.bS,...
                                                        'obs_cS',S_coex,'exp_cS',output.exp_cS);
                                                    tieline_output(i) = struct('aC',output.aC,'bC',output.bC,'fpafpb',output.fpafpb,'fafb',output.fafb,'Kp',output.Kp);
                                                    pack;
                                                    disp(sprintf('%f%% done',(i/npositions)*100));
                                                end
                                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                                % analyze results of search and fitting   
                                                i = find(min(Chi2) == Chi2);
                                                if isvector(i) 
                                                    for ii = 1:length(i)
                                                        best_chisquare(ii,1) = Chi2(i(ii));
                                                        best_indices(ii) = i(ii);
                                                        best_params(ii,:) = fit_output(i(ii)).search_params;
                                                        % plot boundary and alpha and beta phase compositions
                                                        % and plot tielines.
                                                        plot_tielines(tieline_output(i(ii)).aC,tieline_output(i(ii)).bC,C_bdy);
                                                    end 
                                                else
                                                    best_chisquare = Chi2(i);
                                                    best_indices = i;
                                                    best_params = fit_output(i).search_params;
                                                    % plot boundary and alpha and beta phase compositions
                                                    plot_tielines(tieline_output(i).aC,tieline_output(i).bC,C_bdy);
                                                end
                                                display = [best_chisquare';best_indices';best_params'];
                                                disp(sprintf('best chisquare = %f, with indices %d and parameters [%f %f %f]\n',display));
                        %                         disp(sprintf('best slope = %f => angle = %f\n',[best_params;slope2angle(best_params)]));
                                                disp(sprintf('done!'));
                                                return % fitting finished  
                                            otherwise
                                                error('invalid Kp configuration');
                                        end
                                    otherwise
                                        error('invalid boundary search method');
                                end
                            case {'range','unknown'} % for critical point 2, critical point 1 also range or unknown
                                if length(boundary.critpt2parameter) == 2 % picked from figure as 2 points
                                    if boundary.critpt2parameter(1) > boundary.critpt2parameter(2)
                                        start = boundary.critpt2parameter(2);
                                        finish = boundary.critpt2parameter(1);
                                    elseif boundary.critpt2parameter(1) < boundary.critpt2parameter(2)
                                        start = boundary.critpt2parameter(1);
                                        finish = boundary.critpt2parameter(2);
                                    else
                                        error('critical point 2 range is same point');
                                    end
                                    boundary.critpt2parameter = linspace(start,finish,50*(finish-start));
                                end
                                % perform search and fit
                                c1 = boundary.critpt1parameter;
                                c2 = boundary.critpt2parameter;
                                switch boundary.search_method
                                    case 'continuous'
                                        disp(sprintf('starting search and fitting...'));
                                        options = optimset('MaxFunEvals',500,'TolFun',1e-8,'TolX',1e-8,'MaxIter',500,'LargeScale','on','LevenbergMarquardt','on','Display','iter','Jacobian','off');
                                        data = reshape(S_coex,ncB*ncS,1);
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
                                                reply = input(sprintf('enter 3 parameter starting vector as [%4.3f<=cpt1<=%4.3f %4.3f<=cpt2<=%4.3f %4.2f<=kp<=%4.2f]:(enter for default)\n',[c1(1);c1(end);c2(1);c2(end);kp(1);kp(end)]));
                                                if isempty(reply)
                                                    x0 = [c1(1)+(c1(end)-c1(1))*rand(1) c2(1)+(c2(end)-c2(1))*rand(1) kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                                elseif isvector(reply)
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end
                                                lb = [c1(1) c2(1) kp(1)]';
                                                ub = [c1(end) c2(end) kp(end)]';
                                                [x,resnorm] = lsqcurvefit(@tangent_tieline_fit_fxn,x0,C_coex,data,lb,ub,options,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                            case 'function'
                                                reply = input(sprintf('enter 3 parameter starting vector as [%4.3f<=cpt1<=%4.3f %4.3f<=cpt2<=%4.3f c]:(enter for default)\n',[c1(1);c1(end);c2(1);c2(end)]));
                                                if isempty(reply)
                                                    x0 = [c1(1)+(c1(end)-c1(1))*rand(1) c2(1)+(c2(end)-c2(1))*rand(1) 0]';% for 2 critical points, Kp function has 1 parameter 
                                                elseif isvector(reply)
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end
                                                lb = [c1(1) c2(1) -Inf]';
                                                ub = [c1(end) c2(end) Inf]';
                                                [x,resnorm] = lsqcurvefit(@tangent_tieline_fit_fxn,x0,C_coex,data,lb,ub,options,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                            otherwise
                                                error('invalid Kp configuration');
                                        end
                                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                                        [expcS,output] = tangent_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                        plot_tielines(output.aC,output.bC,C_bdy);
                                        display = [resnorm;x];
                                        disp(sprintf('best chisquare = %f, with parameters [%f %f %f]\n',display));
                                        % pack up output structures
                                        fit_output = struct('Chisquare',output.ChiSquare,'search_params',x','condn',output.condn);
                                        spectra_output = struct('magnetic_field',spectra.magnetic_fields,'aS',output.aS,'bS',output.bS,...
                                                    'obs_cS',S_coex,'exp_cS',output.exp_cS);
                                        tieline_output = struct('aC',output.aC,'bC',output.bC,'fpafpb',output.fpafpb,'fafb',output.fafb,'Kp',output.Kp);
                                        disp(sprintf('done!'));
                                        return % fitting finished
                                    case 'grid'
                                        switch part_coeff.config
                                            case 'constant'
                                                % ask user whether to search Kp
                                                % continuously or in a grid
                                                reply = input(sprintf('both critical points are gridded, "continuous" or "grid" Kp search?\n'),'s');
                                                while isempty(reply) || ~any(strcmp(reply,{'continuous','grid'}))
                                                    reply = input(sprintf('both critical points are gridded, "continuous" or "grid" Kp search?\n'),'s');
                                                end
                                                switch reply
                                                    case 'continuous'
                                                        reply = input(sprintf('what is the range of Kp to search? [1/10:10] (if unspecified, press enter)\n'));
                                                        if isempty(reply)
                                                            kp = [0.1 10];
                                                        elseif isvector(reply) & all(reply >= 0.1 & reply <= 10)
                                                            kp = reply;
                                                        else
                                                            error('invalid response');
                                                        end
                                                        disp(sprintf('starting search and fitting...'));
                                                        options = optimset('MaxFunEvals',500,'TolFun',1e-8,'TolX',1e-8,'MaxIter',500,'LargeScale','on','LevenbergMarquardt','on','Display','off','Jacobian','off');
                                                        data = reshape(S_coex,ncB*ncS,1);
                                                        npositions1 = length(c1);
                                                        npositions2 = length(c2);
                                                        for i = 1:npositions1
                                                            for j = 1:npositions2
                                                                x0 = [c1(i) c2(j) kp(1)+(kp(end)-kp(1)).*0.5]';% kp(1)+(kp(end)-kp(1)).*rand(1)
                                                                lb = [c1(i)-0.001 c2(j)-0.001 kp(1)]';
                                                                ub = [c1(i)+0.001 c2(j)+0.001 kp(end)]';
                                                                [x,resnorm] = lsqcurvefit(@tangent_tieline_fit_fxn,x0,C_coex,data,lb,ub,options,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                                                [expcS,output] = tangent_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                                                % output =
                                                                % struct('aC',alphaC,'bC',betaC,'aS',alphaS,'bS',betaS,
                                                                % 'exp_cS',expcS,'fpafpb',fpafpb,'fafb',faf
                                                                % b,'Kp',Kp,'ChiSquare',ChiSquare,'condn',condn);
                                                                Chi2(i,j) = output.ChiSquare;
                                                                % pack up variables for output
                                                                fit_output(i,j) = struct('Chisquare',output.ChiSquare,'search_params',x','condn',output.condn);
                                                                spectra_output(i,j) = struct('magnetic_field',spectra.magnetic_fields,'aS',output.aS,'bS',output.bS,...
                                                                    'obs_cS',S_coex,'exp_cS',output.exp_cS);
                                                                tieline_output(i,j) = struct('aC',output.aC,'bC',output.bC,'fpafpb',output.fpafpb,'fafb',output.fafb,'Kp',output.Kp);
                                                            end
                                                            pack;
                                                            disp(sprintf('%f%% done',(i/npositions1)*100));
                                                        end
                                                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                                                        % analyze results of search and fitting   
                                                        [i,j] = find(min(Chi2(:)) == Chi2);
                                                        if isvector(i) 
                                                            for ii = 1:length(i)
                                                                best_chisquare(ii,1) = Chi2(i(ii),j(ii));
                                                                best_indices(ii,:) = [i(ii) j(ii)];
                                                                best_params(ii,:) = fit_output(i(ii),j(ii)).search_params;
                                                                % plot boundary and alpha and beta phase compositions
                                                                % and plot tielines.
                                                                plot_tielines(tieline_output(i(ii),j(ii)).aC,tieline_output(i(ii),j(ii)).bC,C_bdy);
                                                            end 
                                                        else
                                                            best_chisquare = Chi2(i,j);
                                                            best_indices = [i j];
                                                            best_params = fit_output(i,j).search_params;
                                                            % plot boundary and alpha and beta phase compositions
                                                            plot_tielines(tieline_output(i,j).aC,tieline_output(i,j).bC,C_bdy);
                                                        end
                                                        display = [best_chisquare';best_indices';best_params'];
                                                        disp(sprintf('best chisquare = %f, with indices [%d %d] and parameters [%f %f %f]\n',display));
                                %                         disp(sprintf('best slope = %f => angle = %f\n',[best_params;slope2angle(best_params)]));
                                                        disp(sprintf('done!'));
                                                        return % fitting finished
                                                    case 'grid'
                                                        reply = input(sprintf('what is the range of Kp to search? [1/10:10] (if unspecified, press enter)\n'));
                                                        if isempty(reply)
                                                            % default is 21 points for grid search
                                                            kp = linspace(0.1,10,21);
                                                        elseif isvector(reply)
                                                            kp = reply;
                                                        else
                                                            error('invalid response');
                                                        end
                                                        disp(sprintf('starting search and fitting...'));
                                                        npositions1 = length(c1);
                                                        npositions2 = length(c2);
                                                        nkp = length(kp);
                                                        for i = 1:npositions1
                                                            for j = 1:npositions2
                                                                for k = 1:nkp
                                                                    [expcS,output] = tangent_tieline_fit_fxn([c1(i);c2(j);kp(k)],C_coex,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                                                    % output =
                                                                    % struct('aC',alphaC,'bC',betaC,'aS',alphaS,'bS',betaS,
                                                                    % 'exp_cS',expcS,'fpafpb',fpafpb,'fafb',faf
                                                                    % b,'Kp',Kp,'ChiSquare',ChiSquare,'condn',condn);
                                                                    Chi2(i,j,k) = output.ChiSquare;
                                                                    % pack up variables for output
                                                                    fit_output(i,j,k) = struct('Chisquare',output.ChiSquare,'search_params',[c1(i) c2(j) kp(k)],'condn',output.condn);
                                                                    spectra_output(i,j,k) = struct('magnetic_field',spectra.magnetic_fields,'aS',output.aS,'bS',output.bS,...
                                                                        'obs_cS',S_coex,'exp_cS',output.exp_cS);
                                                                    tieline_output(i,j,k) = struct('aC',output.aC,'bC',output.bC,'fpafpb',output.fpafpb,'fafb',output.fafb,'Kp',output.Kp);
                                                                end
                                                            end
                                                            pack;
                                                            disp(sprintf('%f%% done',(i/npositions1)*100));
                                                        end
                                                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                                                        % analyze results of search and fitting   
                                                        [i,j,k] = find(min(Chi2(:)) == Chi2);
                                                        if isvector(i) 
                                                            for ii = 1:length(i)
                                                                best_chisquare(ii,1) = Chi2(i(ii),j(ii),k(ii));
                                                                best_indices(ii,:) = [i(ii) j(ii) k(ii)];
                                                                best_params(ii,:) = fit_output(i(ii),j(ii),k(ii)).search_params;
                                                                % plot boundary and alpha and beta phase compositions
                                                                % and plot tielines.
                                                                plot_tielines(tieline_output(i(ii),j(ii),k(ii)).aC,tieline_output(i(ii),j(ii),k(ii)).bC,C_bdy);
                                                            end 
                                                        else
                                                            best_chisquare = Chi2(i,j,k);
                                                            best_indices = [i j k];
                                                            best_params = fit_output(i,j,k).search_params;
                                                            % plot boundary and alpha and beta phase compositions
                                                            plot_tielines(tieline_output(i,j,k).aC,tieline_output(i,j,k).bC,C_bdy);
                                                        end
                                                        display = [best_chisquare';best_indices';best_params'];
                                                        disp(sprintf('best chisquare = %f, with indices [%d %d %d] and parameters [%f %f %f]\n',display));
                                %                         disp(sprintf('best slope = %f => angle = %f\n',[best_params;slope2angle(best_params)]));
                                                        disp(sprintf('done!'));
                                                        return % fitting finished
                                                    otherwise
                                                        error('invalid response');
                                                end
                                            case 'function'
                                                disp(sprintf('starting search and fitting...'));
                                                options = optimset('MaxFunEvals',500,'TolFun',1e-8,'TolX',1e-8,'MaxIter',500,'LargeScale','on','LevenbergMarquardt','on','Display','off','Jacobian','off');
                                                data = reshape(S_coex,ncB*ncS,1);
                                                npositions1 = length(c1);
                                                npositions2 = length(c2);
                                                for i = 1:npositions1
                                                    for j = 1:npositions2
                                                        x0 = [c1(i) c2(j) 0]';
                                                        lb = [c1(i)-0.001 c2(j)-0.001 -Inf]';
                                                        ub = [c1(i)+0.001 c2(j)+0.001 Inf]';
                                                        [x,resnorm] = lsqcurvefit(@tangent_tieline_fit_fxn,x0,C_coex,data,lb,ub,options,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                                        [expcS,output] = tangent_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                                        % output =
                                                        % struct('aC',alphaC,'bC',betaC,'aS',alphaS,'bS',betaS,
                                                        % 'exp_cS',expcS,'fpafpb',fpafpb,'fafb',faf
                                                        % b,'Kp',Kp,'ChiSquare',ChiSquare,'condn',condn);
                                                        Chi2(i) = output.ChiSquare;
                                                        % pack up variables for output
                                                        fit_output(i) = struct('Chisquare',output.ChiSquare,'search_params',x','condn',output.condn);
                                                        spectra_output(i) = struct('magnetic_field',spectra.magnetic_fields,'aS',output.aS,'bS',output.bS,...
                                                            'obs_cS',S_coex,'exp_cS',output.exp_cS);
                                                        tieline_output(i) = struct('aC',output.aC,'bC',output.bC,'fpafpb',output.fpafpb,'fafb',output.fafb,'Kp',output.Kp);
                                                    end 
                                                    pack;
                                                    disp(sprintf('%f%% done',(i/npositions1)*100));
                                                end
                                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                                % analyze results of search and fitting   
                                                [i,j] = find(min(Chi2(:)) == Chi2);
                                                if isvector(i) 
                                                    for ii = 1:length(i)
                                                        best_chisquare(ii,1) = Chi2(i(ii),j(ii));
                                                        best_indices(ii,:) = [i(ii) j(ii)];
                                                        best_params(ii,:) = fit_output(i(ii),j(ii)).search_params;
                                                        % plot boundary and alpha and beta phase compositions
                                                        % and plot tielines.
                                                        plot_tielines(tieline_output(i(ii),j(ii)).aC,tieline_output(i(ii),j(ii)).bC,C_bdy);
                                                    end 
                                                else
                                                    best_chisquare = Chi2(i,j);
                                                    best_indices = [i j];
                                                    best_params = fit_output(i,j).search_params;
                                                    % plot boundary and alpha and beta phase compositions
                                                    plot_tielines(tieline_output(i,j).aC,tieline_output(i,j).bC,C_bdy);
                                                end
                                                display = [best_chisquare';best_indices';best_params'];
                                                disp(sprintf('best chisquare = %f, with indices [%d %d] and parameters [%f %f %f]\n',display));
                        %                         disp(sprintf('best slope = %f => angle = %f\n',[best_params;slope2angle(best_params)]));
                                                disp(sprintf('done!'));
                                                return % fitting finished
                                            otherwise
                                                error('invalid Kp configuration');
                                        end
                                    otherwise
                                        error('invalid boundary search method');
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
                % analyze boundary structure to organize search
                switch boundary.critpt1position
                    case 'located' 
                        switch boundary.critpt2position
                            case 'located'
                                if boundary.critpt1parameter == boundary.critpt2parameter
                                    error('both critical points cannot be the same point');
                                end
                                c1 = boundary.critpt1parameter;
                                c2 = boundary.critpt2parameter;
                                % perform search and fit
                                disp(sprintf('starting search and fitting...'));
                                options = optimset('MaxFunEvals',500,'TolFun',1e-8,'TolX',1e-8,'MaxIter',500,'LargeScale','on','LevenbergMarquardt','on','Display','iter','Jacobian','off');
                                data = reshape(S_coex,ncB*ncS,1);
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
                                        reply = input(sprintf('enter 4 parameter starting vector as [%4.3f<=cpt1<=%4.3f %4.3f<=cpt2<=%4.3f a>0 %4.2f<=kp<=%4.2f]:(enter for default)\n',[c1(1);c1(end);c2(1);c2(end);kp(1);kp(end)]));
                                        if isempty(reply)
                                            x0 = [c1(1)+(c1(end)-c1(1))*rand(1) c2(1)+(c2(end)-c2(1))*rand(1) 1 kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                        elseif isvector(reply)
                                            x0 = reply';
                                        else 
                                            error('invalid response');
                                        end
                                        lb = [c1-0.001 c2-0.001 1e-10 kp(1)]';
                                        ub = [c1+0.001 c2+0.001 Inf kp(end)]';
                                        [x,resnorm] = lsqcurvefit(@ruled_tieline_fit_fxn,x0,C_coex,data,lb,ub,options,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                    case 'function'
                                        reply = input(sprintf('enter 4 parameter starting vector as [%4.3f<=cpt1<=%4.3f %4.3f<=cpt2<=%4.3f a>0 c]:(enter for default)\n',[c1(1);c1(end);c2(1);c2(end)]));
                                        if isempty(reply)
                                            x0 = [c1(1)+(c1(end)-c1(1))*rand(1) c2(1)+(c2(end)-c2(1))*rand(1) 1 0]'; % for 2 critical points, Kp function has 1 parameter
                                        elseif isvector(reply)
                                            x0 = reply';
                                        else 
                                            error('invalid response');
                                        end
                                        lb = [c1-0.001 c2-0.001 1e-10 -Inf]';
                                        ub = [c1+0.001 c2+0.001 Inf Inf]';
                                        [x,resnorm] = lsqcurvefit(@ruled_tieline_fit_fxn,x0,C_coex,data,lb,ub,options,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                    otherwise
                                        error('invalid Kp configuration');
                                end
                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                [expcS,output] = ruled_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                plot_tielines(output.aC,output.bC,C_bdy);
                                display = [resnorm;x];
                                disp(sprintf('best chisquare = %f, with parameters [%f %f %f %f]\n',display));
                                % pack up output structures
                                fit_output = struct('Chisquare',output.ChiSquare,'search_params',x','condn',output.condn);
                                spectra_output = struct('magnetic_field',spectra.magnetic_fields,'aS',output.aS,'bS',output.bS,...
                                            'obs_cS',S_coex,'exp_cS',output.exp_cS);
                                tieline_output = struct('aC',output.aC,'bC',output.bC,'fpafpb',output.fpafpb,'fafb',output.fafb,'Kp',output.Kp);
                                disp(sprintf('done!'));
                                return % fitting finished         
                            case {'range','unknown'}
                                if length(boundary.critpt2parameter) == 2 % picked from figure as 2 points or unknown
                                    if boundary.critpt2parameter(1) > boundary.critpt2parameter(2)
                                        start = boundary.critpt2parameter(2);
                                        finish = boundary.critpt2parameter(1);
                                    elseif boundary.critpt2parameter(1) < boundary.critpt2parameter(2)
                                        start = boundary.critpt2parameter(1);
                                        finish = boundary.critpt2parameter(2);
                                    else
                                        error('critical point 2 range is same point');
                                    end
                                    boundary.critpt2parameter = linspace(start,finish,50*(finish-start));
                                end
                                c1 = boundary.critpt1parameter;
                                c2 = boundary.critpt2parameter;
                                % perform search and fit
                                switch boundary.search_method
                                    case 'continuous'
                                        disp(sprintf('starting search and fitting...'));
                                        options = optimset('MaxFunEvals',500,'TolFun',1e-8,'TolX',1e-8,'MaxIter',500,'LargeScale','on','LevenbergMarquardt','on','Display','iter','Jacobian','off');
                                        data = reshape(S_coex,ncB*ncS,1);
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
                                                reply = input(sprintf('enter 4 parameter starting vector as [%4.3f<=cpt1<=%4.3f %4.3f<=cpt2<=%4.3f a>0 %4.2f<=kp<=%4.2f]:(enter for default)\n',[c1(1);c1(end);c2(1);c2(end);kp(1);kp(end)]));
                                                if isempty(reply)
                                                    x0 = [c1(1)+(c1(end)-c1(1))*rand(1) c2(1)+(c2(end)-c2(1))*rand(1) 1 kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                                elseif isvector(reply)
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end
                                                lb = [c1-0.001 c2(1) 1e-10 kp(1)]';
                                                ub = [c1+0.001 c2(end) Inf kp(end)]';
                                                [x,resnorm] = lsqcurvefit(@ruled_tieline_fit_fxn,x0,C_coex,data,lb,ub,options,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                            case 'function'
                                                reply = input(sprintf('enter 4 parameter starting vector as [%4.3f<=cpt1<=%4.3f %4.3f<=cpt2<=%4.3f a>0 c]:(enter for default)\n',[c1(1);c1(end);c2(1);c2(end)]));
                                                if isempty(reply)
                                                    x0 = [c1(1)+(c1(end)-c1(1))*rand(1) c2(1)+(c2(end)-c2(1))*rand(1) 1 0]'; % for 2 critical points, Kp function has 1 parameter
                                                elseif isvector(reply)
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end
                                                lb = [c1-0.001 c2(1) 1e-10 -Inf]';
                                                ub = [c1+0.001 c2(end) Inf Inf]';
                                                [x,resnorm] = lsqcurvefit(@ruled_tieline_fit_fxn,x0,C_coex,data,lb,ub,options,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                            otherwise
                                                error('invalid Kp configuration');
                                        end
                                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                                        [expcS,output] = ruled_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                        plot_tielines(output.aC,output.bC,C_bdy);
                                        display = [resnorm;x];
                                        disp(sprintf('best chisquare = %f, with parameters [%f %f %f %f]\n',display));
                                        % pack up output structures
                                        fit_output = struct('Chisquare',output.ChiSquare,'search_params',x','condn',output.condn);
                                        spectra_output = struct('magnetic_field',spectra.magnetic_fields,'aS',output.aS,'bS',output.bS,...
                                                    'obs_cS',S_coex,'exp_cS',output.exp_cS);
                                        tieline_output = struct('aC',output.aC,'bC',output.bC,'fpafpb',output.fpafpb,'fafb',output.fafb,'Kp',output.Kp);
                                        disp(sprintf('done!'));
                                        return % fitting finished
                                    case 'grid'
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
                                                disp(sprintf('starting search and fitting...'));
                                                options = optimset('MaxFunEvals',500,'TolFun',1e-8,'TolX',1e-8,'MaxIter',500,'LargeScale','on','LevenbergMarquardt','on','Display','off','Jacobian','off');
                                                data = reshape(S_coex,ncB*ncS,1);
                                                npositions = length(c2);
                                                for i = 1:npositions
                                                    x0 = [c1 c2(i) 1 kp(1)+(kp(end)-kp(1)).*0.5]';% kp(1)+(kp(end)-kp(1)).*rand(1)
                                                    lb = [c1-0.001 c2(i)-0.001 1e-10 kp(1)]';
                                                    ub = [c1+0.001 c2(i)+0.001 Inf kp(end)]';
                                                    [x,resnorm] = lsqcurvefit(@ruled_tieline_fit_fxn,x0,C_coex,data,lb,ub,options,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                                    [expcS,output] = ruled_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                                    % output =
                                                    % struct('aC',alphaC,'bC',betaC,'aS',alphaS,'bS',betaS,
                                                    % 'exp_cS',expcS,'fpafpb',fpafpb,'fafb',faf
                                                    % b,'Kp',Kp,'ChiSquare',ChiSquare,'condn',condn);
                                                    Chi2(i) = output.ChiSquare;
                                                    % pack up variables for output
                                                    fit_output(i) = struct('Chisquare',output.ChiSquare,'search_params',x','condn',output.condn);
                                                    spectra_output(i) = struct('magnetic_field',spectra.magnetic_fields,'aS',output.aS,'bS',output.bS,...
                                                        'obs_cS',S_coex,'exp_cS',output.exp_cS);
                                                    tieline_output(i) = struct('aC',output.aC,'bC',output.bC,'fpafpb',output.fpafpb,'fafb',output.fafb,'Kp',output.Kp);
                                                    pack;
                                                    disp(sprintf('%f%% done',(i/npositions)*100));
                                                end
                                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                                % analyze results of search and fitting   
                                                i = find(min(Chi2) == Chi2);
                                                if isvector(i) 
                                                    for ii = 1:length(i)
                                                        best_chisquare(ii,1) = Chi2(i(ii));
                                                        best_indices(ii) = i(ii);
                                                        best_params(ii,:) = fit_output(i(ii)).search_params;
                                                        % plot boundary and alpha and beta phase compositions
                                                        % and plot tielines.
                                                        plot_tielines(tieline_output(i(ii)).aC,tieline_output(i(ii)).bC,C_bdy);
                                                    end 
                                                else
                                                    best_chisquare = Chi2(i);
                                                    best_indices = i;
                                                    best_params = fit_output(i).search_params;
                                                    % plot boundary and alpha and beta phase compositions
                                                    plot_tielines(tieline_output(i).aC,tieline_output(i).bC,C_bdy);
                                                end
                                                display = [best_chisquare';best_indices';best_params'];
                                                disp(sprintf('best chisquare = %f, with index %d and parameters [%f %f %f %f]\n',display));
                        %                         disp(sprintf('best slope = %f => angle = %f\n',[best_params;slope2angle(best_params)]));
                                                disp(sprintf('done!'));
                                                return % fitting finished
                                            case 'function'
                                                disp(sprintf('starting search and fitting...'));
                                                options = optimset('MaxFunEvals',500,'TolFun',1e-8,'TolX',1e-8,'MaxIter',500,'LargeScale','on','LevenbergMarquardt','on','Display','off','Jacobian','off');
                                                data = reshape(S_coex,ncB*ncS,1);
                                                npositions = length(c2);
                                                for i = 1:npositions
                                                    x0 = [c1 c2(i) 1 0]';
                                                    lb = [c1-0.001 c2(i)-0.001 1e-10 -Inf]';
                                                    ub = [c1+0.001 c2(i)+0.001 Inf Inf]';
                                                    [x,resnorm] = lsqcurvefit(@ruled_tieline_fit_fxn,x0,C_coex,data,lb,ub,options,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                                    [expcS,output] = ruled_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                                    % output =
                                                    % struct('aC',alphaC,'bC',betaC,'aS',alphaS,'bS',betaS,
                                                    % 'exp_cS',expcS,'fpafpb',fpafpb,'fafb',faf
                                                    % b,'Kp',Kp,'ChiSquare',ChiSquare,'condn',condn);
                                                    Chi2(i) = output.ChiSquare;
                                                    % pack up variables for output
                                                    fit_output(i) = struct('Chisquare',output.ChiSquare,'search_params',x','condn',output.condn);
                                                    spectra_output(i) = struct('magnetic_field',spectra.magnetic_fields,'aS',output.aS,'bS',output.bS,...
                                                        'obs_cS',S_coex,'exp_cS',output.exp_cS);
                                                    tieline_output(i) = struct('aC',output.aC,'bC',output.bC,'fpafpb',output.fpafpb,'fafb',output.fafb,'Kp',output.Kp);
                                                    pack;
                                                    disp(sprintf('%f%% done',(i/npositions)*100));
                                                end
                                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                                % analyze results of search and fitting   
                                                i = find(min(Chi2) == Chi2);
                                                if isvector(i) 
                                                    for ii = 1:length(i)
                                                        best_chisquare(ii,1) = Chi2(i(ii));
                                                        best_indices(ii) = i(ii);
                                                        best_params(ii,:) = fit_output(i(ii)).search_params;
                                                        % plot boundary and alpha and beta phase compositions
                                                        % and plot tielines.
                                                        plot_tielines(tieline_output(i(ii),k(ii)).aC,tieline_output(i(ii),k(ii)).bC,C_bdy);
                                                    end 
                                                else
                                                    best_chisquare = Chi2(i);
                                                    best_indices = i;
                                                    best_params = fit_output(i).search_params;
                                                    % plot boundary and alpha and beta phase compositions
                                                    plot_tielines(tieline_output(i,k).aC,tieline_output(i,k).bC,C_bdy);
                                                end
                                                display = [best_chisquare';best_indices';best_params'];
                                                disp(sprintf('best chisquare = %f, with index %d and parameters [%f %f %f %f]\n',display));
                        %                         disp(sprintf('best slope = %f => angle = %f\n',[best_params;slope2angle(best_params)]));
                                                disp(sprintf('done!'));
                                                return % fitting finished  
                                            otherwise
                                                error('invalid Kp configuration');
                                        end
                                    otherwise
                                        error('invalid boundary search method');
                                end
                            otherwise
                                error('invalid critical point 2 position on boundary');
                        end
                    case {'range','unknown'} % for critical point 1
                        if length(boundary.critpt1parameter) == 2 % picked from figure as 2 points or unknown
                            if boundary.critpt1parameter > boundary.critpt1parameter(2)
                                start = boundary.critpt1parameter(2);
                                finish = boundary.critpt1parameter(1);
                            elseif boundary.critpt1parameter(1) < boundary.critpt1parameter(2)
                                start = boundary.critpt1parameter(1);
                                finish = boundary.critpt1parameter(2);
                            else
                                error('critical point 1 range is same point');
                            end
                            boundary.critpt1parameter = linspace(start,finish,50*(finish-start));
                        end
                        switch boundary.critpt2position
                            case 'located'
%                                 if any(boundary.critpt1parameter == boundary.critpt2parameter)
%                                     error('both critical points cannot be the same point');
%                                 end
                                % perform search and fit
                                c1 = boundary.critpt1parameter;
                                c2 = boundary.critpt2parameter;
                                switch boundary.search_method
                                    case 'continuous'
                                        disp(sprintf('starting search and fitting...'));
                                        options = optimset('MaxFunEvals',500,'TolFun',1e-8,'TolX',1e-8,'MaxIter',500,'LargeScale','on','LevenbergMarquardt','on','Display','iter','Jacobian','off');
                                        data = reshape(S_coex,ncB*ncS,1);
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
                                                reply = input(sprintf('enter 4 parameter starting vector as [%4.3f<=cpt1<=%4.3f %4.3f<=cpt2<=%4.3f a>0 %4.2f<=kp<=%4.2f]:(enter for default)\n',[c1(1);c1(end);c2(1);c2(end);kp(1);kp(end)]));
                                                if isempty(reply)
                                                    x0 = [c1(1)+(c1(end)-c1(1))*rand(1) c2(1)+(c2(end)-c2(1))*rand(1) 1 kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                                elseif isvector(reply)
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end
                                                lb = [c1(1) c2-0.001 1e-10 kp(1)]';
                                                ub = [c1(end) c2+0.001 Inf kp(end)]';
                                                [x,resnorm] = lsqcurvefit(@ruled_tieline_fit_fxn,x0,C_coex,data,lb,ub,options,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                            case 'function'
                                                reply = input(sprintf('enter 4 parameter starting vector as [%4.3f<=cpt1<=%4.3f %4.3f<=cpt2<=%4.3f a>0 c]:(enter for default)\n',[c1(1);c1(end);c2(1);c2(end)]));
                                                if isempty(reply)
                                                    x0 = [c1(1)+(c1(end)-c1(1))*rand(1) c2(1)+(c2(end)-c2(1))*rand(1) 1 0]'; % for 2 critical points, Kp function has 1 parameter
                                                elseif isvector(reply)
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end
                                                lb = [c1(1) c2-0.001 1e-10 -Inf]';
                                                ub = [c1(end) c2+0.001 Inf Inf]';
                                                [x,resnorm] = lsqcurvefit(@ruled_tieline_fit_fxn,x0,C_coex,data,lb,ub,options,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                            otherwise
                                                error('invalid Kp configuration');
                                        end
                                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                                        [expcS,output] = ruled_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                        plot_tielines(output.aC,output.bC,C_bdy);
                                        display = [resnorm;x];
                                        disp(sprintf('best chisquare = %f, with parameters [%f %f %f %f]\n',display));
                                        % pack up output structures
                                        fit_output = struct('Chisquare',output.ChiSquare,'search_params',x','condn',output.condn);
                                        spectra_output = struct('magnetic_field',spectra.magnetic_fields,'aS',output.aS,'bS',output.bS,...
                                                    'obs_cS',S_coex,'exp_cS',output.exp_cS);
                                        tieline_output = struct('aC',output.aC,'bC',output.bC,'fpafpb',output.fpafpb,'fafb',output.fafb,'Kp',output.Kp);
                                        disp(sprintf('done!'));
                                        return % fitting finished
                                    case 'grid'
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
                                                disp(sprintf('starting search and fitting...'));
                                                options = optimset('MaxFunEvals',500,'TolFun',1e-8,'TolX',1e-8,'MaxIter',500,'LargeScale','on','LevenbergMarquardt','on','Display','off','Jacobian','off');
                                                data = reshape(S_coex,ncB*ncS,1);
                                                npositions = length(c1);
                                                for i = 1:npositions
                                                    x0 = [c1(i) c2 1 kp(1)+(kp(end)-kp(1)).*0.5]';% kp(1)+(kp(end)-kp(1)).*rand(1)
                                                    lb = [c1(i)-0.001 c2-0.001 1e-10 kp(1)]';
                                                    ub = [c1(i)+0.001 c2+0.001 Inf kp(end)]';
                                                    [x,resnorm] = lsqcurvefit(@ruled_tieline_fit_fxn,x0,C_coex,data,lb,ub,options,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                                    [expcS,output] = ruled_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                                    % output =
                                                    % struct('aC',alphaC,'bC',betaC,'aS',alphaS,'bS',betaS,
                                                    % 'exp_cS',expcS,'fpafpb',fpafpb,'fafb',faf
                                                    % b,'Kp',Kp,'ChiSquare',ChiSquare,'condn',condn);
                                                    Chi2(i) = output.ChiSquare;
                                                    % pack up variables for output
                                                    fit_output(i) = struct('Chisquare',output.ChiSquare,'search_params',x','condn',output.condn);
                                                    spectra_output(i) = struct('magnetic_field',spectra.magnetic_fields,'aS',output.aS,'bS',output.bS,...
                                                        'obs_cS',S_coex,'exp_cS',output.exp_cS);
                                                    tieline_output(i) = struct('aC',output.aC,'bC',output.bC,'fpafpb',output.fpafpb,'fafb',output.fafb,'Kp',output.Kp);
                                                    pack;
                                                    disp(sprintf('%f%% done',(i/npositions)*100));
                                                end
                                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                                % analyze results of search and fitting   
                                                i = find(min(Chi2) == Chi2);
                                                if isvector(i) 
                                                    for ii = 1:length(i)
                                                        best_chisquare(ii,1) = Chi2(i(ii));
                                                        best_indices(ii) = i(ii);
                                                        best_params(ii,:) = fit_output(i(ii)).search_params;
                                                        % plot boundary and alpha and beta phase compositions
                                                        % and plot tielines.
                                                        plot_tielines(tieline_output(i(ii)).aC,tieline_output(i(ii)).bC,C_bdy);
                                                    end 
                                                else
                                                    best_chisquare = Chi2(i);
                                                    best_indices = i;
                                                    best_params = fit_output(i).search_params;
                                                    % plot boundary and alpha and beta phase compositions
                                                    plot_tielines(tieline_output(i).aC,tieline_output(i).bC,C_bdy);
                                                end
                                                display = [best_chisquare';best_indices';best_params'];
                                                disp(sprintf('best chisquare = %f, with index %d and parameters [%f %f %f %f]\n',display));
                        %                         disp(sprintf('best slope = %f => angle = %f\n',[best_params;slope2angle(best_params)]));
                                                disp(sprintf('done!'));
                                                return % fitting finished
                                            case 'function'
                                                disp(sprintf('starting search and fitting...'));
                                                options = optimset('MaxFunEvals',500,'TolFun',1e-8,'TolX',1e-8,'MaxIter',500,'LargeScale','on','LevenbergMarquardt','on','Display','off','Jacobian','off');
                                                data = reshape(S_coex,ncB*ncS,1);
                                                npositions = length(c2);
                                                for i = 1:npositions
                                                    x0 = [c1(i) c2 1 0]';
                                                    lb = [c1(i)-0.001 c2-0.001 1e-10 -Inf]';
                                                    ub = [c1(i)+0.001 c2+0.001 Inf Inf]';
                                                    [x,resnorm] = lsqcurvefit(@ruled_tieline_fit_fxn,x0,C_coex,data,lb,ub,options,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                                    [expcS,output] = ruled_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                                    % output =
                                                    % struct('aC',alphaC,'bC',betaC,'aS',alphaS,'bS',betaS,
                                                    % 'exp_cS',expcS,'fpafpb',fpafpb,'fafb',faf
                                                    % b,'Kp',Kp,'ChiSquare',ChiSquare,'condn',condn);
                                                    Chi2(i) = output.ChiSquare;
                                                    % pack up variables for output
                                                    fit_output(i) = struct('Chisquare',output.ChiSquare,'search_params',x','condn',output.condn);
                                                    spectra_output(i) = struct('magnetic_field',spectra.magnetic_fields,'aS',output.aS,'bS',output.bS,...
                                                        'obs_cS',S_coex,'exp_cS',output.exp_cS);
                                                    tieline_output(i) = struct('aC',output.aC,'bC',output.bC,'fpafpb',output.fpafpb,'fafb',output.fafb,'Kp',output.Kp);
                                                    pack;
                                                    disp(sprintf('%f%% done',(i/npositions)*100));
                                                end
                                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                                % analyze results of search and fitting   
                                                i = find(min(Chi2) == Chi2);
                                                if isvector(i) 
                                                    for ii = 1:length(i)
                                                        best_chisquare(ii,1) = Chi2(i(ii));
                                                        best_indices(ii) = i(ii);
                                                        best_params(ii,:) = fit_output(i(ii)).search_params;
                                                        % plot boundary and alpha and beta phase compositions
                                                        % and plot tielines.
                                                        plot_tielines(tieline_output(i(ii)).aC,tieline_output(i(ii)).bC,C_bdy);
                                                    end 
                                                else
                                                    best_chisquare = Chi2(i);
                                                    best_indices = i;
                                                    best_params = fit_output(i).search_params;
                                                    % plot boundary and alpha and beta phase compositions
                                                    plot_tielines(tieline_output(i).aC,tieline_output(i).bC,C_bdy);
                                                end
                                                display = [best_chisquare';best_indices';best_params'];
                                                disp(sprintf('best chisquare = %f, with indices %d and parameters [%f %f %f %f]\n',display));
                        %                         disp(sprintf('best slope = %f => angle = %f\n',[best_params;slope2angle(best_params)]));
                                                disp(sprintf('done!'));
                                                return % fitting finished  
                                            otherwise
                                                error('invalid Kp configuration');
                                        end
                                    otherwise
                                        error('invalid boundary search method');
                                end
                            case {'range','unknown'} % for critical point 2, critical point 1 also range or unknown
                                if length(boundary.critpt2parameter) == 2 % picked from figure as 2 points
                                    if boundary.critpt2parameter(1) > boundary.critpt2parameter(2)
                                        start = boundary.critpt2parameter(2);
                                        finish = boundary.critpt2parameter(1);
                                    elseif boundary.critpt2parameter(1) < boundary.critpt2parameter(2)
                                        start = boundary.critpt2parameter(1);
                                        finish = boundary.critpt2parameter(2);
                                    else
                                        error('critical point 2 range is same point');
                                    end
                                    boundary.critpt2parameter = linspace(start,finish,50*(finish-start));
                                end
                                % perform search and fit
                                c1 = boundary.critpt1parameter;
                                c2 = boundary.critpt2parameter;
                                switch boundary.search_method
                                    case 'continuous'
                                        disp(sprintf('starting search and fitting...'));
                                        options = optimset('OutputFcn',@outfun,'MaxFunEvals',2500,'TolFun',1e-8,'TolX',1e-8,'MaxIter',500,'LargeScale','on','LevenbergMarquardt','on','Display','iter','Jacobian','off');
                                        data = reshape(S_coex,ncB*ncS,1);
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
                                                reply = input(sprintf('enter 4 parameter starting vector as [%4.3f<=cpt1<=%4.3f %4.3f<=cpt2<=%4.3f a>0 %4.2f<=kp<=%4.2f]:(enter for default)\n',[c1(1);c1(end);c2(1);c2(end);kp(1);kp(end)]));
                                                if isempty(reply)
                                                    x0 = [c1(1)+(c1(end)-c1(1))*rand(1) c2(1)+(c2(end)-c2(1))*rand(1) 1 kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                                elseif isvector(reply)
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end
                                                lb = [c1(1) c2(1) 1e-10 kp(1)]';
                                                ub = [c1(end) c2(end) Inf kp(end)]';
                                                [x,resnorm] = lsqcurvefit(@ruled_tieline_fit_fxn,x0,C_coex,data,lb,ub,options,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                            case 'function'
                                                reply = input(sprintf('enter 4 parameter starting vector as [%4.3f<=cpt1<=%4.3f %4.3f<=cpt2<=%4.3f a>0 c]:(enter for default)\n',[c1(1);c1(end);c2(1);c2(end)]));
                                                if isempty(reply)
                                                    x0 = [c1(1)+(c1(end)-c1(1))*rand(1) c2(1)+(c2(end)-c2(1))*rand(1) 1 0]'; % for 2 critical points, Kp function has 1 parameter
                                                elseif isvector(reply)
                                                    x0 = reply';
                                                else 
                                                    error('invalid response');
                                                end
                                                lb = [c1(1) c2(1) 1e-10 -Inf]';
                                                ub = [c1(end) c2(end) Inf Inf]';
                                                [x,resnorm] = lsqcurvefit(@ruled_tieline_fit_fxn,x0,C_coex,data,lb,ub,options,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                            otherwise
                                                error('invalid Kp configuration');
                                        end
                                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                                        [expcS,output] = ruled_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                        plot_tielines(output.aC,output.bC,C_bdy);
                                        display = [resnorm;x];
                                        disp(sprintf('best chisquare = %f, with parameters [%f %f %f %f]\n',display));
                                        % pack up output structures
                                        fit_output = struct('Chisquare',output.ChiSquare,'search_params',x','condn',output.condn);
                                        spectra_output = struct('magnetic_field',spectra.magnetic_fields,'aS',output.aS,'bS',output.bS,...
                                                    'obs_cS',S_coex,'exp_cS',output.exp_cS);
                                        tieline_output = struct('aC',output.aC,'bC',output.bC,'fpafpb',output.fpafpb,'fafb',output.fafb,'Kp',output.Kp);
                                        disp(sprintf('done!'));
                                        return % fitting finished
                                    case 'grid'
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
                                                disp(sprintf('starting search and fitting...'));
                                                options = optimset('MaxFunEvals',500,'TolFun',1e-8,'TolX',1e-8,'MaxIter',500,'LargeScale','on','LevenbergMarquardt','on','Display','off','Jacobian','off');
                                                data = reshape(S_coex,ncB*ncS,1);
                                                npositions1 = length(c1);
                                                npositions2 = length(c2);
                                                for i = 1:npositions1
                                                    for j = 1:npositions2
                                                        x0 = [c1(i) c2(j) 1 kp(1)+(kp(end)-kp(1)).*0.5]';% kp(1)+(kp(end)-kp(1)).*rand(1)
                                                        lb = [c1(i)-0.001 c2(j)-0.001 1e-10 kp(1)]';
                                                        ub = [c1(i)+0.001 c2(j)+0.001 Inf kp(end)]';
                                                        [x,resnorm] = lsqcurvefit(@ruled_tieline_fit_fxn,x0,C_coex,data,lb,ub,options,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                                        [expcS,output] = ruled_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                                        % output =
                                                        % struct('aC',alphaC,'bC',betaC,'aS',alphaS,'bS',betaS,
                                                        % 'exp_cS',expcS,'fpafpb',fpafpb,'fafb',faf
                                                        % b,'Kp',Kp,'ChiSquare',ChiSquare,'condn',condn);
                                                        Chi2(i,j) = output.ChiSquare;
                                                        % pack up variables for output
                                                        fit_output(i,j) = struct('Chisquare',output.ChiSquare,'search_params',x','condn',output.condn);
                                                        spectra_output(i,j) = struct('magnetic_field',spectra.magnetic_fields,'aS',output.aS,'bS',output.bS,...
                                                            'obs_cS',S_coex,'exp_cS',output.exp_cS);
                                                        tieline_output(i,j) = struct('aC',output.aC,'bC',output.bC,'fpafpb',output.fpafpb,'fafb',output.fafb,'Kp',output.Kp);
                                                    end
                                                    pack;
                                                    disp(sprintf('%f%% done',(i/npositions1)*100));
                                                end
                                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                                % analyze results of search and fitting   
                                                [i,j] = find(min(Chi2(:)) == Chi2);
                                                if isvector(i) 
                                                    for ii = 1:length(i)
                                                        best_chisquare(ii,1) = Chi2(i(ii),j(ii));
                                                        best_indices(ii,:) = [i(ii) j(ii)];
                                                        best_params(ii,:) = fit_output(i(ii),j(ii)).search_params;
                                                        % plot boundary and alpha and beta phase compositions
                                                        % and plot tielines.
                                                        plot_tielines(tieline_output(i(ii),j(ii)).aC,tieline_output(i(ii),j(ii)).bC,C_bdy);
                                                    end 
                                                else
                                                    best_chisquare = Chi2(i,j);
                                                    best_indices = [i j];
                                                    best_params = fit_output(i,j).search_params;
                                                    % plot boundary and alpha and beta phase compositions
                                                    plot_tielines(tieline_output(i,j).aC,tieline_output(i,j).bC,C_bdy);
                                                end
                                                display = [best_chisquare';best_indices';best_params'];
                                                disp(sprintf('best chisquare = %f, with indices [%d %d] and parameters [%f %f %f %f]\n',display));
                        %                         disp(sprintf('best slope = %f => angle = %f\n',[best_params;slope2angle(best_params)]));
                                                disp(sprintf('done!'));
                                                return % fitting finished
                                            case 'function'
                                                disp(sprintf('starting search and fitting...'));
                                                options = optimset('MaxFunEvals',500,'TolFun',1e-8,'TolX',1e-8,'MaxIter',500,'LargeScale','on','LevenbergMarquardt','on','Display','off','Jacobian','off');
                                                data = reshape(S_coex,ncB*ncS,1);
                                                npositions1 = length(c1);
                                                npositions2 = length(c2);
                                                for i = 1:npositions1
                                                    for j = 1:npositions2
                                                        x0 = [c1(i) c2(j) 1 0]';
                                                        lb = [c1(i)-0.001 c2(j)-0.001 1e-10 -Inf]';
                                                        ub = [c1(i)+0.001 c2(j)+0.001 Inf Inf]';
                                                        [x,resnorm] = lsqcurvefit(@ruled_tieline_fit_fxn,x0,C_coex,data,lb,ub,options,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                                        [expcS,output] = ruled_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,tieline,part_coeff,spectra);
                                                        % output =
                                                        % struct('aC',alphaC,'bC',betaC,'aS',alphaS,'bS',betaS,
                                                        % 'exp_cS',expcS,'fpafpb',fpafpb,'fafb',faf
                                                        % b,'Kp',Kp,'ChiSquare',ChiSquare,'condn',condn);
                                                        Chi2(i) = output.ChiSquare;
                                                        % pack up variables for output
                                                        fit_output(i) = struct('Chisquare',output.ChiSquare,'search_params',x','condn',output.condn);
                                                        spectra_output(i) = struct('magnetic_field',spectra.magnetic_fields,'aS',output.aS,'bS',output.bS,...
                                                            'obs_cS',S_coex,'exp_cS',output.exp_cS);
                                                        tieline_output(i) = struct('aC',output.aC,'bC',output.bC,'fpafpb',output.fpafpb,'fafb',output.fafb,'Kp',output.Kp);
                                                    end 
                                                    pack;
                                                    disp(sprintf('%f%% done',(i/npositions1)*100));
                                                end
                                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                                % analyze results of search and fitting   
                                                [i,j] = find(min(Chi2(:)) == Chi2);
                                                if isvector(i) 
                                                    for ii = 1:length(i)
                                                        best_chisquare(ii,1) = Chi2(i(ii),j(ii));
                                                        best_indices(ii,:) = [i(ii) j(ii)];
                                                        best_params(ii,:) = fit_output(i(ii),j(ii)).search_params;
                                                        % plot boundary and alpha and beta phase compositions
                                                        % and plot tielines.
                                                        plot_tielines(tieline_output(i(ii),j(ii)).aC,tieline_output(i(ii),j(ii)).bC,C_bdy);
                                                    end 
                                                else
                                                    best_chisquare = Chi2(i,j);
                                                    best_indices = [i j];
                                                    best_params = fit_output(i,j).search_params;
                                                    % plot boundary and alpha and beta phase compositions
                                                    plot_tielines(tieline_output(i,j).aC,tieline_output(i,j).bC,C_bdy);
                                                end
                                                display = [best_chisquare';best_indices';best_params'];
                                                disp(sprintf('best chisquare = %f, with indices [%d %d] and parameters [%f %f %f %f]\n',display));
                        %                         disp(sprintf('best slope = %f => angle = %f\n',[best_params;slope2angle(best_params)]));
                                                disp(sprintf('done!'));
                                                return % fitting finished
                                            otherwise
                                                error('invalid Kp configuration');
                                        end
                                    otherwise
                                        error('invalid boundary search method');
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
        disp(sprintf('with one critical point, there must be one end tieline with two end points:\n'));
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
                    reply = input(sprintf('what is the boundary parameter range to search for critical point %d? [0:1] (if unknown press enter)\n',i));
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
        end % end of loop over number critical points
         
%         reply = input(sprintf('what are the boundary parameters of the end tieline? [0<=#<=1 0<=#<=1] (if unknown press enter)\n'));
%         if isempty(reply)
%             disp(sprintf('click the 2 points on the boundary of the end tieline:\n'));
%             ternary_plot(C_bdy,'-k','linewidth',3);
%             [x,y] = ginput(2);
%             pause(3);
%             close;
%             p = bdypt2b([x y],b,C_bdy);
%             endpts = b2bdypt(b,C_bdy,p);
%             boundary.('endtl') = endpts;
%             boundary.(['endtlparameter']) = p;
%         elseif isvector(reply)
%             endpts = b2bdypt(b,C_bdy,reply);
%             boundary.('endtl') = endpts;
%             boundary.(['endtlparameter']) = reply;
%         else
%             error('invalid response');
%         end
        
        % reparameterize boundary so that zero is at one end point
        % of the end tieline and one is at the other end point.
        cp = boundary.critptparameter;
        e1 = boundary.endpt1parameter;
        e2 = boundary.endpt2parameter;
        
        switch boundary.critptposition
            case 'located'
                if e1 > e2
                    % divide boundary into two sections (same procedure
                    % when going from tangent to ruled in the 2
                    % critical point case above).
                    [u,directrixO,v,directrixA] = tangent2ruled(e1,e2,C_bdy);
                    if cp > e2 & cp < e1
                        Cbdy = directrixO;
                        bdy = [Cbdy;Cbdy(1,:)];
                        boundary.parameterization = u;
                    elseif cp < e2 | cp > e1
                        Cbdy = directrixA;
                        bdy = [Cbdy;Cbdy(1,:)];
                        boundary.parameterization = v;
                    else
                        error('cp cannot equal e1 or e2');
                    end
                    bi = bdypt2b(Cbdy,b,C_bdy);
                    Sbdy = interp1(b,S_bdy',bi);
                    Sbdy = Sbdy';
                    B = repmat(spectra.magnetic_fields,1,size(Sbdy,2));
                    temp = zeros(length(spectra.magnetic_fields),size(Sbdy,2)*2);
                    temp(:,1:2:end) = B;
                    temp(:,2:2:end) = Sbdy;
                    temp = normalize_spectra(temp,'d');
                    Sbdy = temp(:,2:2:end);
                    [IN ON] = inpolygon(C_coex(:,1),C_coex(:,2),bdy(:,1),bdy(:,2));
                    Ccoex = C_coex(IN,:);
                    Scoex = S_coex(:,IN);
                    boundary.critptparameter = bdypt2b(boundary.critpt,boundary.parameterization,Cbdy);
%                         ion = find(ON);
%                         inon = find(~ON);
                    clear temp B;
                elseif e1 < e2
                    % divide boundary into two sections (same procedure
                    % when going from tangent to ruled in the 2
                    % critical point case above).
                    [u,directrixO,v,directrixA] = tangent2ruled(e1,e2,C_bdy);
                    if cp > e1 & cp < e2
                        Cbdy = directrixA;
                        bdy = [Cbdy;Cbdy(1,:)];
                        boundary.parameterization = v;
                    elseif cp < e1 | cp > e2
                        Cbdy = directrixO;
                        bdy = [Cbdy;Cbdy(1,:)];
                        boundary.parameterization = u;
                    else
                        error('cp cannot equal e1 or e2');
                    end
                    bi = bdypt2b(Cbdy,b,C_bdy);
                    Sbdy = interp1(b,S_bdy',bi);
                    Sbdy = Sbdy';
                    B = repmat(spectra.magnetic_fields,1,size(Sbdy,2));
                    temp = zeros(length(spectra.magnetic_fields),size(Sbdy,2)*2);
                    temp(:,1:2:end) = B;
                    temp(:,2:2:end) = Sbdy;
                    temp = normalize_spectra(temp,'d');
                    Sbdy = temp(:,2:2:end);
                    [IN ON] = inpolygon(C_coex(:,1),C_coex(:,2),bdy(:,1),bdy(:,2));
                    Ccoex = C_coex(IN,:);
                    Scoex = S_coex(:,IN);
                    boundary.critptparameter = bdypt2b(boundary.critpt,boundary.parameterization,Cbdy);
%                         ion = find(ON);
%                         inon = find(~ON);
                    clear temp B;
                end
            case 'range'
                if cp(1) > cp(end)
                    start = cp(end);
                    finish = cp(1);
                elseif cp(1) < cp(end)
                    start = cp(1);
                    finish = cp(end);
                else
                    error('critical point range is same point');
                end
                if e1 > e2
                    % divide boundary into two sections (same procedure
                    % when going from tangent to ruled in the 2
                    % critical point case above).
                    [u,directrixO,v,directrixA] = tangent2ruled(e1,e2,C_bdy);
                    if all([start finish] > e2) & all([start finish] < e1)
                        Cbdy = directrixO;
                        bdy = [Cbdy;Cbdy(1,:)];
                        boundary.parameterization = u; 
                    elseif all([start finish] < e2) | all([start finish] > e1)
                        Cbdy = directrixA;
                        bdy = [Cbdy;Cbdy(1,:)];
                        boundary.parameterization = v;
                    else
                        error('cp cannot equal e1 or e2');
                    end
                    bi = bdypt2b(Cbdy,b,C_bdy);
                    Sbdy = interp1(b,S_bdy',bi);
                    Sbdy = Sbdy';
                    B = repmat(spectra.magnetic_fields,1,size(Sbdy,2));
                    temp = zeros(length(spectra.magnetic_fields),size(Sbdy,2)*2);
                    temp(:,1:2:end) = B;
                    temp(:,2:2:end) = Sbdy;
                    temp = normalize_spectra(temp,'d');
                    Sbdy = temp(:,2:2:end);
                    [IN ON] = inpolygon(C_coex(:,1),C_coex(:,2),bdy(:,1),bdy(:,2));
                    Ccoex = C_coex(IN,:);
                    Scoex = S_coex(:,IN);
                    boundary.critptparameter = bdypt2b(boundary.critpt,boundary.parameterization,Cbdy);
%                         ion = find(ON);
%                         inon = find(~ON);
                    clear temp B;
                elseif e1 < e2
                    % divide boundary into two sections (same procedure
                    % when going from tangent to ruled in the 2
                    % critical point case above).
                    [u,directrixO,v,directrixA] = tangent2ruled(e1,e2,C_bdy);
                    if all([start finish] > e1) & all([start finish] < e2) 
                        Cbdy = directrixA;
                        bdy = [Cbdy;Cbdy(1,:)];
                        boundary.parameterization = v;
                    elseif all([start finish] < e1) | all([start finish] > e2)  
                        Cbdy = directrixO;
                        bdy = [Cbdy;Cbdy(1,:)];
                        boundary.parameterization = u;
                    else
                        error('cp cannot equal e1 or e2');
                    end
                    bi = bdypt2b(Cbdy,b,C_bdy);
                    Sbdy = interp1(b,S_bdy',bi);
                    Sbdy = Sbdy';
                    B = repmat(spectra.magnetic_fields,1,size(Sbdy,2));
                    temp = zeros(length(spectra.magnetic_fields),size(Sbdy,2)*2);
                    temp(:,1:2:end) = B;
                    temp(:,2:2:end) = Sbdy;
                    temp = normalize_spectra(temp,'d');
                    Sbdy = temp(:,2:2:end);
                    [IN ON] = inpolygon(C_coex(:,1),C_coex(:,2),bdy(:,1),bdy(:,2));
                    Ccoex = C_coex(IN,:);
                    Scoex = S_coex(:,IN);
                    boundary.critptparameter = bdypt2b(boundary.critpt,boundary.parameterization,Cbdy);
%                         ion = find(ON);
%                         inon = find(~ON);
                    clear temp B;
                end
            case 'unknown'
                if cp(1) > cp(end)
                    start = cp(end);
                    finish = cp(1);
                elseif cp(1) < cp(end)
                    start = cp(1);
                    finish = cp(end);
                else
                    error('critical point range is same point');
                end
                [u,directrixO,v,directrixA] = tangent2ruled(e1,e2,C_bdy);
                if size(directrixO,1) >= size(directrixA,1)
                    Cbdy = directrixO;
                else
                    Cbdy = directrixA;
                end
                bdy = [Cbdy;Cbdy(1,:)];
                boundary.parameterization = v;
                bi = bdypt2b(Cbdy,b,C_bdy);
                Sbdy = interp1(b,S_bdy',bi);
                Sbdy = Sbdy';
                B = repmat(spectra.magnetic_fields,1,size(Sbdy,2));
                temp = zeros(length(spectra.magnetic_fields),size(Sbdy,2)*2);
                temp(:,1:2:end) = B;
                temp(:,2:2:end) = Sbdy;
                temp = normalize_spectra(temp,'d');
                Sbdy = temp(:,2:2:end);
                [IN ON] = inpolygon(C_coex(:,1),C_coex(:,2),bdy(:,1),bdy(:,2));
                Ccoex = C_coex(IN,:);
                Scoex = S_coex(:,IN);
%                         ion = find(ON);
%                         inon = find(~ON);
                clear temp B;
            otherwise
                error('invalid critical point position');
        end
        [ncC,ncd] = size(Ccoex);
        [ncB,ncS] = size(Scoex);
        [nbC,nbd] = size(Cbdy);
        [nbB,nbS] = size(Sbdy);
        switch tieline.config
            case 'parallel'
                
            case 'tangent'
                
            case 'ruled'
                switch boundary.critptposition
                    case 'located'
                        cp = boundary.critptparameter;
                        % perform search and fit
                        disp(sprintf('starting search and fitting...'));
                        options = optimset('MaxFunEvals',500,'TolFun',1e-8,'TolX',1e-8,'MaxIter',500,'LargeScale','on','LevenbergMarquardt','on','Display','iter','Jacobian','off');
                        data = reshape(Scoex,ncB*ncS,1);
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
                                reply = input(sprintf('enter 3 parameter starting vector as [%4.3f a>0 %4.2f<=kp<=%4.2f]:(enter for default)\n',[cp;kp(1);kp(end)]));
                                if isempty(reply)
                                    x0 = [cp 1 kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                elseif isvector(reply)
                                    x0 = reply';
                                else 
                                    error('invalid response');
                                end
                                lb = [cp-0.001 1e-10 kp(1)]';
                                ub = [cp+0.001 Inf kp(end)]';
                                [x,resnorm] = lsqcurvefit(@ruled_tieline_fit_fxn,x0,Ccoex,data,lb,ub,options,Scoex,Cbdy,Sbdy,boundary,tieline,part_coeff,spectra);
                            case 'function'
                                reply = input(sprintf('enter 4 parameter starting vector as [%4.3f a>0 c1 c2]:(enter for default)\n',cp));
                                if isempty(reply)
                                    x0 = [cp 1 0 0]'; % for 2 critical points, Kp function has 1 parameter
                                elseif isvector(reply)
                                    x0 = reply';
                                else 
                                    error('invalid response');
                                end
                                lb = [cp-0.001 1e-10 -Inf -Inf]';
                                ub = [cp+0.001 Inf Inf Inf]';
                                [x,resnorm] = lsqcurvefit(@ruled_tieline_fit_fxn,x0,Ccoex,data,lb,ub,options,Scoex,Cbdy,Sbdy,boundary,tieline,part_coeff,spectra);
                            otherwise
                                error('invalid Kp configuration');
                        end
                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                        [expcS,output] = ruled_tieline_fit_fxn(x,Ccoex,Scoex,Cbdy,Sbdy,boundary,tieline,part_coeff,spectra);
                        plot_tielines(output.aC,output.bC,Cbdy);
                        display = [resnorm;x];
                        disp(sprintf('best chisquare = %f, with parameters [%f %f %f %f]\n',display));
                        % pack up output structures
                        fit_output = struct('Chisquare',output.ChiSquare,'search_params',x','condn',output.condn);
                        spectra_output = struct('magnetic_field',spectra.magnetic_fields,'aS',output.aS,'bS',output.bS,...
                                    'obs_cS',Scoex,'exp_cS',output.exp_cS);
                        tieline_output = struct('aC',output.aC,'bC',output.bC,'fpafpb',output.fpafpb,'fafb',output.fafb,'Kp',output.Kp);
                        disp(sprintf('done!'));
                        return % fitting finished         
                    case {'range','unknown'}
                        cp = boundary.critptparameter;
                        if cp(1) > cp(end)
                            start = cp(end);
                            finish = cp(1);
                        elseif cp(1) < cp(end)
                            start = cp(1);
                            finish = cp(end);
                        else
                            error('critical point range is same point');
                        end
                        cp = [start finish];
                        % perform search and fit
                        switch boundary.search_method
                            case 'continuous'
                                disp(sprintf('starting search and fitting...'));
                                options = optimset('OutputFcn',@outfun,'MaxFunEvals',2500,'TolFun',1e-8,'TolX',1e-8,'MaxIter',500,'LargeScale','on','LevenbergMarquardt','on','Display','iter','Jacobian','off');
                                data = reshape(Scoex,ncB*ncS,1);
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
                                        reply = input(sprintf('enter 3 parameter starting vector as [%4.3f<=cpt<=%4.3f a>0 %4.2f<=kp<=%4.2f]:(enter for default)\n',[cp(1);cp(end);kp(1);kp(end)]));
                                        if isempty(reply)
                                            x0 = [cp(1)+(cp(end)-cp(1))*rand(1) 1 kp(1)+(kp(end)-kp(1)).*rand(1)]'; 
                                        elseif isvector(reply)
                                            x0 = reply';
                                        else 
                                            error('invalid response');
                                        end
                                        lb = [cp(1) 1e-10 kp(1)]';
                                        ub = [cp(end) Inf kp(end)]';
                                        [x,resnorm] = lsqcurvefit(@ruled_tieline_fit_fxn,x0,Ccoex,data,lb,ub,options,Scoex,Cbdy,Sbdy,boundary,tieline,part_coeff,spectra);
                                    case 'function'
                                        reply = input(sprintf('enter 4 parameter starting vector as [%4.3f<=cpt<=%4.3f a>0 c1 c2]:(enter for default)\n',[cp(1);cp(end)]));
                                        if isempty(reply)
                                            x0 = [cp(1)+(cp(end)-cp(1))*rand(1) 1 0 0]'; % for 2 critical points, Kp function has 1 parameter
                                        elseif isvector(reply)
                                            x0 = reply';
                                        else 
                                            error('invalid response');
                                        end
                                        lb = [cp(1) 1e-10 -Inf -Inf]';
                                        ub = [cp(end) Inf Inf Inf]';
                                        [x,resnorm] = lsqcurvefit(@ruled_tieline_fit_fxn,x0,Ccoex,data,lb,ub,options,Scoex,Cbdy,Sbdy,boundary,tieline,part_coeff,spectra);
                                    otherwise
                                        error('invalid Kp configuration');
                                end
                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                [expcS,output] = ruled_tieline_fit_fxn(x,Ccoex,Scoex,Cbdy,Sbdy,boundary,tieline,part_coeff,spectra);
                                plot_tielines(output.aC,output.bC,Cbdy);
                                display = [resnorm;x];
                                disp(sprintf('best chisquare = %f, with parameters [%f %f %f %f]\n',display));
                                % pack up output structures
                                fit_output = struct('Chisquare',output.ChiSquare,'search_params',x','condn',output.condn);
                                spectra_output = struct('magnetic_field',spectra.magnetic_fields,'aS',output.aS,'bS',output.bS,...
                                            'obs_cS',S_coex,'exp_cS',output.exp_cS);
                                tieline_output = struct('aC',output.aC,'bC',output.bC,'fpafpb',output.fpafpb,'fafb',output.fafb,'Kp',output.Kp);
                                disp(sprintf('done!'));
                                return % fitting finished
                            case 'grid'
                                if length(cp) == 2
                                    cp = linspace(start,finish,50*(finish-start));
                                end
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
                                        disp(sprintf('starting search and fitting...'));
                                        options = optimset('MaxFunEvals',500,'TolFun',1e-8,'TolX',1e-8,'MaxIter',500,'LargeScale','on','LevenbergMarquardt','on','Display','off','Jacobian','off');
                                        data = reshape(Scoex,ncB*ncS,1);
                                        npositions = length(cp);
                                        for i = 1:npositions
                                            x0 = [cp(i) 1 kp(1)+(kp(end)-kp(1)).*0.5]';% kp(1)+(kp(end)-kp(1)).*rand(1)
                                            lb = [cp(i)-0.001 1e-10 kp(1)]';
                                            ub = [cp(i)+0.001 Inf kp(end)]';
                                            [x,resnorm] = lsqcurvefit(@ruled_tieline_fit_fxn,x0,Ccoex,data,lb,ub,options,Scoex,Cbdy,Sbdy,boundary,tieline,part_coeff,spectra);
                                            [expcS,output] = ruled_tieline_fit_fxn(x,Ccoex,Scoex,Cbdy,Sbdy,boundary,tieline,part_coeff,spectra);
                                            % output =
                                            % struct('aC',alphaC,'bC',betaC,'aS',alphaS,'bS',betaS,
                                            % 'exp_cS',expcS,'fpafpb',fpafpb,'fafb',faf
                                            % b,'Kp',Kp,'ChiSquare',ChiSquare,'condn',condn);
                                            Chi2(i) = output.ChiSquare;
                                            % pack up variables for output
                                            fit_output(i) = struct('Chisquare',output.ChiSquare,'search_params',x','condn',output.condn);
                                            spectra_output(i) = struct('magnetic_field',spectra.magnetic_fields,'aS',output.aS,'bS',output.bS,...
                                                'obs_cS',Scoex,'exp_cS',output.exp_cS);
                                            tieline_output(i) = struct('aC',output.aC,'bC',output.bC,'fpafpb',output.fpafpb,'fafb',output.fafb,'Kp',output.Kp);
                                            pack;
                                            disp(sprintf('%f%% done',(i/npositions)*100));
                                        end
                                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                                        % analyze results of search and fitting   
                                        i = find(min(Chi2) == Chi2);
                                        if isvector(i) 
                                            for ii = 1:length(i)
                                                best_chisquare(ii,1) = Chi2(i(ii));
                                                best_indices(ii) = i(ii);
                                                best_params(ii,:) = fit_output(i(ii)).search_params;
                                                % plot boundary and alpha and beta phase compositions
                                                % and plot tielines.
                                                plot_tielines(tieline_output(i(ii)).aC,tieline_output(i(ii)).bC,Cbdy);
                                            end 
                                        else
                                            best_chisquare = Chi2(i);
                                            best_indices = i;
                                            best_params = fit_output(i).search_params;
                                            % plot boundary and alpha and beta phase compositions
                                            plot_tielines(tieline_output(i).aC,tieline_output(i).bC,Cbdy);
                                        end
                                        display = [best_chisquare';best_indices';best_params'];
                                        disp(sprintf('best chisquare = %f, with index %d and parameters [%f %f %f]\n',display));
                %                         disp(sprintf('best slope = %f => angle = %f\n',[best_params;slope2angle(best_params)]));
                                        disp(sprintf('done!'));
                                        return % fitting finished
                                    case 'function'
                                        disp(sprintf('starting search and fitting...'));
                                        options = optimset('MaxFunEvals',500,'TolFun',1e-8,'TolX',1e-8,'MaxIter',500,'LargeScale','on','LevenbergMarquardt','on','Display','off','Jacobian','off');
                                        data = reshape(Scoex,ncB*ncS,1);
                                        npositions = length(cp);
                                        for i = 1:npositions
                                            x0 = [c1 c2(i) 1 0 0]';
                                            lb = [c1-0.001 c2(i)-0.001 1e-10 -Inf -Inf]';
                                            ub = [c1+0.001 c2(i)+0.001 Inf Inf Inf]';
                                            [x,resnorm] = lsqcurvefit(@ruled_tieline_fit_fxn,x0,Ccoex,data,lb,ub,options,Scoex,Cbdy,Sbdy,boundary,tieline,part_coeff,spectra);
                                            [expcS,output] = ruled_tieline_fit_fxn(x,Ccoex,Scoex,Cbdy,Sbdy,boundary,tieline,part_coeff,spectra);
                                            % output =
                                            % struct('aC',alphaC,'bC',betaC,'aS',alphaS,'bS',betaS,
                                            % 'exp_cS',expcS,'fpafpb',fpafpb,'fafb',faf
                                            % b,'Kp',Kp,'ChiSquare',ChiSquare,'condn',condn);
                                            Chi2(i) = output.ChiSquare;
                                            % pack up variables for output
                                            fit_output(i) = struct('Chisquare',output.ChiSquare,'search_params',x','condn',output.condn);
                                            spectra_output(i) = struct('magnetic_field',spectra.magnetic_fields,'aS',output.aS,'bS',output.bS,...
                                                'obs_cS',Scoex,'exp_cS',output.exp_cS);
                                            tieline_output(i) = struct('aC',output.aC,'bC',output.bC,'fpafpb',output.fpafpb,'fafb',output.fafb,'Kp',output.Kp);
                                            pack;
                                            disp(sprintf('%f%% done',(i/npositions)*100));
                                        end
                                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                                        % analyze results of search and fitting   
                                        i = find(min(Chi2) == Chi2);
                                        if isvector(i) 
                                            for ii = 1:length(i)
                                                best_chisquare(ii,1) = Chi2(i(ii));
                                                best_indices(ii) = i(ii);
                                                best_params(ii,:) = fit_output(i(ii)).search_params;
                                                % plot boundary and alpha and beta phase compositions
                                                % and plot tielines.
                                                plot_tielines(tieline_output(i(ii),k(ii)).aC,tieline_output(i(ii),k(ii)).bC,Cbdy);
                                            end 
                                        else
                                            best_chisquare = Chi2(i);
                                            best_indices = i;
                                            best_params = fit_output(i).search_params;
                                            % plot boundary and alpha and beta phase compositions
                                            plot_tielines(tieline_output(i,k).aC,tieline_output(i,k).bC,Cbdy);
                                        end
                                        display = [best_chisquare';best_indices';best_params'];
                                        disp(sprintf('best chisquare = %f, with index %d and parameters [%f %f %f %f]\n',display));
                %                         disp(sprintf('best slope = %f => angle = %f\n',[best_params;slope2angle(best_params)]));
                                        disp(sprintf('done!'));
                                        return % fitting finished  
                                    otherwise
                                        error('invalid Kp configuration');
                                end
                            otherwise
                                error('invalid boundary search method');
                        end
                    otherwise
                        error('invalid critical point 2 position on boundary');
                end
            otherwise
                error('invalid tieline configuration');
        end
    case '2endtls'
        switch tieline.config
            case 'parallel'
                % if end tielines parallel, only one configuration
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