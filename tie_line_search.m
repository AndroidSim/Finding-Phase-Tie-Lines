function [fit_output,tl_output] = tie_line_search(S_traj,C_traj,S_bdy,C_bdy,varargin)
% tie_line_search tries to find a set of thermodynamic tie lines for the
% given phase boundary by fitting the ESR spectra within the phase boundary
% as a linear combination of the ESR spectra on the phase boundary.
% 
% input arguments:
% S_traj = ESR spectra within the phase boundary or within the phase
% coexistence region.
% C_traj = the compositions of the spectra in S_traj.
% S_bdy = ESR spectra on the phase boundary.
% C_bdy = the composition of the spectra in S_bdy.
%
% variable input arguments:
%   varargin{1} = the search method as one of the following:
%       'continuous' =  search parameters vary continously.
%       'grid' = grid search.
%   varargin{2} = only applicable for the grid searches.
%   specifies the increment of the grid search.  must be
%   between 0 and 1.  leave empty if continuous variation.

disp(sprintf('checking input arguments:'));

% checking trajectory spectra
if all(size(S_traj) > 1) & ndims(S_traj) == 2 % if spectra is a matrix
    [ntb,ntcol] = size(S_traj);

    if ntcol < 2 || rem(ntcol,2) ~= 0
        error('each trajectory spectrum is two columns: [B-field intensity_values]');
    end
    
    nts = ntcol/2;
else
    error('first argument must be a matrix of the trajectory spectra');
end

% checking trajectory compositions
if  any(size(C_traj) == 1) & any(size(C_traj) > 1) % if compositions is a vector
    [nrow,ncol] = size(C_traj);
    
    if ncol == 1 % if column vector
        C_traj = C_traj'; % convert to row vector
    end
    
    ntc = 1; 
    ntd = length(C_traj);
    
    if ntd == 2
        t_tern = false;
    elseif ntd == 3
        t_tern = true;
        C_traj = tern2cart(C_traj,1);
    end
    
    if ntc ~= nts
        error('number of trajectory compositions must equal number of trajectory spectra');
    end
elseif all(size(C_traj) > 1) & ndims(C_traj) == 2 % if compositions is a matrix 
    [ntc,ntd] = size(C_traj);
    
    if ntd == 2
        t_tern = false;
    elseif ntd == 3
        t_tern = true;
        C_traj = tern2cart(C_traj,1);
    end
    
    if ntc ~= nts
        error('number of trajectory compositions must equal number of trajectory spectra');
    end
end

% checking boundary spectra
if all(size(S_bdy) > 1) & ndims(S_bdy) == 2 % if spectra is a matrix
    [nbb,nbcol] = size(S_bdy);

    if nbcol < 2 || rem(nbcol,2) ~= 0
        error('each boundary spectrum is two columns: [B-field intensity_values]');
    end
    
    nbs = nbcol/2;
else
    error('third argument must be a matrix of the boundary spectra');
end

% checking boundary compositions
if  any(size(C_bdy) == 1) & any(size(C_bdy) > 1) % if compositions is a vector
    [nrow,ncol] = size(C_bdy);
    
    if nrow == 1 % if row vector
        nbc = 1; 
        nbd = length(C_bdy);   
    else
        nbc = length(C_bdy); 
        nbd = 1;
        C_bdy = C_bdy'; % convert to row vector
    end
    
    if nbc ~= nbs
        error('number of boundary compositions must equal number of boundary spectra');
    end
elseif all(size(C_bdy) > 1) & ndims(C_bdy) == 2 % if compositions is a matrix 
    [nbc,nbd] = size(C_bdy);
    
    if nbd == 2
        b_tern = false;
    elseif nbd == 3
        b_tern = true;
        C_bdy = tern2cart(C_bdy,1);
    end
    
    if nbc ~= nbs
        error('number of boundary compositions must equal number of boundary spectra');
    end
end

% process varargin to get search method, grid search parameter (if
% neeeded), and fit method
if isempty(varargin)
    search_method = 'continuous'; % search method
    slope_search_param = []; % if grid search, increment for grid search.
    Kp_search_param = [];
else
    if isempty(varargin{1})
        search_method = 'continuous';
    else
        if ischar(varargin{1})
            if any(strcmp(varargin{1},{'continuous';'grid';'mc';'mcmc'}))
                search_method = varargin{1};
            else
                error('search method must be either "continuous", "grid", "mc", or "mcmc"');
            end
        else
            error('search method must be either "continuous", "grid", "mc", or "mcmc"');
        end
    end
    
    if isempty(varargin{2})
        if strcmp(varargin{1},'grid')
            error('if grid search, second variable argument must be a number between 0 and 1');
        end
    else
        if isnumeric(varargin{2}) & ~isempty(varargin{2})
            if varargin{2} <= 0 | varargin{2} >= 1;
                error('slope (basis) search parameter must be a number between 0 and 1');
            else
                slope_search_param = varargin{2};
            end
        else
            error('slope (basis) search parameter must be a number between 0 and 1');
        end
    end
    
    if isempty(varargin{3})
        if strcmp(varargin{1},'grid')
            error('if grid search, third variable argument must be the Kp search increment');
        end
    else
        if isnumeric(varargin{3}) & ~isempty(varargin{3})
            if varargin{3} <= 0 | varargin{3} >= 1;
                error('Kp search increment parameter must be a number between 0 and 1');
            else
                Kp_search_param = varargin{3};
            end
        else
            error('Kp search increment parameter must be a number between 0 and 1');
        end
    end
end

disp(sprintf('processing boundary:'));
% processing boundary (potential area of improvement)
[bdy_len,bdy_p] = bdy_fxn(C_bdy,'linear');
bdy_p = bdy_p(:,1); % boundary parameter, goes from 0 to 1

disp(sprintf('processing spectra:'));
% processing spectra (potential area of improvement)
alinedS = align_spectra([S_traj S_bdy]); % align all spectra
alinedS = alinedS([300:end-300],:); % cutoff baseline (needs more work here)

% normalize boundary spectra
% alinedS(:,2:2:end) = normalize_vector(alinedS(:,2:2:end));

% normalize by converting to absorbance spectrum, transform so that the
% area of the absorbance spectrum equals one, then convert back to the
% derivative spectrum.
% alinedS = deriv2abs(alinedS);
% alinedS = normalize_spectra(alinedS,'a');
% alinedS = abs2deriv(alinedS);

S_traj = alinedS(:,1:ntcol);
S_bdy = alinedS(:,ntcol+1:ntcol+nbcol);
clear alinedS;
[ntb,nts] = size(S_traj);
[ntc,ntd] = size(C_traj);
[nbb,nbs] = size(S_bdy);
[nbc,nbd] = size(C_bdy);

% smooth composition intensities of boundary spectra
for b = 1:nbb
    S_bdy(b,2:2:end) = smooth(S_bdy(b,2:2:end),3,'moving')';
end

% starting search:
disp(sprintf('starting search:'));

switch search_method
    case 'continuous'
        options = optimset('display','iter');
        ISt = S_traj(:,2:2:end);
        ISt = ISt./repmat(sum(ISt),ntb,1);
        vISt = reshape(ISt',ntb*ntc,1);
    case 'grid'
        slope_angles = (pi/2).*[0:slope_search_param:1];
        slopes = angle2slope(slope_angles);
        nslopes = length(slopes);
        Kp = [0.5:Kp_search_param:1 1+round(Kp_search_param*10):round(Kp_search_param*10):2];
        nKp = length(Kp);
        
        % loop over data pts within phase boundary (trajectory
        % points)
        for t = 1:ntc
            St = S_traj(:,t*2-1:t*2);
            Ct = C_traj(t,:);
            
            for m = 1:nslopes
                for k = 1:nKp
                    fit_output = tie_line_m_kp_fit_fxn(slopes(m),Kp(k),Ct,St,C_bdy,S_bdy,search_method);
                    % fit_output = struct('Ca',Ca,'Cb',Cb,'Sa',Sa,'Sb',Sb,'fafb',fafb,'chisq',chisq);
                    chisq{t}(m,k) = fit_output.chisq;
                    Ca{m,k,t} = fit_output.Ca;
                    Cb{m,k,t} = fit_output.Cb;
                    Sa{m,k,t} = fit_output.Sa;
                    Sb{m,k,t} = fit_output.Sb;
                    fafb{m,k,t} = fit_output.fafb;
                end
            end
            
            minchisq{t} = min(min(chisq{t}));
            [m,k] = find(minchisq{t} == chisq{t});
            x{t} = [[m k];[slope2angle(slopes(m))'*(180/pi) Kp(k)']];
            disp(sprintf('finished fitting %d out of %d trajectory points',t,ntc));
        end
        
        fit_output = struct('best_chisq',minchisq,'best_fit_parameters',x,'chisq_surface',chisq);
        tl_output = struct('Ca',Ca,'Cb',Cb,'Sa',Sa,'Sb',Sb,'fafb',fafb);
    
        % plot boundary and alpha and beta phase compositions
        figure;
        ternary_plot(cart2tern(C_bdy,1),'-k','linewidth',3);
        for t = 1:ntc
            alpha_pt(t,:) = Ca{x{t}(1,1),x{t}(1,2),t};
            beta_pt(t,:) = Cb{x{t}(1,1),x{t}(1,2),t};
            ternary_plot(cart2tern(alpha_pt(t,:),1),'.r','markersize',20);
            ternary_plot(cart2tern(beta_pt(t,:),1),'.b','markersize',20);
            % plot tie line
            ternary_plot(cart2tern([alpha_pt(t,:);beta_pt(t,:)],1),'-m','linewidth',2);
        end    
    case 'mc'
%         slope_angles = (pi/2).*[0:slope_search_param:1];
%         slopes = angle2slope(slope_angles);
%         nslopes = length(slopes);
        Kp = [0.5:Kp_search_param:1 1+round(Kp_search_param*10):round(Kp_search_param*10):2];
        nKp = length(Kp);
        niters = 3000;
        
        % start off with initial tie line configuration (parallel tie
        % lines)
        slopes = zeros(ntc,1);
%         slopes = repmat(angle2slope(pi/2),ntc,1);
        for t = 1:ntc
            St = S_traj(:,t*2-1:t*2);
            Ct = C_traj(t,:);
            
            for k = 1:nKp
                output = tie_line_m_kp_fit_fxn(slopes(t),Kp(k),Ct,St,C_bdy,S_bdy,search_method);
                % fit_output = struct('Ca',Ca,'Cb',Cb,'Sa',Sa,'Sb',Sb,'fafb',fafb,'chisq',chisq);
                chisq(t,k) = output.chisq;
                C_a(t,:) = output.Ca;
                C_b(t,:) = output.Cb;
                Sa(:,t) = output.Sa;
                Sb(:,t) = output.Sb;
                fafb(t,:) = output.fafb;
            end
        end
        
        [M,I] = min(chisq,[],2);
        minchisq0 = mean(M);
        Kps = Kp(I);
        tl_output = struct('Ca',C_a,'Cb',C_b,'Sa',Sa,'Sb',Sb,'fafb',fafb);
        
        n = 1;
        while (n <= niters) 
            % pick random slope for 1st data point
            St = S_traj(:,1:2);
            Ct = C_traj(1,:);
            at = pi*rand(1);
            st = angle2slope(at);
            Cab = tie_line_m_kp_fit_fxn(st,0,Ct,St,C_bdy,S_bdy,search_method); % Cab = [C_alpha;C_beta];
            Ca = Cab(1,:);
            Cb = Cab(2,:);
            Calpha(1,:) = Ca;
            Cbeta(1,:) = Cb;
            slopes(1) = st;
            
            % loop over data points and randomly choose a tie line that does
            % not intersect the tie line through the previous data point
            for t = 2:ntc
                St = S_traj(:,t*2-1:t*2);
                Ct = C_traj(t,:);
	%             aab = acos(dot(Ca-Ct,Cb-Ct)/(norm(Ca-Ct)*norm(Cb-Ct)));
                aa = pt2angle(Ca-Ct);
                ab = pt2angle(Cb-Ct);
                % 8 cases:
                done = false;
                if aa >= 0 & ab < aa-pi
                    while ~done
                        at = pi*rand(1);
                        if (at < aa & at > ab+pi) | (at > ab & at < aa-pi) 
                            done = true;
                        end
                    end
                elseif aa >= 0 & ab > aa-pi
                    while ~done
                        at = pi*rand(1);
                        if (at > aa & at < ab+pi) | (at < ab & at > aa-pi) 
                            done = true;
                        end
                    end
                elseif ab >= 0 &  aa < ab-pi
                    while ~done
                        at = pi*rand(1);
                        if (at > aa & at < ab-pi) | (at < ab & at > aa+pi) 
                            done = true;
                        end
                    end
                elseif ab >= 0 &  aa > ab-pi 
                    while ~done
                        at = pi*rand(1);
                        if (at < aa & at > ab-pi) | (at > ab & at < aa+pi) 
                            done = true;
                        end
                    end
                elseif (aa <= 0 & ab <= 0) & pt2lin_dist(Ca,[Cb;Ct],2) > 0 
                    while ~done
                        at = pi*rand(1);
                        if (at < aa & at > ab+pi) | (at > ab & at < aa+pi) 
                            done = true;
                        end
                    end
                elseif (aa <= 0 & ab <= 0) & pt2lin_dist(Ca,[Cb;Ct],2) < 0 
                    while ~done
                        at = pi*rand(1);
                        if (at > aa & at < ab+pi) | (at < ab & at > aa+pi) 
                            done = true;
                        end
                    end
                elseif (aa >= 0 & ab >= 0) & pt2lin_dist(Ca,[Cb;Ct],2) > 0 
                    while ~done
                        at = pi*rand(1);
                        if (at < aa & at > ab-pi) | (at > ab & at < aa-pi)
                            done = true;
                        end
                    end
                elseif (aa >= 0 & ab >= 0) & pt2lin_dist(Ca,[Cb;Ct],2) < 0 
                    while ~done
                        at = pi*rand(1);
                        if (at > aa & at < ab-pi) | (at < ab & at > aa-pi)
                            done = true;
                        end
                    end
                end
                st = angle2slope(at);
                Cab = tie_line_m_kp_fit_fxn(st,0,Ct,St,C_bdy,S_bdy,search_method); % Cab = [C_alpha;C_beta];
                Ca = Cab(1,:);
                Cb = Cab(2,:);
                Calpha(t,:) = Ca;
                Cbeta(t,:) = Cb;
                slopes(t) = st;
            end
            
            % debugging
            % plot boundary and alpha and beta phase compositions
%             figure;
%             ternary_plot(cart2tern(C_bdy,1),'-k','linewidth',3);
%             for t = 1:ntc
%                 ternary_plot(cart2tern(Calpha(t,:),1),'.r','markersize',20);
%                 ternary_plot(cart2tern(Cbeta(t,:),1),'.b','markersize',20);
%                 % plot tie line
%                 ternary_plot(cart2tern([Calpha(t,:);Cbeta(t,:)],1),'-m','linewidth',2);
%             end    
            
            % loop over data points again where now each data point has a
            % randomly chosen slope but where the tie lines do not
            % intersect
            for t = 1:ntc
                St = S_traj(:,t*2-1:t*2);
                Ct = C_traj(t,:);
                
                for k = 1:nKp
                    output = tie_line_m_kp_fit_fxn(slopes(t),Kp(k),Ct,St,C_bdy,S_bdy,search_method);
                    % fit_output = struct('Ca',Ca,'Cb',Cb,'Sa',Sa,'Sb',Sb,'fafb',fafb,'chisq',chisq);
                    chisq(t,k) = output.chisq;
                    C_a(t,:) = output.Ca;
                    C_b(t,:) = output.Cb;
                    Sa(:,t) = output.Sa;
                    Sb(:,t) = output.Sb;
                    fafb(t,:) = output.fafb;
                end
            end
            
            [M,I] = min(chisq,[],2);
            minchisq = mean(M);
            
            % keep tie line set if minchisq is lower than previous iteration
%             if minchisq < minchisq0
%                 minchisq0 = minchisq;
%                 Kps = Kp(I);
%                 fit_output = struct('best_chisq',minchisq,'Kps',Kps,'chisq_surface',chisq);
%                 tl_output = struct('Ca',C_a,'Cb',C_b,'Sa',Sa,'Sb',Sb,'fafb',fafb);
%             end
            
            if n > 5
                [best_chisqs,index] = sort(best_chisqs);
                fit_output = fit_output(index);
                tl_output = tl_output(index);
                i = find(minchisq < best_chisqs);
                
                if ~isempty(i)
                    i = max(i);
                    best_chisqs(i) = minchisq;
                    Kps = Kp(I);
                    fit_output(i) = struct('min_chisq',minchisq,'Kps',Kps,'chisq_surface',chisq);
                    tl_output(i) = struct('Ca',C_a,'Cb',C_b,'Sa',Sa,'Sb',Sb,'fafb',fafb);
                end
            else
                best_chisqs(n) = minchisq;
                Kps = Kp(I);
                fit_output(n) = struct('min_chisq',minchisq,'Kps',Kps,'chisq_surface',chisq);
                tl_output(n) = struct('Ca',C_a,'Cb',C_b,'Sa',Sa,'Sb',Sb,'fafb',fafb);
            end
            
            n = n+1;
        end
        
        % plot tie lines
        for n = 1:5
            figure;
            ternary_plot(cart2tern(C_bdy,1),'-k','linewidth',3);
            alpha_pts = tl_output(n).Ca;
            beta_pts = tl_output(n).Cb;
			for t = 1:ntc
			    ternary_plot(cart2tern(alpha_pts(t,:),1),'.r','markersize',20);
				ternary_plot(cart2tern(beta_pts(t,:),1),'.b','markersize',20);
				ternary_plot(cart2tern([alpha_pts(t,:);beta_pts(t,:)],1),'-m','linewidth',2);
			end  
			ternary_plot(cart2tern(C_traj,1),'.g','markersize',20);
        end
    case 'mcmc'
        % each tie line configuration is a structure with the following fields:
        % a = slope angle vector
        % Kp = best fit Kp vector
        % pt = data composition point matrix
        % apt = alpha (Ld) phase composition point matrix
        % bpt = beta (Lo) phase composition point matrix
        tie_lines =   struct('a',{},'m',{},'Kp',{},'C',{},'Ca',{},'Cb',{},'S',{},'Sa',{},'Sb',{},'chisq',{});
        
        % start off with initial tie line configuration (parallel tie
        % lines)
%         slopes = repmat(angle2slope(pi/2),ntc,1);
%         slopes = zeros(ntc,1);
        for t = 1:ntc
            tie_lines(t).m = 0;
            tie_lines(t).a = slope2angle(0);
            tie_lines(t).Kp = 1;
            tie_lines(t).C = C_traj(t,:);
            tie_lines(t).S = S_traj(:,t*2-1:t*2);
            output = tie_line_m_kp_fit_fxn(tie_lines(t).m,tie_lines(t).Kp,tie_lines(t).C,tie_lines(t).S,C_bdy,S_bdy,search_method);
            % output = struct('Ca',Ca,'Cb',Cb,'Sa',Sa,'Sb',Sb,'fafb',fafb,'chisq',chisq);
            tie_lines(t).Ca = output.Ca;
            tie_lines(t).Cb = output.Cb;
            tie_lines(t).Sa = output.Sa;
            tie_lines(t).Sb = output.Sb;
            tie_lines(t).chisq = output.chisq;
        end
        
        chisq_v = cat(1,tie_lines.chisq);
        minchisq0 = min(chisq_v); 
        niters = 3000;
        n = 1;
        while (n <= niters) 
            % loop over data points(and thus tie lines)
            % at each point, randomly choose a slope within a
            % certain range contrained by neighboring tie lines and
            % randomly choose a Kp.
            % with slope and Kp, calculate chisquare, and with previous
            % chisquare at this point use a metropolis algorithm to
            % determine whether to keep the new tie line or reject it and
            % keep the old tie line.
            St = S_traj(:,1:2);
            Ct = C_traj(1,:);
            at = (pi/2)*rand(1);
            st = angle2slope(at);
            Cab = tie_line_m_kp_fit_fxn(st,0,Ct,St,C_bdy,S_bdy,search_method); % Cab = [C_alpha;C_beta];
            Ca = Cab(1,:);
            Cb = Cab(2,:);
            Calpha(1,:) = Ca;
            Cbeta(1,:) = Cb;
            slopes(1) = st;
            
            % loop over data points and randomly choose a tie line that does
            % not intersect the tie line through the previous data point
            for t = 2:ntc
                St = S_traj(:,t*2-1:t*2);
                Ct = C_traj(t,:);
	%             aab = acos(dot(Ca-Ct,Cb-Ct)/(norm(Ca-Ct)*norm(Cb-Ct)));
                aa = pt2angle(Ca-Ct);
                ab = pt2angle(Cb-Ct);
                % 8 cases:
                done = false;
                if aa >= 0 & ab < aa-pi
                    while ~done
                        at = (pi/2)*rand(1);
                        if (at < aa & at > ab+pi) | (at > ab & at < aa-pi) 
                            done = true;
                        end
                    end
                elseif aa >= 0 & ab > aa-pi
                    while ~done
                        at = (pi/2)*rand(1);
                        if (at > aa & at < ab+pi) | (at < ab & at > aa-pi) 
                            done = true;
                        end
                    end
                elseif ab >= 0 &  aa < ab-pi
                    while ~done
                        at = (pi/2)*rand(1);
                        if (at > aa & at < ab-pi) | (at < ab & at > aa+pi) 
                            done = true;
                        end
                    end
                elseif ab >= 0 &  aa > ab-pi 
                    while ~done
                        at = (pi/2)*rand(1);
                        if (at < aa & at > ab-pi) | (at > ab & at < aa+pi) 
                            done = true;
                        end
                    end
                elseif (aa <= 0 & ab <= 0) & pt2lin_dist(Ca,[Cb;Ct],2) > 0 
                    while ~done
                        at = (pi/2)*rand(1);
                        if (at < aa & at > ab+pi) | (at > ab & at < aa+pi) 
                            done = true;
                        end
                    end
                elseif (aa <= 0 & ab <= 0) & pt2lin_dist(Ca,[Cb;Ct],2) < 0 
                    while ~done
                        at = (pi/2)*rand(1);
                        if (at > aa & at < ab+pi) | (at < ab & at > aa+pi) 
                            done = true;
                        end
                    end
                elseif (aa >= 0 & ab >= 0) & pt2lin_dist(Ca,[Cb;Ct],2) > 0 
                    while ~done
                        at = (pi/2)*rand(1);
                        if (at < aa & at > ab-pi) | (at > ab & at < aa-pi)
                            done = true;
                        end
                    end
                elseif (aa >= 0 & ab >= 0) & pt2lin_dist(Ca,[Cb;Ct],2) < 0 
                    while ~done
                        at = (pi/2)*rand(1);
                        if (at > aa & at < ab-pi) | (at < ab & at > aa-pi)
                            done = true;
                        end
                    end
                end
                st = angle2slope(at);
                Cab = tie_line_m_kp_fit_fxn(st,0,Ct,St,C_bdy,S_bdy,search_method); % Cab = [C_alpha;C_beta];
                Ca = Cab(1,:);
                Cb = Cab(2,:);
                Calpha(t,:) = Ca;
                Cbeta(t,:) = Cb;
                slopes(t) = st;
            end
            
            % debugging
            % plot boundary and alpha and beta phase compositions
%             figure;
%             ternary_plot(cart2tern(C_bdy,1),'-k','linewidth',3);
%             for t = 1:ntc
%                 ternary_plot(cart2tern(Calpha(t,:),1),'.r','markersize',20);
%                 ternary_plot(cart2tern(Cbeta(t,:),1),'.b','markersize',20);
%                 % plot tie line
%                 ternary_plot(cart2tern([Calpha(t,:);Cbeta(t,:)],1),'-m','linewidth',2);
%             end    
            
            % loop over data points again where now each data point has a
            % randomly chosen slope but where the tie lines do not
            % intersect
            for t = 1:ntc
                St = S_traj(:,t*2-1:t*2);
                Ct = C_traj(t,:);
                
                for k = 1:nKp
                    output = tie_line_m_kp_fit_fxn(slopes(t),Kp(k),Ct,St,C_bdy,S_bdy,search_method);
                    % fit_output = struct('Ca',Ca,'Cb',Cb,'Sa',Sa,'Sb',Sb,'fafb',fafb,'chisq',chisq);
                    chisq(t,k) = output.chisq;
                    C_a(t,:) = output.Ca;
                    C_b(t,:) = output.Cb;
                    Sa(:,t) = output.Sa;
                    Sb(:,t) = output.Sb;
                    fafb(t,:) = output.fafb;
                end
            end
            
            [M,I] = min(chisq,[],2);
            minchisq = mean(M);
            
            % keep tie line set if minchisq is lower than previous iteration
%             if minchisq < minchisq0
%                 minchisq0 = minchisq;
%                 Kps = Kp(I);
%                 fit_output = struct('best_chisq',minchisq,'Kps',Kps,'chisq_surface',chisq);
%                 tl_output = struct('Ca',C_a,'Cb',C_b,'Sa',Sa,'Sb',Sb,'fafb',fafb);
%             end
            
            if n > 5
                [best_chisqs,index] = sort(best_chisqs);
                fit_output = fit_output(index);
                tl_output = tl_output(index);
                i = find(minchisq < best_chisqs);
                
                if ~isempty(i)
                    i = max(i);
                    best_chisqs(i) = minchisq;
                    Kps = Kp(I);
                    fit_output(i) = struct('min_chisq',minchisq,'Kps',Kps,'chisq_surface',chisq);
                    tl_output(i) = struct('Ca',C_a,'Cb',C_b,'Sa',Sa,'Sb',Sb,'fafb',fafb);
                end
            else
                best_chisqs(n) = minchisq;
                Kps = Kp(I);
                fit_output(n) = struct('min_chisq',minchisq,'Kps',Kps,'chisq_surface',chisq);
                tl_output(n) = struct('Ca',C_a,'Cb',C_b,'Sa',Sa,'Sb',Sb,'fafb',fafb);
            end
            
            n = n+1;
        end
        
        % plot tie lines
        for n = 1:5
            figure;
            ternary_plot(cart2tern(C_bdy,1),'-k','linewidth',3);
            alpha_pts = tl_output(n).Ca;
            beta_pts = tl_output(n).Cb;
			for t = 1:ntc
			    ternary_plot(cart2tern(alpha_pts(t,:),1),'.r','markersize',20);
				ternary_plot(cart2tern(beta_pts(t,:),1),'.b','markersize',20);
				ternary_plot(cart2tern([alpha_pts(t,:);beta_pts(t,:)],1),'-m','linewidth',2);
			end  
			ternary_plot(cart2tern(C_traj,1),'.g','markersize',20);
        end
    otherwise
        error('invalid search_method');
end

disp(sprintf('search done!'));

return