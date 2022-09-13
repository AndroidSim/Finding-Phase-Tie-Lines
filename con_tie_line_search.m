function fit_output = con_tie_line_search(S_traj,C_traj,S_bdy,C_bdy,varargin)

disp(sprintf('processing input arguments:'));

% checking trajectory spectra
if all(size(S_traj) > 1) & ndims(S_traj) == 2 % if S is a matrix
    [ntb,ntcol] = size(S_traj);

    if ntcol < 2 || rem(ntcol,2) ~= 0
        error('each trajectory spectrum is two columns: [B-field intensity_values]');
    end
    
    nts = ntcol/2;
else
    error('first argument must be a matrix of the trajectory spectra');
end

% checking trajectory compositions
if  any(size(C_traj) == 1) & any(size(C_traj) > 1) % if C is a vector
    [nrow,ncol] = size(C_traj);
    
    if nrow == 1 % if row vector
        ntc = 1; 
        ntd = length(C_traj);
        C_traj = C_traj'; % convert to column vector
    else
        ntc = length(C_traj); 
        ntd = 1;
    end
    
    if ntc ~= nts
        error('number of trajectory compositions must equal number of spectra');
    end
elseif all(size(C_traj) > 1) & ndims(C_traj) == 2 % if C is a matrix 
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
if all(size(S_bdy) > 1) & ndims(S_bdy) == 2 % if S is a matrix
    [nbb,nbcol] = size(S_bdy);

    if nbcol < 2 || rem(nbcol,2) ~= 0
        error('each boundary spectrum is two columns: [B-field intensity_values]');
    end
    
    nbs = nbcol/2;
else
    error('third argument must be a matrix of the boundary spectra');
end

% checking boundary compositions
if  any(size(C_bdy) == 1) & any(size(C_bdy) > 1) % if C is a vector
    [nrow,ncol] = size(C_bdy);
    
    if nrow == 1 % if row vector
        nbc = 1; 
        nbd = length(C_bdy);
        C_bdy = C_bdy'; % convert to column vector
    else
        nbc = length(C_bdy); 
        nbd = 1;
    end
    
    if nbc ~= nbs
        error('number of boundary compositions must equal number of boundary spectra');
    end
elseif all(size(C_bdy) > 1) & ndims(C_bdy) == 2 % if C is a matrix 
    [nbc,nbd] = size(C_bdy);
    
    if nbd == 2
        b_tern = false;
    elseif nbd == 3
        b_tern = true;
        C_bdy = tern2cart(C_bdy,1);
    end
    
    if nbc ~= nbs
        error('number of boundary compositions must equal number of spectra');
    end
end

% processing boundary
[bdy_len,bdy_p] = bdy_fxn(C_bdy,'linear');
bdy_p = bdy_p(:,1);

% processing spectra
alinedS = align_spectra([S_traj S_bdy]);
alinedS([1:200 end-200:end],:) = []; 
S_traj = alinedS(:,1:ntcol);
S_bdy = alinedS(:,ntcol+1:end);
clear alinedS;

% process varargin to get search parameters
if isempty(varargin)
    crit_search_i = 0.01; % increment for Basis search along boundary (100 points)
    B_eval_method = 'lls'; % default basis evaluation method, 'lls' = linear least squares
else
    if isempty(varargin{1})
        crit_search_i = 0.01;
    else
        if isnumeric(varargin{1})
            if varargin{1} <= 0 & varargin{1} >= 1;
                error('critical point search parameter must be a number between 0 and 1');
            else
                crit_search_i = varargin{1};
            end
        else
            error('critical point search parameter must be a number between 0 and 1');
        end
    end
    
    if isempty(varargin{2})
        B_eval_method = 'lls';
    else
        if ischar(varargin{2})
            if any(strcmp(varargin{2},{'lls';'spectral';'svd'}))
                B_eval_method = varargin{2};
            else
                error('Basis evaluation method must be "lls" or "spectral" or "svd"');
            end
        else
            error('Basis evaluation method must be "lls" or "spectral" or "svd"');
        end
    end
end

% starting search:
% search parameters are two critical points.
% first, do grid search of all possible critical point pairs
ops = optimset('display','off','largescale','off');
options = optimset('display','off','jacobian','on');
crit_p = [0:crit_search_i:1];
ncritp = length(crit_p);

disp(sprintf('starting search:'));
% loop over first critical point
for c1 = 1:ncritp
    cpt1 = interp1(bdy_p,C_bdy,crit_p(c1),'linear');
    % cpt1 = interp_boundary(C_bdy,crit_p(c1),'linear')
    
    for c2 = c1+1:ncritp-1
        cpt2 = interp1(bdy_p,C_bdy,crit_p(c2),'linear');
        % cpt2 = interp_boundary(C_bdy,crit_p(c2),'linear')
        
        if isequal(cpt1,cpt2)
            continue;
        end
        
        % disp(sprintf('critical point 1 parameter = %f, critical point 2 parameter = %f',crit_p(c1),crit_p(c2)));
        
        % calculate the set of tie lines determined from the critical
        % points that pass through the trajectory points.
        
        % first find slope of tangent line at critical points
        temp1 = [crit_p(c1)-0.02:0.02:crit_p(c1)+0.02]';
        if any(temp1 < 0)
            temp1(temp1 < 0) = temp1(temp1 < 0)+1;
        end
        if any(temp1 > 1)
            temp1(temp1 > 1) = temp1(temp1 > 1)-1;
        end
        ptacpt1 = interp1(bdy_p,C_bdy,temp1,'linear')';
        mc1 = ppval(fnder(spline(temp1,ptacpt1)),crit_p(c1));
        warning off;
        mc1 = mc1(2)/mc1(1);
        warning on;
        
        temp2 = [crit_p(c2)-0.02:0.02:crit_p(c2)+0.02]';
        if any(temp2 < 0)
            temp2(temp2 < 0) = temp2(temp2 < 0)+1;
        end
        if any(temp2 > 1)
            temp2(temp2 > 1) = temp2(temp2 > 1)-1;
        end
        ptacpt2 = interp1(bdy_p,C_bdy,temp2,'linear')';
        mc2 = ppval(fnder(spline(temp2,ptacpt2)),crit_p(c2));
        warning off;
        mc2 = mc2(2)/mc2(1);
        warning on;
        clear temp1 temp2;

        if isequal(mc1,Inf) % mc1 == Inf
            lina = [cpt1;cpt1(1) 1];
        elseif mc1 > 0
            lina = [cpt1;1 mc1*(1-cpt1(1))+cpt1(2)];
        elseif isequal(mc1,0) % mc1 == 0
            lina = [cpt1;1 cpt1(2)];
        elseif mc1 < 0 
            lina = [cpt1;0 mc1*(0-cpt1(1))+cpt1(2)];
        end
        
        if isequal(mc2,Inf) % mc2 == Inf
            linb = [cpt2;cpt2(1) 1];
        elseif mc2 > 0
            linb = [cpt2;1 mc2*(1-cpt2(1))+cpt2(2)];
        elseif isequal(mc2,0) % mc2 == 0
            linb = [cpt2;1 cpt2(2)];
        elseif mc2 < 0 
            linb = [cpt2;0 mc2*(0-cpt2(1))+cpt2(2)];
        end

        cpt_m_int = line_line_int(lina,linb);
        clear lina linb;
        %if ~isfinite(cpt_m_int)
        %    keyboard
        %end

        % given the slope of the tie lines at the critical points, find another point on the line
        % that passes throught the trajectory point that is not the
        % trajectory point
        if isequal(mc1,mc2) % mc1 == mc2
            if isequal(mc2,Inf) % mc2 == Inf;
                pt = [C_traj(:,1) ones(ntc,1)];
            elseif mc2 > 0
                pt = [ones(ntc,1) mc2.*(1-C_traj(:,1))+C_traj(:,2)];
            elseif isequal(mc2,0) % mc2 == 0
                pt = [ones(ntc,1) C_traj(:,2)];
            elseif mc2 < 0 
                pt = [zeros(ntc,1) mc2.*(0-C_traj(:,1))+C_traj(:,2)];
            end
        else
            pt = repmat(cpt_m_int,ntc,1);  
        end
        
        % loop over trajectory points finding fraction of phases
        for t = 1:ntc
            St = S_traj(:,t*2-1:t*2);
            Ct = C_traj(t,:);
            
            % special condition
            if isequal(round(Ct.*1000),round(cpt1.*1000)) | isequal(round(Ct.*1000),round(cpt2.*1000))
                falpha(t) = 0.5;
                fbeta(t) = 0.5;
                continue;
            end
    
            % the intersection of the line drawn through the intersection
            % point of the critical point tangent line and the
            % trajectory point with the boundary gives both the composition
            % point of the alpha phase (C_alpha, or Ld phase) and the
            % composition of the beta phase (C_beta, or Lo phase).  the
            % points returned from boundary_line_int are that
            % C_alpha=int_pts(1,:) and C_beta=int_pts(2,:)
            int_pts = boundary_line_int(C_bdy,[pt(t,:);Ct]);
            nipt = size(int_pts,1);
            
            if isempty(int_pts) % no intersection
                error('tie line does not intersect boundary');
            elseif isequal(nipt,1) % tangent intersection
                error('tie line intersects boundary at only 1 point');
            elseif isequal(nipt,2) % 2 intersections
                C_alpha = int_pts(1,:);
                C_beta = int_pts(2,:);
            elseif nipt > 2 
                error('tie line intersects boundary at > 2 points');
            end
           
            % interpolate spectra at alpha point
            S_alpha = interp_spectra(S_bdy,C_bdy,'curve',[],[],C_alpha,'linear');
            % interpolate spectra at beta point
            S_beta = interp_spectra(S_bdy,C_bdy,'curve',[],[],C_beta,'linear');
            
            % calculate fraction of alpha phase and fraction of beta phase
            fbeta(t) = norm(Ct-C_alpha)/norm(C_alpha-C_beta);
            falpha(t) = 1-fbeta(t);
            
            % solve linear least squares problem or spectral linear combination
            % problem to determine Kp
            % kp defined as into beta (Lo) phase:
            %   fpalpha = fa./(fa+kp.*fb);
            %   fpbeta = (kp.*fb)./(fa+kp.*fb);
            Basis = [S_alpha S_beta]; % define Basis
            
            % reduced Basis and trajectory point spectrum to just intensity values
            Basis = Basis(:,2:2:end);
            ISt = St(:,2);
                    
            % normalize spectra so that the sum of intensities equals one
            Basis(:,1) = Basis(:,1)./sum(Basis(:,1));
            Basis(:,2) = Basis(:,2)./sum(Basis(:,2));
            ISt(:,1) = ISt(:,1)./sum(ISt(:,1));
        
            switch B_eval_method  
                case 'lls'
                    % linear least squares with constraints
                    % solve Ax = b for fraction of probe in alpha phase and
                    % then solve Kp equation as a fxn of fraction of probe in
                    % alpha phase and fraction of alpha phase
                    Aeq = [1 1];
                    beq = 1;
                    Aneq = [-1 0;0 -1;1 0;0 1];
                    bneq = [0;0;1;1];
                    [x,resnorm] = lsqlin(Basis,ISt,Aneq,bneq,Aeq,beq,[],[],[],ops);% x = lsqlin(C,d,A,b,Aeq,beq,lb,ub,x0,options)
                    fpalpha(t) = x(1);
                    fpbeta(t)= x(2);
                    chisq(t) = sqrt(resnorm);
                    kp(t) = (falpha(t)*(1-fpalpha(t)))/(fpalpha(t)*fbeta(t));
                case 'spectral'
                    [x,resnorm] = lsqcurvefit(@lincomb_fit_fun,1.01,fbeta(t),ISt,0.1,10,options,'spectral','spectrum_fit',Basis); 
                    chisq(t) = sqrt(resnorm);
                    kp(t) = x;
                case 'svd'
                    % Aeq = [1 1];
                    % beq = 1;
                    % [U,V,X,C,S] = gsvd(Basis,Aeq,0);
                    % W = inv(X');
                    % x = ((U(:,2)'*ISt)/C(2,2))*W(:,2)+((V(:,1)'*beq)/S(1,1))*W(:,1);
  
                    % [U,W,V] = svd(Basis,0);
                    % W = diag(W);
                    % x(1) = sum(((U(:,1)'*ISt)/W(1))*V(:,1));
                    % x(2) = sum(((U(:,2)'*ISt)/W(2))*V(:,2));
                    % xt = sum(x);
                    % if (x(1) > 0 & x(2) > 0) | (x(1) < 0 & x(2) < 0)
                    %     x = x./xt;
                    % elseif x(1) < 0 & x(2) > 0 & xt > 0
                    %     x = [(x(1)+xt)/xt; (x(2)-xt)/xt];
                    % elseif x(1) < 0 & x(2) > 0 & xt < 0
                    %     x = [(x(1)-xt)/xt; (x(2)+xt)/xt];
                    % elseif x(1) > 0 & x(2) < 0 & xt > 0
                    %     x = [(x(1)-xt)/xt; (x(2)+xt)/xt];
                    % elseif x(1) > 0 & x(2) < 0 & xt < 0
                    %     x = [(x(1)+xt)/xt; (x(2)-xt)/xt];
                    % end
                    
                    % convert to derivative spectra to absorbance spectra
                    % Basis = deriv2abs(Basis);
                    % St = deriv2abs(St);
                    % Basis = normalize_spectra(Basis,'a');
                    % St = normalize_spectra(St,'a');
                    
                    % perform singular value decomposition to obtain eigenspectra basis
                    % disp(sprintf('performing singular value decomposition to obtain eigenspectra basis:\n'));
                    [U,W,V] = svd(Basis,0);
                    B = U(:,1:2); % B = 2 component eigenbasis

                    % solve for eigenspectra coefficient matrix M by linear least squares
                    % constraints
                    % Aeq = [sum(U(:,1)) sum(U(:,2))];
                    Aeq = sum(B);
                    beq1 = sum(Basis(:,1));
                    beq2 = sum(Basis(:,2));
                    %A = -B;
                    %b = zeros(size(U(:,1)));
                    temp1 = lsqlin(B,Basis(:,1),[],[],Aeq,beq1,[],[],[],ops);
                    temp2 = lsqlin(B,Basis(:,2),[],[],Aeq,beq2,[],[],[],ops);
                    M = [temp1 temp2];
                    Aeq = [1 1];
                    beq = 1;
                    x = lsqlin(M,[ISt'*B(:,1);ISt'*B(:,2)],[],[],Aeq,beq,[],[],[],ops);
                    % x = inv(M)*[ISt'*B(:,1);ISt'*B(:,2)];
                        
                    fpalpha(t) = x(1);
                    fpbeta(t)= x(2);
                    chisq(t) = norm((Basis*x)-ISt);
                    kp(t) = (falpha(t)*(1-fpalpha(t)))/(fpalpha(t)*fbeta(t));
                otherwise
                    error('invalid fit method');
            end
        
            Calpha(t,:) = C_alpha;
            Cbeta(t,:) = C_beta;
        end
    
        cpts{c1,c2} = [cpt1;cpt2]; 
        switch B_eval_method
            case 'lls'
                Ca{c1,c2} = Calpha;
                Cb{c1,c2} = Cbeta;
                fpa{c1,c2} = fpalpha';
                fpb{c1,c2} = fpbeta';
                fa{c1,c2} = falpha';
                fb{c1,c2} = fbeta';
                Kp{c1,c2} = kp';
                chisquare{c1,c2} = chisq';
            case 'spectral'
                Ca{c1,c2} = Calpha;
                Cb{c1,c2} = Cbeta;
                fpa{c1,c2} = (falpha')./(falpha'+(kp'.*fbeta'));
                fpb{c1,c2} = (kp'.*fbeta')./(falpha'+(kp'.*fbeta'));
                fa{c1,c2} = falpha';
                fb{c1,c2} = fbeta';
                Kp{c1,c2} = kp';
                chisquare{c1,c2} = chisq';
            case 'svd'
                Ca{c1,c2} = Calpha;
                Cb{c1,c2} = Cbeta;
                fpa{c1,c2} = fpalpha';
                fpb{c1,c2} = fpbeta';
                fa{c1,c2} = falpha';
                fb{c1,c2} = fbeta';
                Kp{c1,c2} = kp';
                chisquare{c1,c2} = chisq';
        end
    end
end

disp(sprintf('done!'));
fit_output = struct('crit_pts',cpts,'Ca',Ca,'Cb',Cb,'fpa',fpa,'fpb',fpb,'fa',fa,'fb',fb,'Kp',Kp,'chisq',chisquare);
return