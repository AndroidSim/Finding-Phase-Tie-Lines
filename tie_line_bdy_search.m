function fit_output = tie_line_bdy_search(S_traj,C_traj,S_bdy,C_bdy,varargin)

disp(sprintf('processing input arguments:'));

if all(size(S_traj) > 1) & ndims(S_traj) == 2 % if S is a matrix
    [ntb,ncol] = size(S_traj);

    if ncol < 2 || rem(ncol,2) ~= 0
        error('each trajectory spectrum is two columns: [B-field intensity_values]');
    end
    
    nts = ncol/2;
else
    error('first argument must be a matrix of the trajectory spectra');
end

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

if all(size(S_bdy) > 1) & ndims(S_bdy) == 2 % if S is a matrix
    [nbb,ncol] = size(S_bdy);

    if ncol < 2 || rem(ncol,2) ~= 0
        error('each boundary spectrum is two columns: [B-field intensity_values]');
    end
    
    nbs = ncol/2;
else
    error('third argument must be a matrix of the boundary spectra');
end

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

% for each trajectory spectrum, search around the boundary and find the
% best fit pair of basis spectra.  for each pair of basis spectra searched
% the fraction of each phase can be calculated and just the kp is searched.

if isempty(varargin)
    B_search_i = 0.01; % increment for Basis search along boundary (100 points)
    B_eval_method = 'lls'; % default basis evaluation method, 'lls' = linear least squares
else
    if isempty(varargin{1})
        B_search_i = 0.01;
    else
        if isnumeric(varargin{1})
            if varargin{1} <= 0 & varargin{1} >= 1;
                error('Basis search parameter must be a number between 0 and 1');
            else
                B_search_i = varargin{1};
            end
        else
            error('Basis search parameter must be a number between 0 and 1');
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

% boundary pre-processing
% bdy = make_convex_bdy(bdy);
[len_bdy,int_bdy] = bdy_fxn(C_bdy,'linear');

% spectra pre-processing
S_bdy = align_spectra(S_bdy);

ops = optimset('display','off','largescale','off');
options = optimset('display','off','jacobian','on');
Bp = [0:B_search_i:1];
nBp = length(Bp);

disp(sprintf('starting search:'));
% loop over trajectory spectra
for t = 1:nts
    St = S_traj(:,t*2-1:t*2);
    Ct = C_traj(t,:);
    
    % start with first boundary point and increment a grid search from
    % there using 'B_search_i'.
    for b = 1:nBp
        % first incremented boundary point will be the alpha phase
        % composition (called C_alpha)
        % C_alpha = interp_boundary(C_bdy,b,'pchip')
        C_alpha = interp1(int_bdy(:,1),C_bdy,Bp(b),'linear');
        % interpolate spectra at alpha point
        S_alpha = interp_spectra(S_bdy,C_bdy,'curve',linspace(S_bdy(1,1),S_bdy(end,1),ntb)',[],C_alpha,'linear');
        
        % the intersection of the line drawn through C_alpha and
        % trajectory point with the boundary is the beta phase composition
        % (called the beta point)
        int_pts = boundary_line_int(C_bdy,[C_alpha;Ct]);
        nipt = size(int_pts,1);
           
        if isempty(int_pts) % no intersection
            error('tie line does not intersect boundary');
        elseif nipt == 2 
            if isequal(round(int_pts(1,:).*1000).*1000,round(C_alpha.*1000).*1000)
                C_beta = int_pts(2,:); 
            else
                C_beta = int_pts(1,:);
            end
            
            % interpolate spectra at beta point
            S_beta = interp_spectra(S_bdy,C_bdy,'curve',linspace(S_bdy(1,1),S_bdy(end,1),ntb)',[],C_beta,'linear');
        elseif nipt > 2 
            error('tie line intersects boundary at > 2 points');
        end
    
        % calculate fraction of alpha phase and fraction of beta phase
        fbeta(b) = norm(Ct-C_alpha)/norm(C_alpha-C_beta);
        falpha(b) = 1-fbeta(b);
        
        % solve linear least squares problem or spectral linear combination
        % problem to determine Kp
        % kp defined as into beta (Lo) phase:
        %   fpalpha = fa./(fa+kp.*fb);
        %   fpbeta = (kp.*fb)./(fa+kp.*fb);
        Basis = [S_alpha S_beta]; % define Basis
        alinedS = align_spectra([St Basis]);
        ISt = alinedS(:,1:2);
        Basis = alinedS(:,3:end);
        
        % reduced Basis and trajectory point spectrum to just intensity values
        Basis = Basis(:,2:2:end);
        ISt = ISt(:,2);
        
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
                fpalpha(b) = x(1);
                fpbeta(b)= x(2);
                chisq(b) = sqrt(resnorm);
                kp(b) = (falpha(b)*(1-fpalpha(b)))/(fpalpha(b)*fbeta(b));
            case 'spectral'
                [x,resnorm] = lsqcurvefit(@lincomb_fit_fun,1.01,fbeta(b),ISt,0.1,10,options,'spectral','spectrum_fit',Basis); 
                chisq(b) = sqrt(resnorm);
                kp(b) = x;
            case 'svd'
                % Aeq = [1 1];
                % beq = 1;
                % [U,V,X,C,S] = gsvd(Basis,Aeq,0);
                % W = inv(X');
                % x = ((U(:,2)'*ISt)/C(2,2))*W(:,2)+((V(:,1)'*beq)/S(1,1))*W(:,1);
  
                [U,W,V] = svd(Basis,0);
                W = diag(W);
                x(1) = sum(((U(:,1)'*ISt)/W(1))*V(:,1));
                x(2) = sum(((U(:,2)'*ISt)/W(2))*V(:,2));
                xt = sum(x);
                if (x(1) > 0 & x(2) > 0) | (x(1) < 0 & x(2) < 0)
                    x = x./xt;
                elseif x(1) < 0 & x(2) > 0 & xt > 0
                    x = [(x(1)+xt)/xt; (x(2)-xt)/xt];
                elseif x(1) < 0 & x(2) > 0 & xt < 0
                    x = [(x(1)-xt)/xt; (x(2)+xt)/xt];
                elseif x(1) > 0 & x(2) < 0 & xt > 0
                    x = [(x(1)-xt)/xt; (x(2)+xt)/xt];
                elseif x(1) > 0 & x(2) < 0 & xt < 0
                    x = [(x(1)+xt)/xt; (x(2)-xt)/xt];
                end
                fpalpha(b) = x(1);
                fpbeta(b)= x(2);
                chisq(b) = norm((Basis*x)-ISt);
                kp(b) = (falpha(b)*(1-fpalpha(b)))/(fpalpha(b)*fbeta(b));
        end
        
        Calpha(b,:) = C_alpha;
        Cbeta(b,:) = C_beta;
    end
    
    switch B_eval_method
        case 'lls'
            Ca{t} = Calpha;
            Cb{t} = Cbeta;
            fpa{t} = fpalpha';
            fpb{t} = fpbeta';
            fa{t} = falpha';
            fb{t} = fbeta';
            Kp{t} = kp';
            chisquare{t} = chisq';
        case 'spectral'
            Ca{t} = Calpha;
            Cb{t} = Cbeta;
            fpa{t} = (falpha')./(falpha'+(kp'.*fbeta'));
            fpb{t} = (kp'.*fbeta')./(falpha'+(kp'.*fbeta'));
            fa{t} = falpha';
            fb{t} = fbeta';
            Kp{t} = kp';
            chisquare{t} = chisq';
        case 'svd'
            Ca{t} = Calpha;
            Cb{t} = Cbeta;
            fpa{t} = fpalpha';
            fpb{t} = fpbeta';
            fa{t} = falpha';
            fb{t} = fbeta';
            Kp{t} = kp';
            chisquare{t} = chisq';
    end
end

disp(sprintf('done!'));
fit_output = struct('Ca',Ca,'Cb',Cb,'fpa',fpa,'fpb',fpb,'fa',fa,'fb',fb,'Kp',Kp,'chisq',chisquare);
return 