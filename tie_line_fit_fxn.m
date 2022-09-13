function fit_output = tie_line_fit_fxn(p,C_traj,S_traj,C_bdy,S_bdy,bdy_p,search_method,slope_method,fit_method,refS,metric)
% fxn tie_line_fit_fxn only called from fxn tie_line_search

[ntb,ntcol] = size(S_traj);
[ntc,ntd] = size(C_traj);
[nbb,nbcol] = size(S_bdy);
[nbc,nbd] = size(C_bdy);

reference_S = false;

if ~isempty(refS)
    ref_S_alpha = refS(:,1:2);
    ref_S_beta = refS(:,3:4);
    [rab,rac] = size(ref_S_alpha);
    [rbb,rbc] = size(ref_S_beta);
    reference_S = true;
end

ops = optimset('display','off','largescale','off');
options = optimset('display','off','jacobian','on');

switch slope_method
    case 'tangent'
        cpt1 = interp1(bdy_p,C_bdy,p(1),'linear');
        % cpt1 = interp_boundary(C_bdy,crit_p(c1),'linear')
        cpt2 = interp1(bdy_p,C_bdy,p(2),'linear');
        % cpt2 = interp_boundary(C_bdy,crit_p(c2),'linear')
        % find slope of tangent line at critical points
        temp1 = [p(1)-0.02:0.02:p(1)+0.02]';
        if any(temp1 < 0)
            temp1(temp1 < 0) = temp1(temp1 < 0)+1;
        end
        if any(temp1 > 1)
            temp1(temp1 > 1) = temp1(temp1 > 1)-1;
        end
        ptacpt1 = interp1(bdy_p,C_bdy,temp1,'linear')'; % pts around cpt1
        mc1 = ppval(fnder(spline(temp1,ptacpt1)),p(1)); % derivative of spline through pts
        warning off;
        mc1 = mc1(2)/mc1(1); % slope of tangent line through cpt1
        warning on;
        
        temp2 = [p(2)-0.02:0.02:p(2)+0.02]';
        if any(temp2 < 0)
            temp2(temp2 < 0) = temp2(temp2 < 0)+1;
        end
        if any(temp2 > 1)
            temp2(temp2 > 1) = temp2(temp2 > 1)-1;
        end
        ptacpt2 = interp1(bdy_p,C_bdy,temp2,'linear')'; % pts around cpt2
        mc2 = ppval(fnder(spline(temp2,ptacpt2)),p(2)); % derivative of spline through pts
        warning off;
        mc2 = mc2(2)/mc2(1); % slope of tangent line through cpt2
        warning on;
        clear temp1 temp2;

        % represent tangent line in 2 pt form
        if isequal(mc1,Inf) | isequal(mc1,-Inf) % mc1 == Inf or -Inf
            lina = [cpt1;cpt1(1) 1];
        elseif mc1 > 0
            lina = [cpt1;1 mc1*(1-cpt1(1))+cpt1(2)];
        elseif isequal(mc1,0) % mc1 == 0
            lina = [cpt1;1 cpt1(2)];
        elseif mc1 < 0 
            lina = [cpt1;0 mc1*(0-cpt1(1))+cpt1(2)];
        end
        
        if isequal(mc2,Inf) | isequal(mc2,-Inf) % mc2 == Inf or -Inf
            linb = [cpt2;cpt2(1) 1];
        elseif mc2 > 0
            linb = [cpt2;1 mc2*(1-cpt2(1))+cpt2(2)];
        elseif isequal(mc2,0) % mc2 == 0
            linb = [cpt2;1 cpt2(2)];
        elseif mc2 < 0 
            linb = [cpt2;0 mc2*(0-cpt2(1))+cpt2(2)];
        end

        % find intersection pt of tangent lines
        int_pt = line_line_int(lina,linb);
        clear lina linb;
        
        if isempty(int_pt) % tangent lines parallel
            t_slopes = repmat(mc1,ntc,1); % or mc2
            t_angles = slope2angle(t_slopes);
        else 
            t_angles = pt2angle(C_traj-repmat(int_pt,ntc,1));
            t_slopes = angle2slope(t_angles');
        end
        
    case 'parallel'
        t_slopes = repmat(p,ntc,1); 
        t_angles = slope2angle(t_slopes);
        int_pt = [];
    otherwise
        error('invalid slope method');
end

% loop over trajectory points 
for t = 1:ntc
    St = S_traj(:,t*2-1:t*2);
    Ct = C_traj(t,:);
    
    % the intersection of the tie line drawn through the intersection
    % point of the critical point tangent lines and the
    % trajectory point with the boundary gives both the composition
    % point of the alpha phase (C_alpha, or Ld phase) and the
    % composition of the beta phase (C_beta, or Lo phase).  the
    % points returned from boundary_line_int are that
    % C_alpha=int_pts(1,:) and C_beta=int_pts(2,:)
    int_pts = boundary_line_int(C_bdy,t_slopes(t),Ct,int_pt);
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
           
    % interpolate spectra at alpha point on the boundary
    S_alpha = interp_spectra(S_bdy,C_bdy,'curve',[],[],C_alpha,'linear');
%     if ~isnorm(S_alpha(:,2))
%         S_alpha(:,2) = normalize_vector(S_alpha(:,2));
%     end

    % interpolate spectra at beta point on the boundary
    S_beta = interp_spectra(S_bdy,C_bdy,'curve',[],[],C_beta,'linear');
%     if ~isnorm(S_beta(:,2))
%         S_beta(:,2) = normalize_vector(S_beta(:,2));
%     end
        
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
        
    switch fit_method
        case 'lls'
            % linear least squares with constraints
            % solve Ax = b for fraction of probe in alpha phase and
            % then solve Kp equation as a fxn of fraction of probe in
            % alpha phase and fraction of alpha phase
            Aeq = [1 1];
            beq =1;
            [x,resnorm] = lsqlin(Basis,ISt,[],[],Aeq,beq,0.00001,0.99999,[],ops);% x = lsqlin(C,d,A,b,Aeq,beq,lb,ub,x0,options)
            if any(x < 0) | any(x > 1)
                keyboard;
            end
%             Aeq = [1 1];
%             beq = 1;
%             Aneq = [-1 0;0 -1;1 0;0 1];
%             bneq = [0.00001;0.00001;0.99999;0.99999];
%             [x,resnorm] = lsqlin(Basis,ISt,Aneq,bneq,Aeq,beq,[],[],[],ops);% x = lsqlin(C,d,A,b,Aeq,beq,lb,ub,x0,options)
            fpalpha(t) = x(1);
            fpbeta(t)= x(2);
            IStc(:,t) = Basis*x;
            kp(t) = (falpha(t)*(1-fpalpha(t)))/(fpalpha(t)*fbeta(t));
        case 'spectral'
            [x,resnorm] = lsqcurvefit(@lincomb_fit_fun,1.01,fbeta(t),ISt,0.1,10,options,'spectral','spectrum_fit',Basis);
            IStc(:,t) = feval(@lincomb_fit_fun,x,fbeta(t),'spectral','spectrum_fit',Basis);
            fpalpha(t) = falpha(t)/(falpha(t)+(x*fbeta(t))); % x = Kp
            fpbeta(t) = (x*fbeta(t))/(falpha(t)+(x*fbeta(t)));
            kp(t) = x;
        case 'svd'
%             Aeq = [1 1];
%             beq = 1;
%             [U,V,X,C,S] = gsvd(Basis,Aeq,0);
%             W = inv(X');
%             x = ((U(:,2)'*ISt)/C(2,2))*W(:,2)+((V(:,1)'*beq)/S(1,1))*W(:,1);
%   
%             [U,W,V] = svd(Basis,0);
%             W = diag(W);
%             x(1) = sum(((U(:,1)'*ISt)/W(1))*V(:,1));
%             x(2) = sum(((U(:,2)'*ISt)/W(2))*V(:,2));
%             xt = sum(x);
%                     if (x(1) > 0 & x(2) > 0) | (x(1) < 0 & x(2) < 0)
%                         x = x./xt;
%                     elseif x(1) < 0 & x(2) > 0 & xt > 0
%                         x = [(x(1)+xt)/xt; (x(2)-xt)/xt];
%                     elseif x(1) < 0 & x(2) > 0 & xt < 0
%                         x = [(x(1)-xt)/xt; (x(2)+xt)/xt];
%                     elseif x(1) > 0 & x(2) < 0 & xt > 0
%                         x = [(x(1)-xt)/xt; (x(2)+xt)/xt];
%                     elseif x(1) > 0 & x(2) < 0 & xt < 0
%                         x = [(x(1)+xt)/xt; (x(2)-xt)/xt];
%                     end
%                     
%             convert to derivative spectra to absorbance spectra
%             Basis = deriv2abs(Basis);
%             St = deriv2abs(St);
%             Basis = normalize_spectra(Basis,'a');
%             St = normalize_spectra(St,'a');
%             
%             perform singular value decomposition to obtain eigenspectra basis
%             disp(sprintf('performing singular value decomposition to obtain eigenspectra basis:\n'));

%           normalize spectra so that the sum of intensities equals one?
    %       Basis(:,1) = Basis(:,1)./sum(Basis(:,1));
    %       Basis(:,2) = Basis(:,2)./sum(Basis(:,2));
    %       ISt = ISt./sum(ISt);
    
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
            kp(t) = (falpha(t)*(1-fpalpha(t)))/(fpalpha(t)*fbeta(t));
            IStc(:,t) = Basis*x; 
        otherwise
            error('invalid fit method');
    end % end of switch/case for fit methods

    Calpha(t,:) = C_alpha;
    Cbeta(t,:) = C_beta;
    Salpha(:,t) = S_alpha(:,2);
    Sbeta(:,t) = S_beta(:,2);
end % end of loop over trajectory points 

% calculate distance: each magnetic field for both the data trajectory 
% points and the calculated trajectory points has a intenstiy surface
% as a function of composition.  the chisquare is the square root
% of the difference between all intensity surfaces of the data
% trajectory points (ISt) and the calculated trajectory points
% (IStc).

if reference_S
	IrSa = ref_S_alpha(:,2);
	IrSb = ref_S_beta(:,2);
	ISt = S_traj(:,2:2:end);
	[ndb,ndc] = size(ISt);
	[ncb,ncc] = size(IStc);
	[ab,ac] = size(Salpha);
	[bb,bc] = size(Sbeta);
	vISt = reshape(ISt,ndb*ndc,1);
	vIStc = reshape(IStc,ndb*ncc,1);
	vSa = reshape(Salpha,ab*ac,1);
	vSb = reshape(Sbeta,bb*bc,1);
    NSt = ntc; % number of trajectory spectra = ndc = ncc
    NSa = ac; % number of alpha basis spectra
    NSb = bc; % number of beta basis spectra
    ref_dist = (1/NSa)*pdist([vSa repmat(IrSa,ac,1)]',metric)+(1/NSb)*pdist([vSb repmat(IrSb,bc,1)]',metric);
    data_dist = (1/NSt)*pdist([vISt vIStc]',metric);
    total_dist = data_dist+ref_dist;
else
    ISt = S_traj(:,2:2:end);
	[ndb,ndc] = size(ISt);
	[ncb,ncc] = size(IStc);
	vISt = reshape(ISt,ndb*ndc,1);
	vIStc = reshape(IStc,ndb*ncc,1);
    NSt = ntc; % number of trajectory spectra = ndc = ncc
    ref_dist = [];
    data_dist = (1/NSt)*pdist([vISt vIStc]',metric);
    total_dist = data_dist;
end

% output according to search method
switch search_method
    case 'continuous'
        fit_output = vIStc;
    case 'grid'
        Ca = Calpha;
        Cb = Cbeta;
        Sa = Salpha;
        Sb = Sbeta;
        fpafpb = [fpalpha' fpbeta'];
        fafb = [falpha' fbeta'];
        Kp = kp';
        fit_output = struct('Ca',Ca,'Cb',Cb,'Sa',Sa,'Sb',Sb,'fpafpb',fpafpb,'fafb',fafb,'Kp',kp,'total_dist',total_dist,'data_dist',data_dist,'ref_dist',ref_dist);
    otherwise 
        error('invalid search method');
end

return