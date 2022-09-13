function fit_output = tie_line_pt_m_fit_fxn(slope,Ct,St,C_bdy,S_bdy,bdy_p,search_method,slope_method,fit_method,refS,metric)
% tie_line_pt_m_fit_fxn is only called from tie_line_search fxn

[ntb,ntcol] = size(St);
[ntc,ntd] = size(Ct);
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

% the intersection of the tie line drawn through the intersection
% point of the critical point tangent lines and the
% trajectory point with the boundary gives both the composition
% point of the alpha phase (C_alpha, or Ld phase) and the
% composition of the beta phase (C_beta, or Lo phase).  the
% points returned from boundary_line_int are that
% C_alpha=int_pts(1,:) and C_beta=int_pts(2,:)
int_pts = boundary_line_int(C_bdy,slope,Ct,[]);
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
% if ~isnorm(S_alpha(:,2))
%     S_alpha(:,2) = normalize_vector(S_alpha(:,2));
% end
    
% interpolate spectra at beta point on the boundary
S_beta = interp_spectra(S_bdy,C_bdy,'curve',[],[],C_beta,'linear');
% if ~isnorm(S_beta(:,2))
%     S_beta(:,2) = normalize_vector(S_beta(:,2));
% end
    
% calculate fraction of alpha phase and fraction of beta phase
fbeta = norm(Ct-C_alpha)/norm(C_alpha-C_beta);
falpha = 1-fbeta;
        
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
    case 'pt2lin'
        fpalpha = [];
        fpbeta = [];
        IStc = ISt;
        kp = [];
    case 'lls'
        % linear least squares with constraints
        % solve Ax = b for fraction of probe in alpha phase and
        % then solve Kp equation as a fxn of fraction of probe in
        % alpha phase and fraction of alpha phase
%         Nt = norm(ISt);
%         Na = norm(Basis(:,1));
%         Nb = norm(Basis(:,2));
%         Aeq = [1 1;Na Nb];
%         beq = [1;Nt];
%         Aneq = [-1 0;0 -1;1 0;0 1];
%         bneq = [-0.00001;-0.00001;0.99999;0.99999];
%         z = [1 1;Na Nb]\[1;Nt];
        Aeq = [1 1];
        beq = 1;
        [x,resnorm] = lsqlin(Basis,ISt,[],[],Aeq,beq,0.00001,0.99999,[],ops);% x = lsqlin(C,d,A,b,Aeq,beq,lb,ub,x0,options)
        if any(x < 0) | any(x > 1)
            keyboard;
        end
        fpalpha = x(1);
        fpbeta = x(2);
        IStc = Basis*x;
        kp = (falpha*(1-fpalpha))/(fpalpha*fbeta);
    case 'spectral'
        [x,resnorm] = lsqcurvefit(@lincomb_fit_fun,1.01,fbeta,ISt,0.1,10,options,'spectral','spectrum_fit',Basis);
        IStc = feval(@lincomb_fit_fun,x,fbeta,'spectral','spectrum_fit',Basis);
        fpalpha = falpha/(falpha+(x*fbeta)); % x = Kp
        fpbeta = (x*fbeta)/(falpha+(x*fbeta));
        kp = x;
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
            
        fpalpha = x(1);
        fpbeta = x(2);
        kp = (falpha*(1-fpalpha))/(fpalpha*fbeta);
        IStc = Basis*x; 
    otherwise
        error('invalid fit method');
end % end of switch/case for fit methods

Calpha = C_alpha;
Cbeta = C_beta;
Salpha = S_alpha(:,2);
Sbeta = S_beta(:,2);
    
% calculate distance: each magnetic field for both the data trajectory 
% points and the calculated trajectory points has a intensity surface
% as a function of composition.  the chisquare is the square root
% of the difference between all intensity surfaces of the data
% trajectory points (ISt) and the calculated trajectory points
% (IStc).
switch fit_method
    case 'pt2lin'
        warning off;
        for b = 1:size(Basis,1);
            x = [0 fbeta 1]';
            y = [Basis(b,1) ISt(b) Basis(b,2)]';
            R = corrcoef([x y]);
            r(b) = R(1,2)^2;
        end
        
        warning on;
        ref_dist = [];
        data_dist = mean(r(~isnan(r)));
        total_dist = data_dist;
%         for b = 1:size(Basis,1);
%             lin = [Calpha Basis(b,1);Cbeta Basis(b,2)];
%             pt = [Ct ISt(b)];
%             d(b) = pt2lin_dist(pt,lin,3);
%         end
%         
%         ref_dist = [];
%         data_dist = sum(d);
%         total_dist = data_dist;
    otherwise
		if reference_S
			IrSa = ref_S_alpha(:,2);
			IrSb = ref_S_beta(:,2);
			ISt = St(:,2);
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
            ISt = St(:,2);
			[ndb,ndc] = size(ISt);
			[ncb,ncc] = size(IStc);
			vISt = reshape(ISt,ndb*ndc,1);
			vIStc = reshape(IStc,ndb*ncc,1);
            NSt = ntc; % number of trajectory spectra = ndc = ncc
            ref_dist = [];
            data_dist = (1/NSt)*pdist([vISt vIStc]',metric);
            total_dist = data_dist;
		end
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