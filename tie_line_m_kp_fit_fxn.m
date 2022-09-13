function fit_output = tie_line_m_kp_fit_fxn(slope,Kp,Ct,St,C_bdy,S_bdy,search_method)
% tie_line_pt_m_fit_fxn is only called from tie_line_search fxn

[ntb,ntcol] = size(St);
[ntc,ntd] = size(Ct);
[nbb,nbcol] = size(S_bdy);
[nbc,nbd] = size(C_bdy);

options = optimset('display','off');

% the intersection of the line drawn through the
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
  
if Kp == 0
    fit_output = [C_alpha;C_beta];
    return;
end

% calculate fraction of alpha phase and fraction of beta phase
fbeta = norm(Ct-C_alpha)/norm(C_alpha-C_beta);
falpha = 1-fbeta;

% interpolate spectra at alpha point on the boundary
S_alpha = interp_spectra(S_bdy,C_bdy,'curve',[],[],C_alpha,'linear');
    
% interpolate spectra at beta point on the boundary
S_beta = interp_spectra(S_bdy,C_bdy,'curve',[],[],C_beta,'linear');

% kp defined as into alpha (La) phase:
%   fpalpha = (kp.*fa)./(fb+kp.*fa);
%   fpbeta = fb./(fb+kp.*fa);
Basis = [S_alpha S_beta]; % define Basis
        
% reduced Basis and trajectory point spectrum to just intensity values
Basis = Basis(:,2:2:end);
ISt = St(:,2);

Calpha = C_alpha;
Cbeta = C_beta;
Salpha = S_alpha(:,2);
Sbeta = S_beta(:,2);

% calculate regression: the fit is the square root
% of the difference between the data points and the best fit fpalpha or
% fpbeta curve as a fxn of Kp and falpha and fbeta at each magnetic field.
% warning off;
% % loop over magnetic fields
% for b = 1:size(Basis,1)
%     x = [1 falpha 0]';
%     y = [Basis(b,1) ISt(b) Basis(b,2)]';
%     [p,chisq] = lsqcurvefit(@I_fit_fxn,[1 1],x,y,[],[],options,Kp);
%     r(b) = chisq;
% end
% warning on;
% 
% chisq = norm(r)/length(r);

fpalpha = (Kp.*falpha)./(fbeta+Kp.*falpha);
fpbeta = fbeta./(fbeta+Kp.*falpha);
IStc = fpalpha*Salpha + fpbeta*Sbeta;
% chisq = (1-dot(ISt,IStc)/(norm(ISt)*norm(IStc)));% + dot(Salpha,Sbeta)/(norm(Salpha)*norm(Sbeta));
chisq = norm(ISt-IStc);

% ops = optimset('display','off','largescale','off');
% [U,W,V] = svd(Basis,0);
% B = U(:,1:2); % B = 2 component eigenbasis
% % solve for eigenspectra coefficient matrix M by linear least squares
% % constraints
% % Aeq = [sum(U(:,1)) sum(U(:,2))];
% Aeq = sum(B);
% beq1 = sum(Basis(:,1));
% beq2 = sum(Basis(:,2));
% %A = -B;
% %b = zeros(size(U(:,1)));
% temp1 = lsqlin(B,Basis(:,1),[],[],Aeq,beq1,[],[],[],ops);
% temp2 = lsqlin(B,Basis(:,2),[],[],Aeq,beq2,[],[],[],ops);
% M = [temp1 temp2];
% Aeq = [1 1];
% beq = 1;
% x = lsqlin(M,[ISt'*B(:,1);ISt'*B(:,2)],[],[],Aeq,beq,[],[],[],ops);
% %z = inv(M)*[ISt'*B(:,1);ISt'*B(:,2)];
% 
% % [U,W,V] = svd(Basis,0);
% % W = diag(W);
% % x(1) = sum(((U(:,1)'*ISt)/W(1))*V(:,1));
% % x(2) = sum(((U(:,2)'*ISt)/W(2))*V(:,2));
% 
% fpalphac = x(1);
% fpbetac = x(2);
% chisq = sqrt((fpalpha-fpalphac)^2 + (fpbeta-fpbetac)^2);

% Kp = ((Xsa.^2).*(log(Xsa)-0.5) + (Xca.^2).*(log(Xca)-0.5))./((Xsb.^2).*(log(Xsb)-0.5) + (Xcb.^2).*(log(Xcb)-0.5))

% output according to search method
switch search_method
    case 'continuous'
        fit_output = vIStc;
    case 'grid'
        Ca = Calpha;
        Cb = Cbeta;
        Sa = Salpha;
        Sb = Sbeta;
        fafb = [falpha fbeta];
        fit_output = struct('Ca',Ca,'Cb',Cb,'Sa',Sa,'Sb',Sb,'fafb',fafb,'chisq',chisq);
    case 'mc'
        Ca = Calpha;
        Cb = Cbeta;
        Sa = Salpha;
        Sb = Sbeta;
        fafb = [falpha fbeta];
        fit_output = struct('Ca',Ca,'Cb',Cb,'Sa',Sa,'Sb',Sb,'fafb',fafb,'chisq',chisq);
    case 'mcmc'
        Ca = Calpha;
        Cb = Cbeta;
        Sa = Salpha;
        Sb = Sbeta;
        fafb = [falpha fbeta];
        fit_output = struct('Ca',Ca,'Cb',Cb,'Sa',Sa,'Sb',Sb,'fafb',fafb,'chisq',chisq);
    otherwise 
        error('invalid search method');
end

return