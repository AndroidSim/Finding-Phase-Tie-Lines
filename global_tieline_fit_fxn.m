function fit_output = global_tieline_fit_fxn(cp,et,C_coex,S_coex,C_bdy,S_bdy,Bdy,Bdy_p,Cbdy_p,slope_method,kp_method,fit_method,kp,var_cSbline,varargin)
% fit_output =
% global_tieline_fit_fxn(cp,et,C_coex,S_coex,C_bdy,S_bdy,Bdy,Bdy_p,Cbdy_p,slope_method,kp_method,fit_method,kp,var_cSbline,varargin)

[ncB,ncScol] = size(S_coex);
[ncC,ncd] = size(C_coex);
[nbB,nbScol] = size(S_bdy);
[nbC,nbd] = size(C_bdy);

B = S_coex(:,1); % or S_bdy, magnetic field values the same for all spectra because of alignment
% reduced all spectra to a vector of intensity values (needed for fit)
S_coex = S_coex(:,2:2:end);
S_bdy = S_bdy(:,2:2:end);
[ncB,ncS] = size(S_coex);
[nbB,nbS] = size(S_bdy);

if ncS ~= ncC
    error('the number of coex spectra does not equal the number of coex comps');
end

[nBp,Bpd] = size(Bdy);
if length(Bdy_p) ~= nBp
    error('length of parameter vector for the boundary points does not equal number boundary points in global_tieline_fit_fxn');
end

switch slope_method
    case 'tangent'
        spx = varargin{1};
        spy = varargin{2};
        if isempty(cp)
            if iscell(et)
                eta = et{1};
                etb = et{2};
                % 2 end tielines
                intersection = line_line_int(eta,etb);
                if isempty(intersection) % lines parallel
                    m = diff(eta); % or etb, should be the same
                    m = m(2)/m(1);
                    slopes = repmat(m,ncC,1); 
                else 
                    angles = pt2angle(C_coex-repmat(intersection,ntc,1));
                    slopes = angle2slope(angles);
                end
            else
                error('no critical point and end tielines not passed as cell array to global_tieline_fit_fxn');
            end
        else
            if isscalar(cp)
                cpt = param2bdypt(Bdy,cp,'spline',spx,spy);
                warning off;
                mc = fnval(fnder(spy),cp)/fnval(fnder(spx),cp);
                warning on;

                % represent tangent line in 2 pt form
                if isequal(mc,Inf) | isequal(mc,-Inf) % mc1 == Inf or -Inf
                    lin = [cpt;cpt(1) 1];
                elseif mc > 0
                    lin = [cpt;1 mc*(1-cpt(1))+cpt(2)];
                elseif isequal(mc,0) % mc1 == 0
                    lin = [cpt;1 cpt(2)];
                elseif mc < 0 
                    lin = [cpt;0 mc*(0-cpt(1))+cpt(2)];
                end

                if isnumeric(et) & size(et) == [2 2]
                    % find intersection pt of tangent line and end tieline
                    intersection = line_line_int(lin,et);
                    if isempty(intersection) % lines parallel
                        slopes = repmat(mc,ncC,1); % or mc2
                    else 
                        angles = pt2angle(C_coex-repmat(intersection,ncC,1));
                        slopes = angle2slope(angles);
                    end
                end
                clear lin mc;
            elseif isvector(cp)
                cpt1 = param2bdypt(Bdy,cp(1),'spline',spx,spy);
                warning off;
                mc1 = fnval(fnder(spy),cp(1))/fnval(fnder(spx),cp(1));
                warning on;
                cpt2 = param2bdypt(Bdy,cp(2),'spline',spx,spy);
                warning off;
                mc2 = fnval(fnder(spy),cp(2))/fnval(fnder(spx),cp(2));
                warning on;

                % represent tangent line in 2 pt form
                if isequal(mc1,Inf) | isequal(mc1,-Inf) % mc1 == Inf or -Inf
                    lina = [cpt1;cpt1(1) 1];
                elseif mc1 > 0
                    lina = [cpt1;1 mc1*(1-cpt1(1))+cpt1(2)];
                elseif isequal(mc1,0) % mc1 == 0
                    lina = [cpt1;1 cpt1(2)];
                elseif mc1 < 0 
                    lina = [cpt1;0 mc1*(0-cpt1(1))+cpt1(2)];
                elseif isnan(mc1)
                    if fnval(fnder(spy),cp(1)) == Inf | fnval(fnder(spy),cp(1)) == -Inf
                        lina = [cpt1;cpt1(1) 1];
                    end
                    if fnval(fnder(spy),cp(1)) == 0
                        lina = [cpt1;1 cpt1(2)];
                    end 
                end

                if isequal(mc2,Inf) | isequal(mc2,-Inf) % mc2 == Inf or -Inf
                    linb = [cpt2;cpt2(1) 1];
                elseif mc2 > 0
                    linb = [cpt2;1 mc2*(1-cpt2(1))+cpt2(2)];
                elseif isequal(mc2,0) % mc2 == 0
                    linb = [cpt2;1 cpt2(2)];
                elseif mc2 < 0 
                    linb = [cpt2;0 mc2*(0-cpt2(1))+cpt2(2)];
                elseif isnan(mc2)
                    if fnval(fnder(spy),cp(2)) == Inf | fnval(fnder(spy),cp(2)) == -Inf
                        linb = [cpt2;cpt2(1) 1];
                    end
                    if fnval(fnder(spy),cp(2)) == 0
                        linb = [cpt2;1 cpt2(2)];
                    end 
                end

                % find intersection pt of tangent lines
                intersection = line_line_int(lina,linb);
                if isempty(intersection) % lines parallel
                    slopes = repmat(mc1,ncC,1); % or mc2
                else 
                    angles = pt2angle(C_coex-repmat(intersection,ncC,1));
                    slopes = angle2slope(angles);
                end
                clear lina linb mc1 mc2;
            else
                error('critical point boundary parameters passed to global_tieline_fit_fxn has to be a scalar or vector');
            end
        end
    case 'parallel'
        slopes = repmat(cp,ncC,1); 
    otherwise
        error('invalid slope method in the global tieline fitting function');
end

options = optimset('display','off','largescale','on');
% loop over phase coexistence points 
for c = 1:ncC
    cS = S_coex(:,c);
    cC = C_coex(c,:);
    
    % the intersection of the tie line drawn through the intersection
    % point of the critical point tangent line(s) and/or end tielines and the
    % coexistence point with the boundary gives both the composition
    % point of the alpha phase (C_alpha, or Ld phase) and the
    % composition of the beta phase (C_beta, or Lo phase).  the
    % points returned from boundary_line_int are that
    % C_alpha=int_pts(1,:) and C_beta=int_pts(2,:)
    int_pts = boundary_line_int(Bdy,slopes(c),cC);
    nipt = size(int_pts,1);
            
    if isempty(int_pts) % no intersection
        error('tie line does not intersect boundary');
    elseif isequal(nipt,1) % tangent intersection
        error('tie line intersects boundary at only 1 point');
    elseif isequal(nipt,2) % 2 intersections
        aC = int_pts(1,:);
        bC = int_pts(2,:);
    elseif nipt > 2 
        error('tie line intersects boundary at > 2 points');
    end
    
    switch slope_method
        case 'tangent'
            % interpolate spectra at alpha point on the boundary
            ap = bdypt2param(Bdy,aC,'spline',spx,spy);
        %     ap = bdypt2b(Bdy,aC);
            aS = interp1(Cbdy_p,S_bdy',ap,'linear');
            aS = aS';   
        %     S_alpha = interp_spectra(S_bdy,C_bdy,'curve',[],[],C_alpha,'linear');
            % normalize spectrum
            aS = normalize_spectra([B aS],'d');
            aS = aS(:,2);
            % interpolate spectra at beta point on the boundary
            bp = bdypt2param(Bdy,bC,'spline',spx,spy);
        %     bp = bdypt2b(Bdy,bC);
            bS = interp1(Cbdy_p,S_bdy',bp,'linear');
            bS = bS'; 
        %     S_beta = interp_spectra(S_bdy,C_bdy,'curve',[],[],C_beta,'linear');
            % normalize spectrum
            bS = normalize_spectra([B bS],'d');
            bS = bS(:,2);
        case 'parallel'
            % interpolate spectra at alpha point on the boundary
            ap = bdypt2param(Bdy,aC,'polygon',Bdy_p);
        %     ap = bdypt2b(Bdy,aC);
            aS = interp1(Cbdy_p,S_bdy',ap,'linear');
            aS = aS';   
        %     S_alpha = interp_spectra(S_bdy,C_bdy,'curve',[],[],C_alpha,'linear');
            % normalize spectrum
            aS = normalize_spectra([B aS],'d');
            aS = aS(:,2);
            % interpolate spectra at beta point on the boundary
            bp = bdypt2param(Bdy,bC,'polygon',Bdy_p);
        %     bp = bdypt2b(Bdy,bC);
            bS = interp1(Cbdy_p,S_bdy',bp,'linear');
            bS = bS'; 
        %     S_beta = interp_spectra(S_bdy,C_bdy,'curve',[],[],C_beta,'linear');
            % normalize spectrum
            bS = normalize_spectra([B bS],'d');
            bS = bS(:,2);
    end
    
    % define Basis
    Basis = [aS bS]; 
    condn(c,1) = cond(Basis);
    
    % solve linear system of eqns by svd or spectral linear combination
    % kp defined as into beta (Lo) phase:
    %   fpalpha = fa./(fa+kp.*fb);
    %   fpbeta = (kp.*fb)./(fa+kp.*fb);
    
    % calculate fraction of alpha phase and fraction of beta phase
    fb = norm(cC-aC)/norm(aC-bC);
    fa = 1-fb;
    
    % fit coexistence spectrum with basis spectra according to the
    % fit_method and the kp_method
    switch fit_method
        case 'chisq'
            switch kp_method
                case 'constant'
                    % if Kp is given in some form, explicit or function,
                    % then calculate chisq between the convex combination of the basis
                    % spectra and the experimental spectrum.
                    fpb = (kp*fb)/(fa + kp*fb);
                    fpa = fa/(fa + kp*fb);
                    eSc = fpa*aS + fpb*bS;
                    oSc = cS;
                    Kp(c,1) = kp; 
                case 'search'
                    % if Kp not given, then it is a variable parameter in the fitting between 
                    % the convex combination of the basis spectra and the experimental spectrum.
                    try
                        [x,chisq] = lsqcurvefit(@lincomb_fit_fun,1.5,fb,cS,1,2,options,'spectrum','lsqcurvefit',Basis);
                    catch
                        keyboard;
                    end
                    eSc = feval(@lincomb_fit_fun,x,fb,'spectrum','lsqcurvefit',Basis);
                    oSc = cS;
                    fpb = (x*fb)/(fa + x*fb);% x = Kp
                    fpa = fa/(fa + x*fb); 
                    Kp(c,1) = x;
                case 'fxn'
                    fpb = (kp*fb)/(fa + kp*fb);
                    fpa = fa/(fa + kp*fb);
                    eSc = fpa*aS + fpb*bS;
                    oSc = cS;
                    Kp(c,1) = kp;
                otherwise
                    error('invalid kp method');
            end
            fpalpha(c,1) = fpa;
            fpbeta(c,1) = fpb;
%             var_eSbline = var(eSc([1:100 end-100:end]'));
%             covmatoe = cov(oSc([1:100 end-100:end]'),eSc([1:100 end-100:end]'));
%             chisquare(c,1) = ((oSc-eSc)'*(oSc-eSc))/((ncB-1)*sum(diag(covmatoe))); 
%             chisquare(c,1) = ((oSc-eSc)'*(oSc-eSc))/((ncB-1)*prod(sqrt(diag(covmatoe))));
            chisquare(c,1) = ((oSc-eSc)'*(oSc-eSc))/((ncB-1)*mean(var_cSbline));
%             chisquare(c,1) = ((oSc-eSc)'*(oSc-eSc));
%             chisquare(c,1) = ((oSc-eSc)'*(oSc-eSc))/((ncB-1)*var_cSbline(c));
%             var_aSbline(c) = var(aS([1:100 end-100:end]'));
%             var_bSbline(c) = var(bS([1:100 end-100:end]'));
%             covmatab = cov(aS([1:100 end-100:end]'),bS([1:100 end-100:end]'));
%             chisquare(c,1) = ((oSc-eSc)'*(oSc-eSc))/((ncB-1)*prod(sqrt(diag(covmatoe)))) - ((aS-bS)'*(aS-bS))/((ncB-1)*prod(sqrt(diag(covmatab))));
        case 'svd'
            % if Kp is given in some form, explicit or function,
            % the calculate chisq between the svd calculated fraction of probe and
            % the expected from the equation fp(fa,fb,Kp).
            options = optimset('largescale','off');
            
%             [U,W,V] = svd(Basis,0);
%             Bu = U(:,1:2); % B = 2 component eigenbasis
%             % solve for eigenspectra coefficient matrix M by linear least squares
%             % constraints
%             % Aeq = [sum(U(:,1)) sum(U(:,2))];
%             Aeq = sum(Bu);
%             beq1 = sum(Basis(:,1));
%             beq2 = sum(Basis(:,2));
%             %A = -B;
%             %b = zeros(size(U(:,1)));
%             temp1 = lsqlin(Bu,Basis(:,1),[],[],Aeq,beq1,[],[],[],options);
%             temp2 = lsqlin(Bu,Basis(:,2),[],[],Aeq,beq2,[],[],[],options);
%             M = [temp1 temp2];
%             Aeq = [1 1];
%             beq = 1;
%             x = lsqlin(M,[cS'*Bu(:,1);cS'*Bu(:,2)],[],[],Aeq,beq,[],[],[],options);
%             % x = inv(M)*[ISt'*B(:,1);ISt'*B(:,2)];
%             fpa = x(1);
%             fpb = x(2);
            
            x = lsqlin(Basis,cS,[-1 0;0 -1;1 0;0 1],[0;0;1;1],[1 1],1,[],[],[0.5;0.5],options);
            fpa = x(1);
            fpb = x(2);
%             switch kp_method
%                 case 'constant'
%                 case 'search'    
%                 case 'fxn'
%                 otherwise
%                     error('invalid kp method');
%             end
            kp = (fa*fpb)/(fpa*fb);
            eSc = fpa*aS + fpb*bS;
            oSc = cS;
            fpalpha(c,1) = fpa;
            fpbeta(c,1) = fpb;
            Kp(c,1) = kp; 
            chisquare(c,1) = ((oSc-eSc)'*(oSc-eSc))/((ncB-1)*mean(var_cSbline));
        otherwise
            switch kp_method
                case 'constant'
                    fpb = (kp*fb)/(fa + kp*fb);
                    fpa = fa/(fa + kp*fb);
                    eSc = fpa*aS + fpb*bS;
                    oSc = cS;
                    Kp(c,1) = kp;
                case 'search'
                    try
                        [x,chisq] = lsqcurvefit(@lincomb_fit_fun,1.0,fb,cS,0.2,5,options,'spectrum','lsqcurvefit',Basis);
                    catch
                        keyboard;
                    end
                    eSc = feval(@lincomb_fit_fun,x,fb,'spectrum','lsqcurvefit',Basis);
                    oSc = cS;
                    fpb = (x*fb)/(fa + x*fb);% x = Kp
                    fpa = fa/(fa + x*fb); 
                    Kp(c,1) = x;
                case 'fxn'
                otherwise
                    error('invalid kp method');
            end
            fpalpha(c,1) = fpa;
            fpbeta(c,1) = fpb;
            chisquare(c,1) = pdist([eSc';oSc'],fit_method);
%             chisquare(c,1) = pdist([eSc';oSc'],fit_method) - pdist([aS';bS'],fit_method);
    end  
    alphaC(c,:) = aC;
    betaC(c,:) = bC;
    alphaS(:,c) = aS;
    betaS(:,c) = bS;
    expSc(:,c) = eSc;
    obsSc(:,c) = oSc;
    fbeta(c,1) = fb;
    falpha(c,1) = fa;
end % end of loop over phase coexistence points 

% analyze Kp surface
[temp,maC,naC] = unique(alphaC,'rows');
[temp,mbC,nbC] = unique(betaC,'rows');
m = intersect(maC,mbC);
ntls = length(m);
AC = alphaC(m,:);
BC = betaC(m,:);
KP = Kp(m);
fa = linspace(0,1,10)';
for i = 1:ntls
    xtl(:,i) = fa.*AC(i,1) + (1-fa).*BC(i,1);
    ytl(:,i) = fa.*AC(i,2) + (1-fa).*BC(i,2);
end
Kptl = griddata([C_coex(m,1);AC(:,1);BC(:,1)],[C_coex(m,2);AC(:,2);BC(:,2)],[KP;KP;KP],xtl,ytl);
Kpstd = nanmean(nanstd(Kptl));
% calculate distance or fit for all coexistence points: 
ChiSquare = mean(chisquare);

% output fit results 
fpafpb = [fpalpha fpbeta];    
fafb = [falpha fbeta];
fit_output = struct('mag_field',B,'aC',alphaC,'bC',betaC,'aS',alphaS,'bS',betaS,'obs_cS',obsSc,'exp_cS',expSc,'fpafpb',fpafpb,'fafb',fafb,'Kp',Kp,'Kptl',Kptl,'Kpstd',Kpstd,'ChiSquare',ChiSquare,'condn',condn);

return