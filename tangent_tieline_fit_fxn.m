function varargout = tangent_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method)
% varargout =
% tangent_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method)

[ncC,ncd] = size(C_coex);
[ncB,ncS] = size(S_coex);
[nbC,nbd] = size(C_bdy);
[nbB,nbS] = size(S_bdy);

if ncS ~= ncC
    error('the number of coexistence spectra does not equal the number of coexistence compositions');
end

% shorten names from spectra data structure
% spectra = struct('magnetic_fields',magnetic_fields,'cSbline_var',var_cSbline,'cSbline_std',std_cSbline,...
%     'bSC_std',stdbS,'bSC_var',varbS,'weight_fxn_S',weight_fxn_S,'weight_fxn_fp',weight_fxn_fp);
varSc = spectra.cSbline_var;
stdSc = spectra.cSbline_std;
wfxnS = spectra.weight_fxn_S;
wfxnfp = spectra.weight_fxn_fp;
varbS = spectra.bSC_var;
stdbS = spectra.bSC_std;
Nfit = length(find(wfxnS));

% preallocate data arrays
theorycS = zeros(ncB,ncS); % theorycS(:,c) = tcS;
expcS = zeros(ncB,ncS); % expcS(:,c) = ecS;
obscS = zeros(ncB,ncS); % obscS(:,c) = cS;
alphaS = zeros(ncB,ncS); % alphaS(:,c) = aS;
betaS = zeros(ncB,ncS); % betaS(:,c) = bS;
alphaC = zeros(ncC,ncd); % alphaC(c,:) = aC;
betaC = zeros(ncC,ncd); % betaC(c,:) = bC;
coexC = zeros(ncC,ncd); % coexC(c,:) = cC;
fafb = zeros(ncC,2); % fafb(c,:) = [fa fb];
efpafpb = zeros(ncC,2); % efpafpb(c,:) = y' or [efpa efpb];
tfpafpb = zeros(ncC,2); % tfpafpb(c,:) = [fpa fpb];
varefp = zeros(ncC,2); % varfp(c,:) = [covar(1,1) covar(2,2)];
Kp = zeros(ncC,1); % Kp(c,1) = kp;
eKp = zeros(ncC,1); % eKp(c,1) = ekp;
uC = zeros(ncC,1); % uC(c) = uc;
vC = zeros(ncC,1); % vC(c) = vc;
residual = zeros(ncB,ncS); % residual(:,c) = (cS-ecS)./(stdSc(c)*sqrt(ncB*ncS));
chisq_TD = zeros(ncC,1); % chisq_TD(c,1) = ((cS-tcS)'*(cS-tcS))/(varSc(c)*ncB);
chisq_TE = zeros(ncC,1); % chisq_TE(c,1) = ((ecS-tcS)'*(ecS-tcS))/(varSc(c)*ncB);
chisq_ED = zeros(ncC,1); % chisq_ED(c,1) = ((cS-ecS)'*(cS-ecS))/(varSc(c)*ncB);
chisq_AD = zeros(ncC,1); % chisq_AD(c,1) = ((cS-aS)'*(cS-aS))/(varSc(c)*ncB);
chisq_BD = zeros(ncC,1); % chisq_BD(c,1) = ((cS-bS)'*(cS-bS))/(varSc(c)*ncB);
chisq_AB = zeros(ncC,1); % chisq_AB(c,1) = ((aS-bS)'*(aS-bS))/(varSc(c)*ncB);
chisq_AE = zeros(ncC,1); % chisq_AE(c,1) = ((aS-ecS)'*(aS-ecS))/(varSc(c)*ncB);
chisq_BE = zeros(ncC,1); % chisq_BE(c,1) = ((bS-ecS)'*(bS-ecS))/(varSc(c)*ncB);

switch boundary.config
    case '2critpts'
        b = boundary.parameterization;
        [u,directrixO,v,directrixA] = divide_bdy(x(1),x(2),C_bdy);
        % x(1) = the boundary parameter of critical point 1
        % x(2) = the boundary parameter of critical point 2
        % the La directrix is the section of the boundary from x(1) to x(2)
        % in increasing order and the Lo directrix is from x(1) to x(2) in
        % decreasing order.
        
        % the tangent line to a critical point is calculated from a
        % quadratic spline fit to an appropriate section of the piece-wise
        % linear boundary around the critical point.
        [cpt1,braket1] = b2bdypt(b,C_bdy,x(1));
        [cpt2,braket2] = b2bdypt(b,C_bdy,x(2));
        bracket1 = braket1{:};
        bracket2 = braket2{:};
        if isscalar(bracket1)
            if bracket1 == 1
                cpt1sectionb = b([nbC-1;bracket1;bracket1+1]);
            elseif bracket1 == nbC
                cpt1sectionb = b([bracket1-1;bracket1;2]);
            else
                cpt1sectionb = b([bracket1-1;bracket1;bracket1+1]);
            end
            cpt1section = b2bdypt(b,C_bdy,cpt1sectionb);
        elseif isvector(bracket1)
            if bracket1(1) == 1
                cpt1sectionb = b([nbC-1;bracket1(1);bracket1(2);bracket1(2)+1]);
            elseif bracket1(2) == nbC
                cpt1sectionb = b([bracket1(1)-1;bracket1(1);bracket1(2);2]);
            else
                cpt1sectionb = b([bracket1(1)-1;bracket1(1);bracket1(2);bracket1(2)+1]);
            end
            cpt1section = b2bdypt(b,C_bdy,cpt1sectionb);
        else
            error('bracket1 cannot be a matrix');
        end
        if isscalar(bracket2)
            if bracket2 == 1
                cpt2sectionb = b([nbC-1;bracket2;bracket2+1]);
            elseif bracket2 == nbC
                cpt2sectionb = b([bracket2-1;bracket2;2]);
            else
                cpt2sectionb = b([bracket2-1;bracket2;bracket2+1]);
            end
            cpt2section = b2bdypt(b,C_bdy,cpt2sectionb);
        elseif isvector(bracket2)
            if bracket2(1) == 1
                cpt2sectionb = b([nbC-1;bracket2(1);bracket2(2);bracket2(2)+1]);
            elseif bracket2(2) == nbC
                cpt2sectionb = b([bracket2(1)-1;bracket2(1);bracket2(2);2]);
            else
                cpt2sectionb = b([bracket2(1)-1;bracket2(1);bracket2(2);bracket2(2)+1]);
            end
            cpt2section = b2bdypt(b,C_bdy,cpt2sectionb);
        else
            error('bracket2 cannot be a matrix');
        end
%         if all(cpt1sectionb == cpt2sectionb)
%             mc1 = 
%         else
%             
%         end
        % fit each coordinate of the boundary section around critical point 1 with a parabola (quadratic)
        cpt1splinex = spap2(1,3,cpt1sectionb,cpt1section(:,1)); % knot sequence chosen for me 
        cpt1spliney = spap2(1,3,cpt1sectionb,cpt1section(:,2)); % knot sequence chosen for me
        % fit each coordinate of the boundary section around critical point 2 with a parabola (quadratic)
        cpt2splinex = spap2(1,3,cpt2sectionb,cpt2section(:,1)); % knot sequence chosen for me
        cpt2spliney = spap2(1,3,cpt2sectionb,cpt2section(:,2)); % knot sequence chosen for me
        % find slope of tangent lines
        warning off;
        mc1 = fnval(fnder(cpt1spliney),x(1))/fnval(fnder(cpt1splinex),x(1));
        mc2 = fnval(fnder(cpt2spliney),x(2))/fnval(fnder(cpt2splinex),x(2));
        warning on;
            
        % represent tangent lines in 2 pt form to find their intersection
        if isequal(mc1,Inf) | isequal(mc1,-Inf) % mc1 == Inf or -Inf
            lina = [cpt1;cpt1(1) 1];
        elseif mc1 > 0
            lina = [cpt1;1 mc1*(1-cpt1(1))+cpt1(2)];
        elseif isequal(mc1,0) % mc1 == 0
            lina = [cpt1;1 cpt1(2)];
        elseif mc1 < 0 
            lina = [cpt1;0 mc1*(0-cpt1(1))+cpt1(2)];
        elseif isnan(mc1)
            keyboard;
%             if fnval(fnder(spy),cp(1)) == Inf | fnval(fnder(spy),cp(1)) == -Inf
%                 lina = [cpt1;cpt1(1) 1];
%             end
%             if fnval(fnder(spy),cp(1)) == 0
%                 lina = [cpt1;1 cpt1(2)];
%             end 
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
            keyboard;
%             if fnval(fnder(spy),cp(2)) == Inf | fnval(fnder(spy),cp(2)) == -Inf
%                 linb = [cpt2;cpt2(1) 1];
%             end
%             if fnval(fnder(spy),cp(2)) == 0
%                 linb = [cpt2;1 cpt2(2)];
%             end 
        end

        % find intersection pt of tangent lines
        intersection = line_line_int(lina,linb);
        if isempty(intersection) % lines parallel
            slopes = repmat(mc1,ncC,1); % or mc2
        else 
            angles = pt2angle(C_coex-repmat(intersection,ncC,1));
            slopes = angle2slope(angles);
        end
        m_at_ipts = [mc1 mc2];
        clear lina linb mc1 mc2;
        invariant_pts = struct('cpt1',cpt1,'cpt2',cpt2);
        
        switch part_coeff.config
            case 'constant' % 2critpts
                % x = parameters, length(x) should = 3 = [cp1 cp2 Kp]
                if length(x) ~= 3
                    error('need 3 parameters, [cp1 cp2 Kp]');
                end
                % loop over phase coexistence points and calculate spectra at each point
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
                    int_pts = boundary_line_int(C_bdy,slopes(c),cC);
                    nipt = size(int_pts,1);

                    if isempty(int_pts) % no intersection
                        error('tie line does not intersect boundary');
                    elseif isequal(nipt,1) % tangent intersection
                        error('tie line intersects boundary at only 1 point');
                    elseif isequal(nipt,2) % 2 intersections
                        pt1 = int_pts(1,:);
                        pt2 = int_pts(2,:);
                    elseif nipt > 2 
                        error('tie line intersects boundary at > 2 points');
                    end

                    % determine which directrix pt1 and pt2 are located
                    u1 = bdypt2b(pt1,u,directrixO);
                    if isnan(u1)
                        v1 = bdypt2b(pt1,v,directrixA);
                        if isnan(v1)
                            error('pt1 not on either boundary section (Lo or La)');
                        else
                            aC = pt1;
                        end
                    else
                        bC = pt1;
                    end
                    u2 = bdypt2b(pt2,u,directrixO);
                    if isnan(u2)
                        v2 = bdypt2b(pt2,v,directrixA);
                        if isnan(v2)
                            error('pt2 not on either boundary section (Lo or La)');
                        else
                            aC = pt2;
                        end
                    else
                        bC = pt2;
                    end
                    uc = bdypt2b(bC,u,directrixO);
                    vc = bdypt2b(aC,v,directrixA);

%                     if exist('aC') == 0
%                         keyboard;
%                     end
                    % interpolate spectra at alpha point on the boundary
                    ab = bdypt2b(aC,b,C_bdy); % bi = bdypt2b(bdypt,b,bdy) 
                %     ap = bdypt2param(C_bdy,aC,'polygon',b);
                    aS = interp1(b,S_bdy',ab);
                    aS = aS';   
                    % normalize spectrum
                    aS = normalize_spectra([spectra.magnetic_fields aS],'d','area',1000);
                    aS = aS(:,2);
                    % interpolate spectra at beta point on the boundary
                    bb = bdypt2b(bC,b,C_bdy);
                    bS = interp1(b,S_bdy',bb);
                    bS = bS'; 
                    % normalize spectrum
                    bS = normalize_spectra([spectra.magnetic_fields bS],'d','area',1000);
                    bS = bS(:,2);

                    % kp defined as into beta (Lo) phase:
                    %   fpalpha = fa./(fa + kp.*fb);
                    %   fpbeta = (kp.*fb)./(fa + kp.*fb);

                    % calculate fraction of alpha phase and fraction of beta phase
                    if all(aC == bC) % only at critical point
                        fb = 1;
                        fa = 1;
                    else
                        fb = norm(cC-aC)/norm(aC-bC);
                        fa = 1-fb;
                    end

                    % calculate spectrum at coexistence point
                    kp = x(3);
                    fpb = (kp*fb)/(fa + kp*fb);
                    fpa = fa/(fa + kp*fb);
                    tcS = fpa*aS + fpb*bS;
                    
                    options = optimset('LargeScale','off','Display','off');
                    y = lsqlin([aS.*wfxnfp bS.*wfxnfp],cS.*wfxnfp,[1 0;0 1],[1-(1e-10);1-(1e-10)],[1 1],1,[],[],[],options);
%                     y = lsqlin([aS.*wfxnfp bS.*wfxnfp],cS.*wfxnfp,[],[],[1 1],1,[],[],[],options);
                    matrix = ([aS bS]'*[aS bS]);
                    covar = inv(matrix);
                    efpa = y(1);
                    efpb = y(2);
                    ecS = efpa*aS + efpb*bS;
                    ekp = (efpb*fa)/(efpa*fb);
                    
                    uC(c) = uc; 
                    vC(c) = vc; 
                    fafb(c,:) = [fa fb];
                    efpafpb(c,:) = [efpa efpb];
                    tfpafpb(c,:) = [fpa fpb];
                    varefp(c,:) = [covar(1,1) covar(2,2)];
                    Kp(c,1) = kp;
                    eKp(c,1) = ekp;
                    alphaC(c,:) = aC;
                    betaC(c,:) = bC;
                    coexC(c,:) = cC;
                    alphaS(:,c) = aS;
                    betaS(:,c) = bS;
                    obscS(:,c) = cS;
                    theorycS(:,c) = tcS;
                    expcS(:,c) = ecS;
                    chisq_TD(c,1) = ((cS-tcS)'*(cS-tcS))/(varSc(c)*ncB);
                    chisq_TE(c,1) = ((ecS-tcS)'*(ecS-tcS))/(varSc(c)*ncB);
                    chisq_ED(c,1) = ((cS-ecS)'*(cS-ecS))/(varSc(c)*ncB);
                    chisq_AD(c,1) = ((cS-aS)'*(cS-aS))/(varSc(c)*ncB);
                    chisq_BD(c,1) = ((cS-bS)'*(cS-bS))/(varSc(c)*ncB);
                    chisq_AB(c,1) = ((aS-bS)'*(aS-bS))/(varSc(c)*ncB);
                    residual(:,c) = ((cS-tcS).*wfxnS)./(stdSc(c)*sqrt(ncS*Nfit));
%                     chisquare(c,1) = ((((cS-ecS).*wfxnS)'*((cS-ecS).*wfxnS))/(varSc(c)*ncB*ncS));
%                     residual(:,c) = (((cS-ecS).*wfxnS)./(stdSc(c)*sqrt(ncB*ncS)));
                end % end of loop over phase coexistence points
            case 'function' % 2critpts
                % x = parameters, length(x) should = 3 = [cp1 cp2 c]
                if length(x) ~= 3
                    error('need 3 parameters, [cp1 cp2 c]');
                end
%               % loop over phase coexistence points and calculate spectra at each point
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
                    int_pts = boundary_line_int(C_bdy,slopes(c),cC);
                    nipt = size(int_pts,1);

                    if isempty(int_pts) % no intersection
                        error('tie line does not intersect boundary');
                    elseif isequal(nipt,1) % tangent intersection
                        error('tie line intersects boundary at only 1 point');
                    elseif isequal(nipt,2) % 2 intersections
                        pt1 = int_pts(1,:);
                        pt2 = int_pts(2,:);
                    elseif nipt > 2 
                        error('tie line intersects boundary at > 2 points');
                    end

                    % determine which directrix pt1 and pt2 are located
                    u1 = bdypt2b(pt1,u,directrixO);
                    if isnan(u1)
                        v1 = bdypt2b(pt1,v,directrixA);
                        if isnan(v1)
                            error('pt1 not on either boundary section (Lo or La)');
                        else
                            aC = pt1;
                        end
                    else
                        bC = pt1;
                    end
                    u2 = bdypt2b(pt2,u,directrixO);
                    if isnan(u2)
                        v2 = bdypt2b(pt2,v,directrixA);
                        if isnan(v2)
                            error('pt2 not on either boundary section (Lo or La)');
                        else
                            aC = pt2;
                        end
                    else
                        bC = pt2;
                    end
                    uc = bdypt2b(bC,u,directrixO);
                    vc = bdypt2b(aC,v,directrixA);

%                     if exist('aC') == 0
%                         keyboard;
%                     end
                    % interpolate spectra at alpha point on the boundary
                    ab = bdypt2b(aC,b,C_bdy); % bi = bdypt2b(bdypt,b,bdy) 
                %     ap = bdypt2param(C_bdy,aC,'polygon',b);
                    aS = interp1(b,S_bdy',ab);
                    aS = aS';   
                    % normalize spectrum
                    aS = normalize_spectra([spectra.magnetic_fields aS],'d','area',1000);
                    aS = aS(:,2);
                    % interpolate spectra at beta point on the boundary
                    bb = bdypt2b(bC,b,C_bdy);
                    bS = interp1(b,S_bdy',bb);
                    bS = bS'; 
                    % normalize spectrum
                    bS = normalize_spectra([spectra.magnetic_fields bS],'d','area',1000);
                    bS = bS(:,2);

                    % kp defined as into beta (Lo) phase:
                    %   fpalpha = fa./(fa + kp.*fb);
                    %   fpbeta = (kp.*fb)./(fa + kp.*fb);

                    % calculate fraction of alpha phase and fraction of beta phase
                    if all(aC == bC)
                        fb = 1;
                        fa = 1;
                    else
                        fb = norm(cC-aC)/norm(aC-bC);
                        fa = 1-fb;
                    end

                    % calculate spectrum at coexistence point
                    kp = Kp_fxn(x(3),uc,boundary);
                    fpb = (kp*fb)/(fa + kp*fb);
                    fpa = fa/(fa + kp*fb);
                    tcS = fpa*aS + fpb*bS;
                    
                    options = optimset('LargeScale','off','Display','off');
                    y = lsqlin([aS.*wfxnfp bS.*wfxnfp],cS.*wfxnfp,[1 0;0 1],[1-(1e-10);1-(1e-10)],[1 1],1,[],[],[],options);
%                     y = lsqlin([aS.*wfxnfp bS.*wfxnfp],cS.*wfxnfp,[],[],[1 1],1,[],[],[],options);
                    matrix = ([aS bS]'*[aS bS]);
                    covar = inv(matrix);
                    efpa = y(1);
                    efpb = y(2);
                    ecS = efpa*aS + efpb*bS;
                    ekp = (efpb*fa)/(efpa*fb);
                    
                    uC(c) = uc; 
                    vC(c) = vc; 
                    fafb(c,:) = [fa fb];
                    efpafpb(c,:) = [efpa efpb];
                    tfpafpb(c,:) = [fpa fpb];
                    varefp(c,:) = [covar(1,1) covar(2,2)];
                    Kp(c,1) = kp;
                    eKp(c,1) = ekp;
                    alphaC(c,:) = aC;
                    betaC(c,:) = bC;
                    coexC(c,:) = cC;
                    alphaS(:,c) = aS;
                    betaS(:,c) = bS;
                    obscS(:,c) = cS;
                    theorycS(:,c) = tcS;
                    expcS(:,c) = ecS;
                    chisq_TD(c,1) = (((cS-tcS).*wfxnS)'*((cS-tcS).*wfxnS))/(varSc(c)*Nfit); %((cS-tcS)'*(cS-tcS))/(varSc(c)*ncB);
                    chisq_TE(c,1) = (((ecS-tcS).*wfxnS)'*((ecS-tcS).*wfxnS))/(varSc(c)*Nfit); %((ecS-tcS)'*(ecS-tcS))/(varSc(c)*ncB);
                    chisq_ED(c,1) = (((cS-ecS).*wfxnS)'*((cS-ecS).*wfxnS))/(varSc(c)*Nfit); %((cS-ecS)'*(cS-ecS))/(varSc(c)*ncB);
                    chisq_AD(c,1) = (((cS-aS).*wfxnS)'*((cS-aS).*wfxnS))/(varSc(c)*Nfit); %((cS-aS)'*(cS-aS))/(varSc(c)*ncB);
                    chisq_BD(c,1) = (((cS-bS).*wfxnS)'*((cS-bS).*wfxnS))/(varSc(c)*Nfit); %((cS-bS)'*(cS-bS))/(varSc(c)*ncB);
                    chisq_AB(c,1) = (((aS-bS).*wfxnS)'*((aS-bS).*wfxnS))/(varSc(c)*Nfit); %((aS-bS)'*(aS-bS))/(varSc(c)*ncB);
                    chisq_AE(c,1) = (((aS-ecS).*wfxnS)'*((aS-ecS).*wfxnS))/(varSc(c)*Nfit);
                    
%                     chisq_TD(c,1) = ((cS-tcS)'*(cS-tcS))/(varSc(c)*ncB);
%                     chisq_TE(c,1) = ((ecS-tcS)'*(ecS-tcS))/(varSc(c)*ncB);
%                     chisq_ED(c,1) = ((cS-ecS)'*(cS-ecS))/(varSc(c)*ncB);
%                     chisq_AD(c,1) = ((cS-aS)'*(cS-aS))/(varSc(c)*ncB);
%                     chisq_BD(c,1) = ((cS-bS)'*(cS-bS))/(varSc(c)*ncB);
%                     chisq_AB(c,1) = ((aS-bS)'*(aS-bS))/(varSc(c)*ncB);
                    residual(:,c) = ((cS-ecS).*wfxnS)./(stdSc(c)*sqrt(ncS*Nfit));
%                     chisquare(c,1) = ((((cS-ecS).*wfxnS)'*((cS-ecS).*wfxnS))/(varSc(c)*ncB*ncS));
%                     residual(:,c) = (((cS-ecS).*wfxnS)./(stdSc(c)*sqrt(ncB*ncS)));
                end % end of loop over phase coexistence points
            otherwise
                error('invalid partition coefficient configuration');
        end    
    case '1critpt/1endtl'
        b = boundary.parameterization;
        cp = x(1);
        e1 = x(2);
        e2 = x(3);
        if e1 >= e2
            % divide boundary into two sections (same procedure
            % when going from tangent to ruled in the 2
            % critical point case above).
            [s1,section1,s2,section2] = divide_bdy(e1,e2,C_bdy);
%             [u,directrixO,v,directrixA] = tangent2ruled(e1,e2,C_bdy);
            if cp < e1 & cp > e2 % configuration D in notes
                Cbdy = section1;
                bdy = [Cbdy;Cbdy(1,:)];
                cpt = b2bdypt(b,C_bdy,x(1));
                cps1 = bdypt2b(cpt,s1,section1);
                [cpt,braket] = b2bdypt(s1,section1,cps1);
                bracket = braket{:};
                if isscalar(bracket)
                    i = bracket;
                    directrixO = flipud(section1(1:i,:));
                    directrixA = section1(i:end,:);
                elseif isvector(bracket)
                    il = bracket(1);
                    iu = bracket(2);
                    directrixO = [cpt;flipud(section1(1:il,:))];
                    directrixA = [cpt;section1(iu:end,:)];
                else
                    error('bracket cannot be a matrix');
                end
            elseif cp > e1 | cp < e2 % configuration B in notes
                Cbdy = section2;
                bdy = [Cbdy;Cbdy(1,:)];
                cpt = b2bdypt(b,C_bdy,x(1));
                cps2 = bdypt2b(cpt,s2,section2);
                [cpt,braket] = b2bdypt(s2,section2,cps2);
                bracket = braket{:};
                if isscalar(bracket)
                    i = bracket;
                    directrixA = flipud(section2(1:i,:));
                    directrixO = section2(i:end,:);
                elseif isvector(bracket)
                    il = bracket(1);
                    iu = bracket(2);
                    directrixA = [cpt;flipud(section2(1:il,:))];
                    directrixO = [cpt;section2(iu:end,:)];
                else
                    error('bracket cannot be a matrix');
                end
            else % cp == e1 | cp == e2
                error('cp cannot equal e1 or e2');
            end
            [Olength,u] = chord_length(directrixO);
            [Alength,v] = chord_length(directrixA);
            [IN ON] = inpolygon(C_coex(:,1),C_coex(:,2),bdy(:,1),bdy(:,2));
            Ccoex = C_coex(IN,:);
            Scoex = S_coex(:,IN);
        else % e1 < e2
            % divide boundary into two sections (same procedure
            % when going from tangent to ruled in the 2
            % critical point case above).
            [s1,section1,s2,section2] = divide_bdy(x(2),x(3),C_bdy);% x(2),x(3) = e1,e2
%             [u,directrixO,v,directrixA] = tangent2ruled(e1,e2,C_bdy);
            if cp > e1 & cp < e2 % configuration C in notes
                Cbdy = section2;
                bdy = [Cbdy;Cbdy(1,:)];
                cpt = b2bdypt(b,C_bdy,x(1));
                cps2 = bdypt2b(cpt,s2,section2);
                [cpt,braket] = b2bdypt(s2,section2,cps2);
                bracket = braket{:};
                if isscalar(bracket)
                    i = bracket;
                    directrixO = flipud(section2(1:i,:));
                    directrixA = section2(i:end,:);
                elseif isvector(bracket)
                    il = bracket(1);
                    iu = bracket(2);
                    directrixO = [cpt;flipud(section2(1:il,:))];
                    directrixA = [cpt;section2(iu:end,:)];
                else
                    error('bracket cannot be a matrix');
                end
            elseif cp < e1 | cp > e2 % configuration A in notes
                Cbdy = section1;
                bdy = [Cbdy;Cbdy(1,:)];
                cpt = b2bdypt(b,C_bdy,x(1));
                cps1 = bdypt2b(cpt,s1,section1);
                [cpt,braket] = b2bdypt(s1,section1,cps1);
                bracket = braket{:};
                if isscalar(bracket)
                    i = bracket;
                    directrixA = flipud(section1(1:i,:));
                    directrixO = section1(i:end,:);
                elseif isvector(bracket)
                    il = bracket(1);
                    iu = bracket(2);
                    directrixA = [cpt;flipud(section1(1:il,:))];
                    directrixO = [cpt;section1(iu:end,:)];
                else
                    error('bracket cannot be a matrix');
                end
            else % cp == e1 | cp == e2
                error('cp cannot equal e1 or e2');
            end
            [Olength,u] = chord_length(directrixO);
            [Alength,v] = chord_length(directrixA);
            [IN ON] = inpolygon(C_coex(:,1),C_coex(:,2),bdy(:,1),bdy(:,2));
            Ccoex = C_coex(IN,:);
            Scoex = S_coex(:,IN);
        end
            
        cin = find(IN);
        nin = length(cin);
        cout = find(~IN);
        nout = length(cout);
        
        % the tangent line to a critical point is calculated from a
        % quadratic spline fit to an appropriate section of the piece-wise
        % linear boundary around the critical point.
        [cpt,braket] = b2bdypt(b,C_bdy,x(1));
        bracket = braket{:};
        if isscalar(bracket)
            if bracket == 1
                cptsectionb = b([nbC-1;bracket;bracket+1]);
            elseif bracket == nbC
                cptsectionb = b([bracket-1;bracket;2]);
            else
                cptsectionb = b([bracket-1;bracket;bracket+1]);
            end
            cptsection = b2bdypt(b,C_bdy,cptsectionb);
        elseif isvector(bracket)
            if bracket(1) == 1
                cptsectionb = b([nbC-1;bracket(1);bracket(2);bracket(2)+1]);
            elseif bracket(2) == nbC
                cptsectionb = b([bracket(1)-1;bracket(1);bracket(2);2]);
            else
                cptsectionb = b([bracket(1)-1;bracket(1);bracket(2);bracket(2)+1]);
            end
            cptsection = b2bdypt(b,C_bdy,cptsectionb);
        else
            error('bracket1 cannot be a matrix');
        end
        % fit each coordinate of the boundary section around critical point 1 with a parabola (quadratic)
        cptsplinex = spap2(1,3,cptsectionb,cptsection(:,1)); % knot sequence chosen for me 
        cptspliney = spap2(1,3,cptsectionb,cptsection(:,2)); % knot sequence chosen for me
        % find slope of tangent lines
        warning off;
        mc = fnval(fnder(cptspliney),x(1))/fnval(fnder(cptsplinex),x(1));
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
        elseif isnan(mc)
            keyboard;
%             if fnval(fnder(spy),cp(1)) == Inf | fnval(fnder(spy),cp(1)) == -Inf
%                 lina = [cpt1;cpt1(1) 1];
%             end
%             if fnval(fnder(spy),cp(1)) == 0
%                 lina = [cpt1;1 cpt1(2)];
%             end 
        end

        ept1 = b2bdypt(b,C_bdy,e1);
        ept2 = b2bdypt(b,C_bdy,e2);
        % find intersection pt of tangent line and end tieline
        intersection = line_line_int(lin,[ept1;ept2]);
        if isempty(intersection) % lines parallel
            slopes = repmat(mc,ncC,1); 
        else 
            angles = pt2angle(C_coex-repmat(intersection,ncC,1));
            slopes = angle2slope(angles);
        end
        m_at_ipts = [mc (ept1(2)-ept2(2))/(ept1(1)-ept2(1))];
        clear lin mc;
        invariant_pts = struct('cpt',cpt,'ept1',ept1,'ept2',ept2);
        
        switch part_coeff.config
            case 'constant' % '1critpt/1endtl'
                % x = parameters, length(x) should = 4 = [cp ep1 ep2 Kp]
                if length(x) ~= 4
                    error('need 4 parameters, [cp ep1 ep2 Kp]');
                end 
                % loop over phase coexistence points and calculate spectra at each point
                for c = cin' % 1:ncC
                    cS = S_coex(:,c);
                    cC = C_coex(c,:);

                    % the intersection of the tie line drawn through the intersection
                    % point of the critical point tangent line(s) and/or end tielines and the
                    % coexistence point with the boundary gives both the composition
                    % point of the alpha phase (C_alpha, or Ld phase) and the
                    % composition of the beta phase (C_beta, or Lo phase).  the
                    % points returned from boundary_line_int are that
                    % C_alpha=int_pts(1,:) and C_beta=int_pts(2,:)
                    int_pts = boundary_line_int(C_bdy,slopes(c),cC);
                    nipt = size(int_pts,1);

                    if isempty(int_pts) % no intersection
                        error('tie line does not intersect boundary');
                    elseif isequal(nipt,1) % tangent intersection
                        error('tie line intersects boundary at only 1 point');
                    elseif isequal(nipt,2) % 2 intersections
                        pt1 = int_pts(1,:);
                        pt2 = int_pts(2,:);
                    elseif nipt > 2 
                        error('tie line intersects boundary at > 2 points');
                    end
                    
                    % determine which directrix pt1 and pt2 are located
                    u1 = bdypt2b(pt1,u,directrixO);
                    if isnan(u1)
                        v1 = bdypt2b(pt1,v,directrixA);
                        if isnan(v1)
                            error('pt1 not on either boundary section (Lo or La)');
                        else
                            aC = pt1;
                        end
                    else
                        bC = pt1;
                    end
                    u2 = bdypt2b(pt2,u,directrixO);
                    if isnan(u2)
                        v2 = bdypt2b(pt2,v,directrixA);
                        if isnan(v2)
                            error('pt2 not on either boundary section (Lo or La)');
                        else
                            aC = pt2;
                        end
                    else
                        bC = pt2;
                    end
                    uc = bdypt2b(bC,u,directrixO);
                    vc = bdypt2b(aC,v,directrixA);

%                     if exist('aC') == 0
%                         keyboard;
%                     end
                    % interpolate spectra at alpha point on the boundary
                    ab = bdypt2b(aC,b,C_bdy); % bi = bdypt2b(bdypt,b,bdy) 
                %     ap = bdypt2param(C_bdy,aC,'polygon',b);
                    aS = interp1(b,S_bdy',ab);
                    aS = aS';   
                    % normalize spectrum
                    aS = normalize_spectra([spectra.magnetic_fields aS],'d','area',1000);
                    aS = aS(:,2);
                    % interpolate spectra at beta point on the boundary
                    bb = bdypt2b(bC,b,C_bdy);
                    bS = interp1(b,S_bdy',bb);
                    bS = bS'; 
                    % normalize spectrum
                    bS = normalize_spectra([spectra.magnetic_fields bS],'d','area',1000);
                    bS = bS(:,2);

                    % kp defined as into beta (Lo) phase:
                    %   fpalpha = fa./(fa + kp.*fb);
                    %   fpbeta = (kp.*fb)./(fa + kp.*fb);

                    % calculate fraction of alpha phase and fraction of beta phase
                    if all(aC == bC)
                        fb = 1;
                        fa = 1;
                    else
                        fb = norm(cC-aC)/norm(aC-bC);
                        fa = 1-fb;
%                         fo = det([coexpt(3) coexpt(1);aC(3) aC(1)])/det([aC(1) aC(3);oC(1) oC(3)]);
%                         fa = 1-fo;
                    end

                    % calculate spectrum at coexistence point
                    kp = x(4);
                    fpb = (kp*fb)/(fa + kp*fb);
                    fpa = fa/(fa + kp*fb);
                    tcS = fpa*aS + fpb*bS;
                    
                    options = optimset('LargeScale','off','Display','off');
                    y = lsqlin([aS.*wfxnfp bS.*wfxnfp],cS.*wfxnfp,[1 0;0 1],[1-(1e-10);1-(1e-10)],[1 1],1,[],[],[],options);
%                     y = lsqlin([aS.*wfxnfp bS.*wfxnfp],cS.*wfxnfp,[],[],[1 1],1,[],[],[],options);
                    matrix = ([aS bS]'*[aS bS]);
                    covar = inv(matrix);
                    efpa = y(1);
                    efpb = y(2);
                    ecS = efpa*aS + efpb*bS;
                    ekp = (efpb*fa)/(efpa*fb);
                    
                    uC(c) = uc; 
                    vC(c) = vc; 
                    fafb(c,:) = [fa fb];
                    efpafpb(c,:) = [efpa efpb];
                    tfpafpb(c,:) = [fpa fpb];
                    varefp(c,:) = [covar(1,1) covar(2,2)];
                    Kp(c,1) = kp;
                    eKp(c,1) = ekp;
                    alphaC(c,:) = aC;
                    betaC(c,:) = bC;
                    coexC(c,:) = cC;
                    alphaS(:,c) = aS;
                    betaS(:,c) = bS;
                    obscS(:,c) = cS;
                    theorycS(:,c) = tcS;
                    expcS(:,c) = ecS;
                    chisq_TD(c,1) = ((cS-tcS)'*(cS-tcS))/(varSc(c)*ncB);
                    chisq_TE(c,1) = ((ecS-tcS)'*(ecS-tcS))/(varSc(c)*ncB);
                    chisq_ED(c,1) = ((cS-ecS)'*(cS-ecS))/(varSc(c)*ncB);
                    chisq_AD(c,1) = ((cS-aS)'*(cS-aS))/(varSc(c)*ncB);
                    chisq_BD(c,1) = ((cS-bS)'*(cS-bS))/(varSc(c)*ncB);
                    chisq_AB(c,1) = ((aS-bS)'*(aS-bS))/(varSc(c)*ncB);
                    residual(:,c) = ((cS-tcS).*wfxnS)./(stdSc(c)*sqrt(ncS*Nfit));
%                     chisquare(c,1) = ((((cS-ecS).*wfxnS)'*((cS-ecS).*wfxnS))/(varSc(c)*ncB*ncS));
%                     residual(:,c) =
                end % end of loop over phase coexistence points
                for c = cout'
                    uC(c) = NaN; 
                    vC(c) = NaN; 
                    fafb(c,:) = [NaN NaN];
                    efpafpb(c,:) = [NaN NaN];
                    tfpafpb(c,:) = [NaN NaN];
                    varefp(c,:) = [NaN NaN];
                    Kp(c,1) = NaN;
                    eKp(c,1) = NaN;
                    alphaC(c,:) = [NaN NaN];
                    betaC(c,:) = [NaN NaN];
                    coexC(c,:) = [NaN NaN];
                    alphaS(:,c) = repmat(NaN,ncB,1);
                    betaS(:,c) = repmat(NaN,ncB,1);
                    obscS(:,c) = repmat(NaN,ncB,1);
                    theorycS(:,c) = repmat(NaN,ncB,1);
                    expcS(:,c) = repmat(NaN,ncB,1);
                    chisq_TD(c,1) = NaN;
                    chisq_TE(c,1) = NaN;
                    chisq_ED(c,1) = NaN;
                    chisq_AD(c,1) = NaN;
                    chisq_BD(c,1) = NaN;
                    chisq_AB(c,1) = NaN;
                    residual(:,c) = 0;
                end
            case 'function' % '1critpt/1endtl'
                % x = parameters, length(x) should = 5 = [cp ep1 ep2 c1 c2]
                if length(x) ~= 5
                    error('need 5 parameters, [cp ep1 ep2 c1 c2]');
                end
                % loop over phase coexistence points and calculate spectra at each point
                for c = cin' % 1:ncC
                    cS = S_coex(:,c);
                    cC = C_coex(c,:);

                    % the intersection of the tie line drawn through the intersection
                    % point of the critical point tangent line(s) and/or end tielines and the
                    % coexistence point with the boundary gives both the composition
                    % point of the alpha phase (C_alpha, or Ld phase) and the
                    % composition of the beta phase (C_beta, or Lo phase).  the
                    % points returned from boundary_line_int are that
                    % C_alpha=int_pts(1,:) and C_beta=int_pts(2,:)
                    int_pts = boundary_line_int(C_bdy,slopes(c),cC);
                    nipt = size(int_pts,1);

                    if isempty(int_pts) % no intersection
                        error('tie line does not intersect boundary');
                    elseif isequal(nipt,1) % tangent intersection
                        error('tie line intersects boundary at only 1 point');
                    elseif isequal(nipt,2) % 2 intersections
                        pt1 = int_pts(1,:);
                        pt2 = int_pts(2,:);
                    elseif nipt > 2 
                        error('tie line intersects boundary at > 2 points');
                    end

                    % determine which directrix pt1 and pt2 are located
                    u1 = bdypt2b(pt1,u,directrixO);
                    if isnan(u1)
                        v1 = bdypt2b(pt1,v,directrixA);
                        if isnan(v1)
                            error('pt1 not on either boundary section (Lo or La)');
                        else
                            aC = pt1;
                        end
                    else
                        bC = pt1;
                    end
                    u2 = bdypt2b(pt2,u,directrixO);
                    if isnan(u2)
                        v2 = bdypt2b(pt2,v,directrixA);
                        if isnan(v2)
                            error('pt2 not on either boundary section (Lo or La)');
                        else
                            aC = pt2;
                        end
                    else
                        bC = pt2;
                    end
                    uc = bdypt2b(bC,u,directrixO);
                    vc = bdypt2b(aC,v,directrixA);

%                     if exist('aC') == 0
%                         keyboard;
%                     end
                    % interpolate spectra at alpha point on the boundary
                    ab = bdypt2b(aC,b,C_bdy); % bi = bdypt2b(bdypt,b,bdy) 
                %     ap = bdypt2param(C_bdy,aC,'polygon',b);
                    aS = interp1(b,S_bdy',ab);
                    aS = aS';   
                    % normalize spectrum
                    aS = normalize_spectra([spectra.magnetic_fields aS],'d','area',1000);
                    aS = aS(:,2);
                    % interpolate spectra at beta point on the boundary
                    bb = bdypt2b(bC,b,C_bdy);
                    bS = interp1(b,S_bdy',bb);
                    bS = bS'; 
                    % normalize spectrum
                    bS = normalize_spectra([spectra.magnetic_fields bS],'d','area',1000);
                    bS = bS(:,2);

                    % kp defined as into beta (Lo) phase:
                    %   fpalpha = fa./(fa + kp.*fb);
                    %   fpbeta = (kp.*fb)./(fa + kp.*fb);

                    % calculate fraction of alpha phase and fraction of beta phase
                    if all(aC == bC) % only at critical point
                        fb = 1;
                        fa = 1;
                    else
                        fb = norm(cC-aC)/norm(aC-bC);
                        fa = 1-fb;
                    end

                    % calculate spectrum at coexistence point
                    kp = Kp_fxn(x(4:5),uc,boundary);
                    fpb = (kp*fb)/(fa + kp*fb);
                    fpa = fa/(fa + kp*fb);
                    tcS = fpa*aS + fpb*bS;
                    
                    options = optimset('LargeScale','off','Display','off');
                    y = lsqlin([aS.*wfxnfp bS.*wfxnfp],cS.*wfxnfp,[1 0;0 1],[1-(1e-10);1-(1e-10)],[1 1],1,[],[],[],options);
%                     y = lsqlin([aS.*wfxnfp bS.*wfxnfp],cS.*wfxnfp,[],[],[1 1],1,[],[],[],options);
                    matrix = ([aS bS]'*[aS bS]);
                    covar = inv(matrix);
                    efpa = y(1);
                    efpb = y(2);
                    ecS = efpa*aS + efpb*bS;
                    ekp = (efpb*fa)/(efpa*fb);
                    
                    uC(c) = uc; 
                    vC(c) = vc; 
                    fafb(c,:) = [fa fb];
                    efpafpb(c,:) = [efpa efpb];
                    tfpafpb(c,:) = [fpa fpb];
                    varefp(c,:) = [covar(1,1) covar(2,2)];
                    Kp(c,1) = kp;
                    eKp(c,1) = ekp;
                    alphaC(c,:) = aC;
                    betaC(c,:) = bC;
                    coexC(c,:) = cC;
                    alphaS(:,c) = aS;
                    betaS(:,c) = bS;
                    obscS(:,c) = cS;
                    theorycS(:,c) = tcS;
                    expcS(:,c) = ecS;
                    chisq_TD(c,1) = (((cS-tcS).*wfxnS)'*((cS-tcS).*wfxnS))/(varSc(c)*Nfit); %((cS-tcS)'*(cS-tcS))/(varSc(c)*ncB);
                    chisq_TE(c,1) = (((ecS-tcS).*wfxnS)'*((ecS-tcS).*wfxnS))/(varSc(c)*Nfit); %((ecS-tcS)'*(ecS-tcS))/(varSc(c)*ncB);
                    chisq_ED(c,1) = (((cS-ecS).*wfxnS)'*((cS-ecS).*wfxnS))/(varSc(c)*Nfit); %((cS-ecS)'*(cS-ecS))/(varSc(c)*ncB);
                    chisq_AD(c,1) = (((cS-aS).*wfxnS)'*((cS-aS).*wfxnS))/(varSc(c)*Nfit); %((cS-aS)'*(cS-aS))/(varSc(c)*ncB);
                    chisq_BD(c,1) = (((cS-bS).*wfxnS)'*((cS-bS).*wfxnS))/(varSc(c)*Nfit); %((cS-bS)'*(cS-bS))/(varSc(c)*ncB);
                    chisq_AB(c,1) = (((aS-bS).*wfxnS)'*((aS-bS).*wfxnS))/(varSc(c)*Nfit); %((aS-bS)'*(aS-bS))/(varSc(c)*ncB);
                    
%                     chisq_TD(c,1) = ((cS-tcS)'*(cS-tcS))/(varSc(c)*ncB);
%                     chisq_TE(c,1) = ((ecS-tcS)'*(ecS-tcS))/(varSc(c)*ncB);
%                     chisq_ED(c,1) = ((cS-ecS)'*(cS-ecS))/(varSc(c)*ncB);
%                     chisq_AD(c,1) = ((cS-aS)'*(cS-aS))/(varSc(c)*ncB);
%                     chisq_BD(c,1) = ((cS-bS)'*(cS-bS))/(varSc(c)*ncB);
%                     chisq_AB(c,1) = ((aS-bS)'*(aS-bS))/(varSc(c)*ncB);
                    residual(:,c) = ((cS-ecS).*wfxnS)./(stdSc(c)*sqrt(ncS*Nfit));
%                     chisquare(c,1) = ((((cS-ecS).*wfxnS)'*((cS-ecS).*wfxnS))/(varSc(c)*ncB*ncS));
%                     residual(:,c) =
%                     (((cS-ecS).*wfxnS)./(stdSc(c)*sqrt(ncB*ncS)));
                end % end of loop over phase coexistence points
                for c = cout'
                    uC(c) = NaN; 
                    vC(c) = NaN; 
                    fafb(c,:) = [NaN NaN];
                    efpafpb(c,:) = [NaN NaN];
                    tfpafpb(c,:) = [NaN NaN];
                    varefp(c,:) = [NaN NaN];
                    Kp(c,1) = NaN;
                    eKp(c,1) = NaN;
                    alphaC(c,:) = [NaN NaN];
                    betaC(c,:) = [NaN NaN];
                    coexC(c,:) = [NaN NaN];
                    alphaS(:,c) = repmat(NaN,ncB,1);
                    betaS(:,c) = repmat(NaN,ncB,1);
                    obscS(:,c) = repmat(NaN,ncB,1);
                    theorycS(:,c) = repmat(NaN,ncB,1);
                    expcS(:,c) = repmat(NaN,ncB,1);
                    chisq_TD(c,1) = NaN;
                    chisq_TE(c,1) = NaN;
                    chisq_ED(c,1) = NaN;
                    chisq_AD(c,1) = NaN;
                    chisq_BD(c,1) = NaN;
                    chisq_AB(c,1) = NaN;
                    residual(:,c) = 0;
                end
            otherwise
                error('invalid partition coefficient configuration');
        end
    case '2endtls'
        switch part_coeff.config
            case 'constant'
                % x = parameters, length(x) should = 3 = [cp1 cp2 Kp]
                if length(x) ~= 3
                    error('need 3 parameters, [cp1 cp2 Kp]');
                end
            case 'function'
                % x = parameters, length(x) should = 3 = [cp1 cp2 c]
                if length(x) ~= 3
                    error('need 3 parameters, [cp1 cp2 c]');
                end
            otherwise
                error('invalid partition coefficient configuration');
        end
        eta = et{1};
        etb = et{2};
        % 2 end tielines
        intersection = line_line_int(eta,etb);
        if isempty(intersection) % lines parallel
            m = diff(eta); % or etb, slope should be the same
            m = m(2)/m(1);
            slopes = repmat(m,ncC,1); 
        else 
            angles = pt2angle(C_coex-repmat(intersection,ntc,1));
            slopes = angle2slope(angles);
        end
    otherwise
        error('invalid boundary configuration');
end

% td = chisq_TD;
% te = chisq_TE;
% ed = chisq_ED;
% X2area_triDET = (0.25).*sqrt(2*(td.*te + td.*ed + te.*ed) - (td.^2 + te.^2 + ed.^2));
% ad = chisq_AD;
% bd = chisq_BD;
% ab = chisq_AB;
% X2area_triDAB = (0.25).*sqrt(2*(ad.*bd + ad.*ab + bd.*ab) - (ad.^2 + bd.^2 + ab.^2)); 
X2_TD = nanmean(chisq_TD);
X2_TE = nanmean(chisq_TE);
X2_ED = nanmean(chisq_ED);
% Area_DET = nanmean(X2area_triDET);
X2_AD = nanmean(chisq_AD);
X2_BD = nanmean(chisq_BD);
X2_AB = nanmean(chisq_AB);
% Area_DAB = nanmean(X2area_triDAB);
projT = diag((obscS-theorycS)'*(obscS-theorycS));
% avg_projT = nanmean(projT);
projA = diag((obscS-theorycS)'*(alphaS-theorycS));
% avg_projA = nanmean(projA);
projB = diag((obscS-theorycS)'*(betaS-theorycS));
% avg_projB = nanmean(projB);
chi2fpa = ((efpafpb(:,1)-tfpafpb(:,1)).^2)./varefp(:,1);
chi2fpb = ((efpafpb(:,2)-tfpafpb(:,2)).^2)./varefp(:,2);
chi2fp = chi2fpa+chi2fpb;
Chi2fp = nanmean(chi2fp);
normfp = norm(efpafpb(~isnan(efpafpb(:)))-tfpafpb(~isnan(tfpafpb(:))));
normfp2 = sqrt(sum((efpafpb-tfpafpb).^2,2));
Normfp2 = nanmean(normfp2);

% output depends on the number of output arguments
if nargout > 1
    % output fit result data structure  
%     varargout{1} = X2_TD;
%     varargout{1} = X2_TD*normfp;
%     varargout{1} = X2_ED*normfp;
%     varargout{1} = X2_ED.*Normfp2;
    fit_measure = X2_TD;
    
    varargout{1} = fit_measure;
    
    configs = struct('bdyconfig',boundary.config,'tlconfig',{'tangent'},'Kpconfig',part_coeff.config);
    fit_data = struct('x',x','fit_measure',fit_measure,'chisq_TD',chisq_TD,'chisq_TE',chisq_TE,'chisq_ED',chisq_ED,...
        'chisq_AD',chisq_AD,'chisq_BD',chisq_BD,'chisq_AB',chisq_AB,'Nfit',Nfit,...
        'X2_TD',X2_TD,'X2_TE',X2_TE,'X2_ED',X2_ED,'X2_AD',X2_AD,'X2_BD',X2_BD,'X2_AB',X2_AB,...
        'projA',projA,'projB',projB,'projT',projT,'Chi2fp',Chi2fp,'chi2fp',chi2fp,'chi2fpa',chi2fpa,'chi2fpb',chi2fpb);
    comps_data = struct('invariant_pts',invariant_pts,'m_at_ipts',m_at_ipts,'aC',alphaC,'bC',betaC,'cC',coexC,'bdy',C_bdy);
    spectra_data = struct('magnetic_fields',spectra.magnetic_fields,...
        'aS',alphaS,'bS',betaS,'exp_cS',expcS,'obs_cS',obscS,'theory_cS',theorycS,...
        'cSbline_var',spectra.cSbline_var','cSbline_std',spectra.cSbline_std','bSC_std',spectra.bSC_std,'bSC_var',spectra.bSC_var,...
        'weight_fxn_S',spectra.weight_fxn_S,'weight_fxn_fp',spectra.weight_fxn_fp);
    tieline_data = struct('efpafpb',efpafpb,'tfpafpb',tfpafpb,'varefp',varefp,'normfp',normfp,'normfp2',normfp2,'Normfp2',Normfp2,...
        'fafb',fafb,'Kp',Kp,'eKp',eKp,'uC',uC,'vC',vC,'directrixO',directrixO,'directrixA',directrixA);
    varargout{2} = struct('configs',configs,'fit',fit_data,'comps',comps_data,'spectra',spectra_data,'tieline',tieline_data);
     
%     varargout{2} = struct('bdyconfig',boundary.config,'tlconfig',{'tangent'},'Kpconfig',part_coeff.config,'x',x',...
%         'chisq_TD',chisq_TD,'chisq_TE',chisq_TE,'chisq_ED',chisq_ED,'X2area_triDET',X2area_triDET,'Area_DET',Area_DET,...
%         'chisq_AD',chisq_AD,'chisq_BD',chisq_BD,'chisq_AB',chisq_AB,'X2area_triDAB',X2area_triDAB,'Area_DAB',Area_DAB,...
%         'projA',projA,'avg_projA',avg_projA,'projB',projB,'avg_projB',avg_projB,'projT',projT,'avg_projT',avg_projT,...
%         'Chi2fp',Chi2fp,'chi2fpa',chi2fpa,'chi2fpb',chi2fpb,'normfp',normfp,'normfp2',normfp2,'Normfp2',Normfp2,...
%         'invariant_pts',invariant_pts,'m_at_ipts',m_at_ipts,'aC',alphaC,'bC',betaC,'cC',coexC,...
%         'magnetic_fields',spectra.magnetic_fields,'aS',alphaS,'bS',betaS,'exp_cS',expcS,'obs_cS',obscS,'theory_cS',theorycS,...
%         'cSbline_var',spectra.cSbline_var','cSbline_std',spectra.cSbline_std','bSC_std',spectra.bSC_std,'bSC_var',spectra.bSC_var,...
%         'weight_fxn_S',spectra.weight_fxn_S,'weight_fxn_fp',spectra.weight_fxn_fp,...
%         'efpafpb',efpafpb,'tfpafpb',tfpafpb,'varefp',varefp,'fafb',fafb,'Kp',Kp,'eKp',eKp,'uC',uC,'vC',vC,'bdy',C_bdy);
else
    % output within fitting function depends on fit method
    switch fit_method
        case 'lsqnonlin'
%             varargout{1} = reshape(residual,ncB*ncS,1); % spectrum chi2 residual, or spectrum chi2 * fp chi2
%             varargout{1} = reshape(residual*sqrt(normfp),ncB*ncS,1);
            temp = reshape(residual*sqrt(normfp),ncB*ncS,1);
            varargout{1} = temp(~isnan(temp));
        case 'levmar'
        case 'simplex'
%             varargout{1} = X2_TD*X2_ED*(X2_TE/X2_AB);
%             varargout{1} = ((X2_TD*X2_ED*X2_TE)^1/3)/X2_AB;
%             varargout{1} = (X2_TD*X2_ED*X2_TE)^1/3;
%             varargout{1} = X2_TD*(1+normfp);
%             varargout{1} = X2_TD*normfp;
            varargout{1} = X2_TD;
%             varargout{1} = normfp;
%             varargout{1} = X2_ED*normfp;
%             varargout{1} = X2_ED.*Normfp2;
        otherwise
            error('invalid fit method');
    end
end

return