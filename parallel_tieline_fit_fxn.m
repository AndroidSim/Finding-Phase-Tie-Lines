function varargout = parallel_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method)
% varargout =
% parallel_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method)

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
regcS = zeros(ncB,ncS);
theorycS = zeros(ncB,ncS); % theorycS(:,c) = tcS;
expcS = zeros(ncB,ncS); % expcS(:,c) = ecS;
obscS = zeros(ncB,ncS); % obscS(:,c) = cS;
alphaS = zeros(ncB,ncS); % alphaS(:,c) = aS;
betaS = zeros(ncB,ncS); % betaS(:,c) = bS;
alphaC = zeros(ncC,ncd); % alphaC(c,:) = aC;
betaC = zeros(ncC,ncd); % betaC(c,:) = bC;
coexC = zeros(ncC,ncd); % coexC(c,:) = cC;
fafb = zeros(ncC,2); % fafb(c,:) = [fa fb];
zfpafpb = zeros(ncC,2); % zfpafpb(c,:) = [1-x_lambda x_lambda];
efpafpb = zeros(ncC,2); % efpafpb(c,:) = y' or [efpa efpb];
tfpafpb = zeros(ncC,2); % tfpafpb(c,:) = [fpa fpb];
varefp = zeros(ncC,4); % varfp(c,:) = [covar(1,1) covar(2,2)];
Kp = zeros(ncC,1); % Kp(c,1) = kp;
eKp = zeros(ncC,1); % eKp(c,1) = ekp;
uC = zeros(ncC,1); % uC(c) = uc;
vC = zeros(ncC,1); % vC(c) = vc;
% lambda = zeros(ncC,1);
reg_lambda = zeros(ncC,1); % reg_lambda(c) = reg_corner;
residual = zeros(ncB,ncS); % residual(:,c) = (cS-ecS)./(stdSc(c)*sqrt(ncB*ncS));
chisq_TD = zeros(ncC,1); % chisq_TD(c,1) = ((cS-tcS)'*(cS-tcS))/(varSc(c)*ncB);
% chisq_TE = zeros(ncC,1); % chisq_TE(c,1) = ((ecS-tcS)'*(ecS-tcS))/(varSc(c)*ncB);
chisq_ED = zeros(ncC,1); % chisq_ED(c,1) = ((cS-ecS)'*(cS-ecS))/(varSc(c)*ncB);
% chisq_AD = zeros(ncC,1); % chisq_AD(c,1) = ((cS-aS)'*(cS-aS))/(varSc(c)*ncB);
% chisq_BD = zeros(ncC,1); % chisq_BD(c,1) = ((cS-bS)'*(cS-bS))/(varSc(c)*ncB);
% chisq_AB = zeros(ncC,1); % chisq_AB(c,1) = ((aS-bS)'*(aS-bS))/(varSc(c)*ncB);
% chisq_AE = zeros(ncC,1); % chisq_AE(c,1) = ((aS-ecS)'*(aS-ecS))/(varSc(c)*ncB);
% chisq_BE = zeros(ncC,1); % chisq_BE(c,1) = ((bS-ecS)'*(bS-ecS))/(varSc(c)*ncB);
chisq_ZD = zeros(ncC,1); % chisq_ZD(c,1) = (((cS-zcS).*wfxnS)'*((cS-zcS).*wfxnS))/(varSc(c)*Nfit);

switch boundary.config
    case '2critpts'
        b = boundary.parameterization;
        angle = pi.*x(1);
        slope = angle2slope(angle);
        [u,directrixO,v,directrixA] = parallel2ruled(x,C_bdy,C_coex(25,:),boundary);
        cpt1 = directrixO(1,:);
        cpt2 = directrixO(end,:);
        invariant_pts = struct('cpt1',cpt1,'cpt2',cpt2);
        
        switch part_coeff.config
            case 'constant' % 2critpts
                % x = parameters, length(x) should = 2 = [slope Kp]
                if length(x) ~= 2
                    error('need 2 parameters, [slope Kp]');
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
                    int_pts = boundary_line_int(C_bdy,slope,cC);
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
                    kp = x(2);
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
                    varefp(c,:) = [covar(1,1) covar(2,2) covar(1,2) covar(2,1)];
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
                    chisq_BE(c,1) = (((bS-ecS).*wfxnS)'*((bS-ecS).*wfxnS))/(varSc(c)*Nfit);
                    residual(:,c) = ((cS-tcS).*wfxnS)./(stdSc(c)*sqrt(ncS*Nfit)); %(cS-ecS)./(stdSc(c)*sqrt(ncB*ncS));
%                     chisquare(c,1) = ((((cS-ecS).*wfxnS)'*((cS-ecS).*wfxnS))/(varSc(c)*ncB*ncS));
%                     residual(:,c) = (((cS-ecS).*wfxnS)./(stdSc(c)*sqrt(ncB*ncS)));
                end % end of loop over phase coexistence points
            case 'function' % 2critpts
                % x = parameters, length(x) should = 2 = [slope c]
                if length(x) ~= 2
                    error('need 2 parameters, [slope c]');
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
                    int_pts = boundary_line_int(C_bdy,slope,cC);
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
                    kp = Kp_fxn(x(2),uc,boundary);
                    fpb = (kp*fb)/(fa + kp*fb);
                    fpa = fa/(fa + kp*fb);
                    tcS = fpa*aS + fpb*bS;
                    
                    options = optimset('LargeScale','off','Display','off');
%                     y = lsqlin([aS.*wfxnfp bS.*wfxnfp],cS.*wfxnfp,[-1 0;0 -1],[1e-10;1e-10],[1 1],1,[],[],[],options);
                    y = lsqlin([aS.*wfxnfp bS.*wfxnfp],cS.*wfxnfp,[1 0;0 1],[1-(1e-10);1-(1e-10)],[1 1],1,[],[],[],options);
%                     y = lsqlin([aS.*wfxnfp bS.*wfxnfp],cS.*wfxnfp,[],[],[1 1],1,[],[],[],options);
                    efpa = y(1);
                    efpb = y(2);
                    
%                     efpb = ((bS-aS)'*(cS-aS))/((bS-aS)'*(bS-aS));
%                     efpa = 1-efpb;
                    
                    matrix = ([aS bS]'*[aS bS]);
                    covar = inv(matrix);
                    ecS = efpa*aS + efpb*bS;
                    ekp = (efpb*fa)/(efpa*fb);
                    
                    A = bS-aS;
                    B = cS-aS;
                    [U,s,V] = csvd(A);
                    reg_corner = l_curve(U,s,B);
                    x_0 = fpb;
                    z = tikhonov(U,s,V,B,reg_corner,x_0);
                    zcS = (1-z)*aS + z*bS;
                    reg_lambda(c) = reg_corner;
                    
%                     lambda(c) = (1/(2*(efpb-fpb)))*((bS-aS)'*(cS-aS) - fpb*norm(bS-aS)^2);
                    
                    uC(c) = uc; 
                    vC(c) = vc; 
                    fafb(c,:) = [fa fb];
                    zfpafpb(c,:) = [1-z z];
                    efpafpb(c,:) = [efpa efpb];
                    tfpafpb(c,:) = [fpa fpb];
                    varefp(c,:) = [covar(1,1) covar(2,2) covar(1,2) covar(2,1)];
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
                    regcS(:,c) = zcS;
                    chisq_ZD(c,1) = (((cS-zcS).*wfxnS)'*((cS-zcS).*wfxnS))/(varSc(c)*Nfit);
                    chisq_TD(c,1) = (((cS-tcS).*wfxnS)'*((cS-tcS).*wfxnS))/(varSc(c)*Nfit); %((cS-tcS)'*(cS-tcS))/(varSc(c)*ncB);
%                     chisq_TE(c,1) = (((ecS-tcS).*wfxnS)'*((ecS-tcS).*wfxnS))/(varSc(c)*Nfit); %((ecS-tcS)'*(ecS-tcS))/(varSc(c)*ncB);
                    chisq_ED(c,1) = (((cS-ecS).*wfxnS)'*((cS-ecS).*wfxnS))/(varSc(c)*Nfit); %((cS-ecS)'*(cS-ecS))/(varSc(c)*ncB);
%                     chisq_AD(c,1) = (((cS-aS).*wfxnS)'*((cS-aS).*wfxnS))/(varSc(c)*Nfit); %((cS-aS)'*(cS-aS))/(varSc(c)*ncB);
%                     chisq_BD(c,1) = (((cS-bS).*wfxnS)'*((cS-bS).*wfxnS))/(varSc(c)*Nfit); %((cS-bS)'*(cS-bS))/(varSc(c)*ncB);
%                     chisq_AB(c,1) = (((aS-bS).*wfxnS)'*((aS-bS).*wfxnS))/(varSc(c)*Nfit); %((aS-bS)'*(aS-bS))/(varSc(c)*ncB);
%                     chisq_AE(c,1) = (((aS-ecS).*wfxnS)'*((aS-ecS).*wfxnS))/(varSc(c)*Nfit);
%                     chisq_BE(c,1) = (((bS-ecS).*wfxnS)'*((bS-ecS).*wfxnS))/(varSc(c)*Nfit);
%                     residual(:,c) = ((cS-tcS).*wfxnS)./(stdSc(c)*sqrt(ncS*Nfit)); %(cS-ecS)./(stdSc(c)*sqrt(ncB*ncS));
                    residual(:,c) = ((cS-ecS).*wfxnS)./(stdSc(c)*sqrt(ncS*Nfit)); %(cS-ecS)./(stdSc(c)*sqrt(ncB*ncS));
%                     chisquare(c,1) = ((((cS-ecS).*wfxnS)'*((cS-ecS).*wfxnS))/(varSc(c)*ncB*ncS));
%                     residual(:,c) = (((cS-ecS).*wfxnS)./(stdSc(c)*sqrt(ncB*ncS)));
                end % end of loop over phase coexistence points
            otherwise
                error('invalid partition coefficient configuration');
        end  
    case '1critpt/1endtl'
        b = boundary.parameterization;
        e1 = x(1);
        e2 = x(2);
        % only one configuration
        ept1 = b2bdypt(b,C_bdy,e1);
        ept2 = b2bdypt(b,C_bdy,e2);
        slope = diff([ept1;ept2]);
        slope = slope(2)/slope(1);
        
        [s1,section1,s2,section2] = divide_bdy(e1,e2,C_bdy);% x(1),x(2) = e1,e2
        [s1length,temp1] = chord_length(section1);
        [s2length,temp2] = chord_length(section2);
        clear temp1 temp2;
        if s1length >= s2length
            bdy = [section1;section1(1,:)];
        else
            bdy = [section2;section2(1,:)];
        end
        [IN ON] = inpolygon(C_coex(:,1),C_coex(:,2),bdy(:,1),bdy(:,2));
        Ccoex = C_coex(IN,:);
        Scoex = S_coex(:,IN); 
        
        cin = find(IN);
        nin = length(cin);
        cout = find(~IN);
        nout = length(cout);
        
        [u,directrixO,v,directrixA] = parallel2ruled(x,C_bdy,[],boundary);
        cpt = directrixO(1,:);
        invariant_pts = struct('cpt',cpt,'ept1',ept1,'ept2',ept2);
        
        switch part_coeff.config
            case 'constant' % '1critpt/1endtl'
                % x = parameters, length(x) should = 3 = [ep1 ep2 Kp]
                if length(x) ~= 3
                    error('need 3 parameters, [ep1 ep2 Kp]');
                end
                % loop over phase coexistence points and calculate spectra at each point
                for c = cin'
                    cS = S_coex(:,c);
                    cC = C_coex(c,:);

                    % the intersection of the tie line drawn through the intersection
                    % point of the critical point tangent line(s) and/or end tielines and the
                    % coexistence point with the boundary gives both the composition
                    % point of the alpha phase (C_alpha, or Ld phase) and the
                    % composition of the beta phase (C_beta, or Lo phase).  the
                    % points returned from boundary_line_int are that
                    % C_alpha=int_pts(1,:) and C_beta=int_pts(2,:)
                    int_pts = boundary_line_int(C_bdy,slope,cC);
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
                    varefp(c,:) = [covar(1,1) covar(2,2) covar(1,2) covar(2,1)];
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
                    chisq_BE(c,1) = (((bS-ecS).*wfxnS)'*((bS-ecS).*wfxnS))/(varSc(c)*Nfit);
                    residual(:,c) = ((cS-tcS).*wfxnS)./(stdSc(c)*sqrt(ncS*Nfit)); %(cS-ecS)./(stdSc(c)*sqrt(ncB*ncS));
%                     chisquare(c,1) = ((((cS-ecS).*wfxnS)'*((cS-ecS).*wfxnS))/(varSc(c)*ncB*ncS));
%                     residual(:,c) =
%                     (((cS-ecS).*wfxnS)./(stdSc(c)*sqrt(ncB*ncS)));

%                     uC(c) = uc; 
%                     vC(c) = vc; 
%                     expcS(:,c) = ecS;
%                     efpafpb(c,:) = y';
%                     ofpafpb(c,:) = [fpa fpb];
%                     Kp(c,1) = kp;
%                     alphaC(c,:) = aC;
%                     betaC(c,:) = bC;
%                     coexC(c,:) = cC;
%                     alphaS(:,c) = aS;
%                     betaS(:,c) = bS;
%                     fafb(c,:) = [fa fb];
%                     obscS(:,c) = cS;
%                     varfp(c,:) = [covar(1,1) covar(2,2)];
%                     chisquare(c,1) = ((((cS-ecS).*wfxnS)'*((cS-ecS).*wfxnS))/(varSc(c)*ncB*nin));
%                     residual(:,c) = (((cS-ecS).*sqrt(wfxnS))./(stdSc(c)*sqrt(ncB*nin)));
                end % end of loop over phase coexistence points
                for c = cout'
                    uC(c) = NaN; 
                    vC(c) = NaN; 
                    fafb(c,:) = [NaN NaN];
                    efpafpb(c,:) = [NaN NaN];
                    tfpafpb(c,:) = [NaN NaN];
                    varefp(c,:) = [NaN NaN NaN NaN];
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
                    chisq_AE(c,1) = NaN;
                    chisq_BE(c,1) = NaN;
                    residual(:,c) = 0;
%                     cS = S_coex(:,c);
%                     uC(c) = NaN; 
%                     vC(c) = NaN; 
%                     expcS(:,c) = cS;
%                     efpafpb(c,:) = [NaN NaN];
%                     ofpafpb(c,:) = [NaN NaN];
%                     Kp(c,1) = NaN;
%                     alphaC(c,:) = [NaN NaN];
%                     betaC(c,:) = [NaN NaN];
%                     coexC(c,:) = [NaN NaN];
%                     alphaS(:,c) = repmat(NaN,ncB,1);
%                     betaS(:,c) = repmat(NaN,ncB,1);
%                     fafb(c,:) = [NaN NaN];
%                     obscS(:,c) = cS;
%                     varfp(c,:) = [NaN NaN];
%                     chisquare(c,1) = 0;
%                     residual(:,c) = 0;
                end
            case 'function' % 1critpt/1endtl
                % x = parameters, length(x) should = 4 = [ep1 ep2 c1 c2]
                if length(x) ~= 4
                    error('need 4 parameters, [ep1 ep2 c1 c2]');
                end
                % loop over phase coexistence points and calculate spectra at each point
                for c = cin'
                    cS = S_coex(:,c);
                    cC = C_coex(c,:);

                    % the intersection of the tie line drawn through the intersection
                    % point of the critical point tangent line(s) and/or end tielines and the
                    % coexistence point with the boundary gives both the composition
                    % point of the alpha phase (C_alpha, or Ld phase) and the
                    % composition of the beta phase (C_beta, or Lo phase).  the
                    % points returned from boundary_line_int are that
                    % C_alpha=int_pts(1,:) and C_beta=int_pts(2,:)
                    int_pts = boundary_line_int(C_bdy,slope,cC);
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
                    kp = Kp_fxn(x(3:4),uc,boundary);
                    fpb = (kp*fb)/(fa + kp*fb);
                    fpa = fa/(fa + kp*fb);
                    tcS = fpa*aS + fpb*bS;
                    
                    options = optimset('LargeScale','off','Display','off');
                    y = lsqlin([aS.*wfxnfp bS.*wfxnfp],cS.*wfxnfp,[1 0;0 1],[1-(1e-10);1-(1e-10)],[1 1],1,[],[],[],options);
%                     y = lsqlin([aS.*wfxnfp bS.*wfxnfp],cS.*wfxnfp,[],[],[1 1],1,[],[],[],options);
                    efpa = y(1);
                    efpb = y(2);
                    
%                     efpb = ((bS-aS)'*(cS-aS))/((bS-aS)'*(bS-aS));
%                     efpa = 1-efpb;
                    
                    matrix = ([aS bS]'*[aS bS]);
                    covar = inv(matrix);
                    ecS = efpa*aS + efpb*bS;
                    ekp = (efpb*fa)/(efpa*fb);
                    
                    A = bS-aS;
                    B = cS-aS;
                    [U,s,V] = csvd(A);
                    reg_corner = l_curve(U,s,B);
                    x_0 = fpb;
                    z = tikhonov(U,s,V,B,reg_corner,x_0);
                    zcS = (1-z)*aS + z*bS;
                    reg_lambda(c) = reg_corner;
                    
                    uC(c) = uc; 
                    vC(c) = vc; 
                    fafb(c,:) = [fa fb];
                    zfpafpb(c,:) = [1-z z];
                    efpafpb(c,:) = [efpa efpb];
                    tfpafpb(c,:) = [fpa fpb];
                    varefp(c,:) = [covar(1,1) covar(2,2) covar(1,2) covar(2,1)];
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
                    regcS(:,c) = zcS;
                    chisq_ZD(c,1) = (((cS-zcS).*wfxnS)'*((cS-zcS).*wfxnS))/(varSc(c)*Nfit);
                    chisq_TD(c,1) = (((cS-tcS).*wfxnS)'*((cS-tcS).*wfxnS))/(varSc(c)*Nfit); %((cS-tcS)'*(cS-tcS))/(varSc(c)*ncB);
%                     chisq_TE(c,1) = (((ecS-tcS).*wfxnS)'*((ecS-tcS).*wfxnS))/(varSc(c)*Nfit); %((ecS-tcS)'*(ecS-tcS))/(varSc(c)*ncB);
                    chisq_ED(c,1) = (((cS-ecS).*wfxnS)'*((cS-ecS).*wfxnS))/(varSc(c)*Nfit); %((cS-ecS)'*(cS-ecS))/(varSc(c)*ncB);
%                     chisq_AD(c,1) = (((cS-aS).*wfxnS)'*((cS-aS).*wfxnS))/(varSc(c)*Nfit); %((cS-aS)'*(cS-aS))/(varSc(c)*ncB);
%                     chisq_BD(c,1) = (((cS-bS).*wfxnS)'*((cS-bS).*wfxnS))/(varSc(c)*Nfit); %((cS-bS)'*(cS-bS))/(varSc(c)*ncB);
%                     chisq_AB(c,1) = (((aS-bS).*wfxnS)'*((aS-bS).*wfxnS))/(varSc(c)*Nfit); %((aS-bS)'*(aS-bS))/(varSc(c)*ncB);
%                     chisq_AE(c,1) = (((aS-ecS).*wfxnS)'*((aS-ecS).*wfxnS))/(varSc(c)*Nfit);
%                     chisq_BE(c,1) = (((bS-ecS).*wfxnS)'*((bS-ecS).*wfxnS))/(varSc(c)*Nfit);
%                     residual(:,c) = ((cS-tcS).*wfxnS)./(stdSc(c)*sqrt(ncS*Nfit)); %(cS-ecS)./(stdSc(c)*sqrt(ncB*ncS));
                    residual(:,c) = ((cS-ecS).*wfxnS)./(stdSc(c)*sqrt(ncS*Nfit)); %(cS-ecS)./(stdSc(c)*sqrt(ncB*ncS));
%                     chisquare(c,1) = ((((cS-ecS).*wfxnS)'*((cS-ecS).*wfxnS))/(varSc(c)*ncB*ncS));
%                     residual(:,c) =
%                     (((cS-ecS).*wfxnS)./(stdSc(c)*sqrt(ncB*ncS)));
                end % end of loop over phase coexistence points
                for c = cout'
                    uC(c) = NaN; 
                    vC(c) = NaN; 
                    fafb(c,:) = [NaN NaN];
                    zfpafpb(c,:) = [NaN NaN];
                    efpafpb(c,:) = [NaN NaN];
                    tfpafpb(c,:) = [NaN NaN];
                    varefp(c,:) = [NaN NaN NaN NaN];
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
                    regcS(:,c) = repmat(NaN,ncB,1);
                    chisq_ZD(c,1) = NaN;
                    chisq_TD(c,1) = NaN;
%                     chisq_TE(c,1) = NaN;
%                     chisq_ED(c,1) = NaN;
%                     chisq_AD(c,1) = NaN;
%                     chisq_BD(c,1) = NaN;
%                     chisq_AB(c,1) = NaN;
%                     chisq_AE(c,1) = NaN;
%                     chisq_BE(c,1) = NaN;
                    residual(:,c) = 0;
                end
            otherwise
                error('invalid partition coefficient configuration');
        end 
    case '2endtls'
        switch part_coeff.config
            case 'constant'
                angle = pi.*x(1);
                slope = angle2slope(angle);
                b = boundary.parameterization;
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
                    int_pts = boundary_line_int(C_bdy,slope,cC);
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
                    kp = x(2);
                    fpb = (kp*fb)/(fa + kp*fb);
                    fpa = fa/(fa + kp*fb);
                    ecS = fpa*aS + fpb*bS;
                    
                    options = optimset('LargeScale','off','Display','off');
                    y = lsqlin([aS.*wfxnfp bS.*wfxnfp],cS.*wfxnfp,[-1 0;0 -1],[0;0],[1 1],1,[],[],[],options);
                    matrix = ([aS bS]'*[aS bS]);
                    covar = inv(matrix);
                    
                    if nargout > 1
                        condn(c,1) = cond([aS bS]);
                        expcS(:,c) = ecS;
                        expcS2fit(:,c) = (ecS./(spectra.cSbline_std(c)*sqrt(ncB)*sqrt(nin))).*sqrt(spectra.weight_fxn);
                %         chisquare(c,1) = ((cS-ecS)'*(cS-ecS))/((ncB-1)*mean(spectra.cSbline_var));
%                         chisquare(c,1) = ((cS-ecS)'*(cS-ecS))/((ncB-1)*spectra.cSbline_var(c));
                %         chisquare(c,1) = ((cS-ecS)'*(cS-ecS));
                        obscS2fit(:,c) = (cS./(spectra.cSbline_std(c)*sqrt(ncB)*sqrt(nin))).*sqrt(spectra.weight_fxn);
                        chisquare(c,1) = (obscS2fit(:,c)-expcS2fit(:,c))'*(obscS2fit(:,c)-expcS2fit(:,c));
                        efpafpb(c,:) = y';
                        ofpafpb(c,:) = [fpa fpb];
                        Kp(c,1) = kp;
                        alphaC(c,:) = aC;
                        betaC(c,:) = bC;
                        coexC(c,:) = cC;
                        alphaS(:,c) = aS;
                        betaS(:,c) = bS;
                        fafb(c,:) = [fa fb];
                        obscS(:,c) = cS;
                    else
                        residual(:,c) = ((cS./(spectra.cSbline_std(c)*sqrt(ncB)*sqrt(nin))).*sqrt(spectra.weight_fxn))...
                            - ((ecS./(spectra.cSbline_std(c)*sqrt(ncB)*sqrt(nin))).*sqrt(spectra.weight_fxn));
                    end
                end % end of loop over phase coexistence points
            case 'function' % 2endtls
                angle = pi.*x(1);
                slope = angle2slope(angle);
                b = boundary.parameterization;
                [u,directrixO] = parallel2ruled(x,C_bdy,C_coex(25,:),boundary);
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
                    int_pts = boundary_line_int(C_bdy,slope,cC);
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
                    kp = Kp_fxn(x(2:end),uc,boundary);
                    fpb = (kp*fb)/(fa + kp*fb);
                    fpa = fa/(fa + kp*fb);
                    ecS = fpa*aS + fpb*bS;
                    
                    options = optimset('LargeScale','off','Display','off');
                    y = lsqlin([aS.*wfxnfp bS.*wfxnfp],cS.*wfxnfp,[-1 0;0 -1],[0;0],[1 1],1,[],[],[],options);
                    matrix = ([aS bS]'*[aS bS]);
                    covar = inv(matrix);
                    
                    if nargout > 1
                        condn(c,1) = cond([aS bS]);
                        expcS(:,c) = ecS;
                        expcS2fit(:,c) = (ecS./(spectra.cSbline_std(c)*sqrt(ncB)*sqrt(ncS))).*sqrt(spectra.weight_fxn);
                %         chisquare(c,1) = ((cS-ecS)'*(cS-ecS))/((ncB-1)*mean(spectra.cSbline_var));
%                         chisquare(c,1) = ((cS-ecS)'*(cS-ecS))/((ncB-1)*spectra.cSbline_var(c));
                %         chisquare(c,1) = ((cS-ecS)'*(cS-ecS));
                        obscS2fit(:,c) = (cS./(spectra.cSbline_std(c)*sqrt(ncB)*sqrt(ncS))).*sqrt(spectra.weight_fxn);
                        chisquare(c,1) = (obscS2fit(:,c)-expcS2fit(:,c))'*(obscS2fit(:,c)-expcS2fit(:,c));
                        efpafpb(c,:) = y';
                        ofpafpb(c,:) = [fpa fpb];
                        Kp(c,1) = kp;
                        alphaC(c,:) = aC;
                        betaC(c,:) = bC;
                        coexC(c,:) = cC;
                        alphaS(:,c) = aS;
                        betaS(:,c) = bS;
                        fafb(c,:) = [fa fb];
                        obscS(:,c) = cS;
                    else
                        residual(:,c) = ((cS./(spectra.cSbline_std(c)*sqrt(ncB)*sqrt(nin))).*sqrt(spectra.weight_fxn))...
                            - ((ecS./(spectra.cSbline_std(c)*sqrt(ncB)*sqrt(nin))).*sqrt(spectra.weight_fxn));
                    end
                end % end of loop over phase coexistence points
            otherwise
                error('invalid partition coefficient configuration');
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
% X2_TE = nanmean(chisq_TE);
X2_ED = nanmean(chisq_ED);
% Area_DET = nanmean(X2area_triDET);
% X2_AD = nanmean(chisq_AD);
% X2_BD = nanmean(chisq_BD);
% X2_AB = nanmean(chisq_AB);
% X2_AE = nanmean(chisq_AE);
% X2_BE = nanmean(chisq_BE);
X2_ZD = nanmean(chisq_ZD);
% Area_DAB = nanmean(X2area_triDAB);
% projT = diag((obscS-theorycS)'*(obscS-theorycS));
% avg_projT = nanmean(projT);
% projA = diag((obscS-theorycS)'*(alphaS-theorycS));
% avg_projA = nanmean(projA);
% projB = diag((obscS-theorycS)'*(betaS-theorycS));
% avg_projB = nanmean(projB);
% chi2fpa = ((efpafpb(:,1)-tfpafpb(:,1)).^2)./varefp(:,1);
% chi2fpb = ((efpafpb(:,2)-tfpafpb(:,2)).^2)./varefp(:,2);
% chi2fp = chi2fpa+chi2fpb;
% Chi2fp = nanmean(chi2fp);
normfpz = norm(zfpafpb(~isnan(zfpafpb(:)))-tfpafpb(~isnan(tfpafpb(:))));
normfp = norm(efpafpb(~isnan(efpafpb(:)))-tfpafpb(~isnan(tfpafpb(:))));
normfp2 = sqrt(sum((efpafpb-tfpafpb).^2,2));
normfp2z = sqrt(sum((zfpafpb-tfpafpb).^2,2));
Normfp2 = nanmean(normfp2);
Normfp2z = nanmean(normfp2z);
corrcoeff = varefp(:,3)./(sqrt(varefp(:,1).*varefp(:,2)));
% CC = nanmean(corrcoeff);

% output depends on the number of output arguments
if nargout > 1           
    % output fit result data structure     
%     varargout{1} = X2_TD*normfp;
%     varargout{1} = X2_TD;
%     varargout{1} = X2_TD*normfp./abs(CC);
%     varargout{1} = nanmean(chisq_ED.*normfp2./abs(corrcoeff));
%     varargout{1} = nanmean(chisq_ZD + (reg_lambda.*sqrt(sum((zfpafpb(:,2)-tfpafpb(:,2)).^2,2))).^2);
%     varargout{1} = X2_ED;
%     varargout{1} = X2_ED*normfp;
%     varargout{1} = nanmean(chisq_ED.*normfp2);
%     varargout{1} = X2_ED.*Normfp2;
    fit_measure = X2_TD;
    
    varargout{1} = fit_measure;
    
    configs = struct('bdyconfig',boundary.config,'tlconfig',{'parallel'},'Kpconfig',part_coeff.config);
    fit_data = struct('x',x','fit_measure',fit_measure,'chisq_TD',chisq_TD,'chisq_ED',chisq_ED,'chisq_ZD',chisq_ZD,...
        'X2_TD',X2_TD,'X2_ED',X2_ED,'X2_ZD',X2_ZD,'Nfit',Nfit);
    comps_data = struct('invariant_pts',invariant_pts,'aC',alphaC,'bC',betaC,'cC',coexC,'bdy',C_bdy);
    spectra_data = struct('magnetic_fields',spectra.magnetic_fields,...
        'aS',alphaS,'bS',betaS,'exp_cS',expcS,'obs_cS',obscS,'theory_cS',theorycS,'regu_cS',regcS,...
        'cSbline_var',spectra.cSbline_var','cSbline_std',spectra.cSbline_std','bSC_std',spectra.bSC_std,'bSC_var',spectra.bSC_var,...
        'weight_fxn_S',spectra.weight_fxn_S,'weight_fxn_fp',spectra.weight_fxn_fp);
    tieline_data = struct('efpafpb',efpafpb,'tfpafpb',tfpafpb,'varefp',varefp,...
        'normfp',normfp,'normfpz',normfpz,'normfp2',normfp2,'normfp2z',normfp2z,'Normfp2',Normfp2,'Normfp2z',Normfp2z,...
        'fafb',fafb,'Kp',Kp,'eKp',eKp,'uC',uC,'vC',vC,'corrcoeff',corrcoeff,'zfpafpb',zfpafpb,'reg_lambda',reg_lambda,...
        'directrixO',directrixO,'directrixA',directrixA);
    
    varargout{2} = struct('configs',configs,'fit',fit_data,'comps',comps_data,'spectra',spectra_data,'tieline',tieline_data);
    
%     configs = struct('bdyconfig',boundary.config,'tlconfig',{'parallel'},'Kpconfig',part_coeff.config);
%     fit_data = struct('x',x','chisq_TD',chisq_TD,'chisq_TE',chisq_TE,'chisq_ED',chisq_ED,'chisq_ZD',chisq_ZD,...
%         'chisq_AD',chisq_AD,'chisq_BD',chisq_BD,'chisq_AB',chisq_AB,'chisq_AE',chisq_AE,'chisq_BE',chisq_BE,...
%         'X2_TD',X2_TD,'X2_TE',X2_TE,'X2_ED',X2_ED,'X2_AD',X2_AD,'X2_BD',X2_BD,'X2_AB',X2_AB,'X2_AE',X2_AE,'X2_BE',X2_BE,'X2_ZD',X2_ZD,...
%         'projA',projA,'projB',projB,'projT',projT,'Chi2fp',Chi2fp,'chi2fp',chi2fp,'chi2fpa',chi2fpa,'chi2fpb',chi2fpb);
%     comps_data = struct('invariant_pts',invariant_pts,'aC',alphaC,'bC',betaC,'cC',coexC,'bdy',C_bdy);
%     spectra_data = struct('magnetic_fields',spectra.magnetic_fields,...
%         'aS',alphaS,'bS',betaS,'exp_cS',expcS,'obs_cS',obscS,'theory_cS',theorycS,'regu_cS',regcS,...
%         'cSbline_var',spectra.cSbline_var','cSbline_std',spectra.cSbline_std','bSC_std',spectra.bSC_std,'bSC_var',spectra.bSC_var,...
%         'weight_fxn_S',spectra.weight_fxn_S,'weight_fxn_fp',spectra.weight_fxn_fp);
%     tieline_data = struct('efpafpb',efpafpb,'tfpafpb',tfpafpb,'varefp',varefp,'normfp',normfp,'normfpz',normfpz,'normfp2',normfp2,'normfp2z',normfp2z,...
%         'fafb',fafb,'Kp',Kp,'eKp',eKp,'uC',uC,'vC',vC,'corrcoeff',corrcoeff,'lambda',lambda,'zfpafpb',zfpafpb,'reg_lambda',reg_lambda);
%     varargout{2} = struct('configs',configs,'fit',fit_data,'comps',comps_data,'spectra',spectra_data,'tieline',tieline_data);
    
%     varargout{2} = struct('bdyconfig',boundary.config,'tlconfig',{'parallel'},'Kpconfig',part_coeff.config,'x',x',...
%         'chisq_TD',chisq_TD,'chisq_TE',chisq_TE,'chisq_ED',chisq_ED,'X2area_triDET',X2area_triDET,'Area_DET',Area_DET,...
%         'chisq_AD',chisq_AD,'chisq_BD',chisq_BD,'chisq_AB',chisq_AB,'X2area_triDAB',X2area_triDAB,'Area_DAB',Area_DAB,...
%         'projA',projA,'avg_projA',avg_projA,'projB',projB,'avg_projB',avg_projB,'projT',projT,'avg_projT',avg_projT,...
%         'Chi2fp',Chi2fp,'chi2fpa',chi2fpa,'chi2fpb',chi2fpb,'normfp',normfp,'normfp2',normfp2,'Normfp2',Normfp2,...
%         'invariant_pts',invariant_pts,'aC',alphaC,'bC',betaC,'cC',coexC,...
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
%             varargout{1} = reshape(residual.*repmat(sqrt(normfp2),ncB,1),ncB*ncS,1)
            temp = reshape(residual*sqrt(normfp),ncB*ncS,1);
            varargout{1} = temp(~isnan(temp));
        case 'levmar'
        case 'simplex'
%             varargout{1} = nanmean(chisq_TE.*chisq_TD.*chisq_ED.*(normfp./chisq_AB));
%             varargout{1} = ((X2_TD*X2_ED*X2_TE)^1/3)/X2_AB; % 48 degrees
%             varargout{1} = (X2_TD*X2_ED*X2_TE)^1/3; % 63 degrees
%             varargout{1} = nanmean(((chisq_TE.*chisq_TD.*chisq_ED).^1/3)); % 113 degrees
%             varargout{1} = nanmean(((chisq_TE.*chisq_TD.*chisq_ED).^1/3))./X2_AB; % 50 degrees
%             varargout{1} = nanmean(((chisq_TE.*chisq_TD.*chisq_ED).^1/3)./chisq_AB);
%             varargout{1} = nanmean(chisq_TD.*(chisq_TE./chisq_AB));
%             varargout{1} = X2_TD*X2_ED*(X2_TE/X2_AB);
%             varargout{1} = X2_ED;
%             varargout{1} = X2_ED*normfp;
%             varargout{1} = nanmean(chisq_ED.*normfp2);
%             varargout{1} = X2_ED.*Normfp2;
%             varargout{1} = X2_TD*normfp;
%             varargout{1} = X2_TD*(1+normfp);
%             varargout{1} = Normfp2;
%             varargout{1} = sqrt(X2_TE*X2_TD);
%             varargout{1} = Area_DET;
            varargout{1} = X2_TD;
%             varargout{1} = normfp;
%             varargout{1} = X2_TD*normfp./abs(CC);
%             varargout{1} = nanmean(chisq_ED.*normfp2./abs(corrcoeff));
%             varargout{1} = nanmean(chisq_ZD +(reg_lambda.*sqrt(sum((zfpafpb(:,2)-tfpafpb(:,2)).^2,2))).^2);
        otherwise
            error('invalid fit method');
    end
end

return