function varargout = ruled_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method)
% varargout =
% ruled_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,fit_method)

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
varefp = zeros(ncC,4); % varfp(c,:) = [covar(1,1) covar(2,2)];
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
        cpt1 = b2bdypt(b,C_bdy,x(1));
        cpt2 = b2bdypt(b,C_bdy,x(2));
        invariant_pts = struct('cpt1',cpt1,'cpt2',cpt2);
        
        switch part_coeff.config
            case 'constant'
                % x = parameters, length(x) should = 4 = [cp1 cp2 gamma Kp]
                if length(x) ~= 4
                    error('need 4 parameters, [cp1 cp2 gamma Kp]');
                end
                % loop over phase coexistence points and calculate spectra at each point
                for c = 1:ncC
                    cS = S_coex(:,c);
                    cC = C_coex(c,:);

                    [uc,vc] = coexpt2ruledp(x(3),cC,[u directrixO],[v directrixA]); % v = gamma
                    bC = interp1(u,directrixO,uc);
                    aC = interp1(v,directrixA,vc);

                    % interpolate spectra at alpha point on the boundary
                    ab = bdypt2b(aC,b,C_bdy); % bi = bdypt2b(bdypt,b,bdy) 
                %     ap = bdypt2param(C_bdy,aC,'polygon',b);
                    aS = interp1(b,S_bdy',ab);
                    aS = aS';   
                    % normalize spectrum
                    aS = normalize_spectra([spectra.magnetic_fields aS],'d','area',1000);
                    aS = aS(:,2);
                    % interpolate spectra at beta point on the boundary
                    ob = bdypt2b(bC,b,C_bdy);
                    bS = interp1(b,S_bdy',ob);
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
                % x = parameters, length(x) should = 4 = [cp1 cp2 gamma c]
                if length(x) ~= 4
                    error('need 4 parameters, [cp1 cp2 gamma c]');
                end
                % loop over phase coexistence points and calculate spectra at each point
                for c = 1:ncC
                    cS = S_coex(:,c);
                    cC = C_coex(c,:);

                    [uc,vc] = coexpt2ruledp(x(3),cC,[u directrixO],[v directrixA]); % v = gamma
                    bC = interp1(u,directrixO,uc);
                    aC = interp1(v,directrixA,vc);

                    % interpolate spectra at alpha point on the boundary
                    ab = bdypt2b(aC,b,C_bdy); % bi = bdypt2b(bdypt,b,bdy) 
                %     ap = bdypt2param(C_bdy,aC,'polygon',b);
                    aS = interp1(b,S_bdy',ab);
                    aS = aS';   
                    % normalize spectrum
                    aS = normalize_spectra([spectra.magnetic_fields aS],'d','area',1000);
                    aS = aS(:,2);
                    % interpolate spectra at beta point on the boundary
                    ob = bdypt2b(bC,b,C_bdy);
                    bS = interp1(b,S_bdy',ob);
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
                    kp = Kp_fxn(x(4),uc,boundary);
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
                    residual(:,c) = ((cS-ecS).*wfxnS)./(stdSc(c)*sqrt(ncS*Nfit)); %(cS-ecS)./(stdSc(c)*sqrt(ncB*ncS));
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
            [s1,section1,s2,section2] = divide_bdy(x(2),x(3),C_bdy);% x(2),x(3) = e1,e2
%             [u,directrixO,v,directrixA] = tangent2ruled(e1,e2,C_bdy);
            if cp < e1 && cp > e2 % configuration D in notes
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
            elseif cp > e1 || cp < e2 % configuration B in notes
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
            if cp > e1 && cp < e2 % configuration C in notes
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
            elseif cp < e1 || cp > e2 % configuration A in notes
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
        
        ept1 = b2bdypt(b,C_bdy,e1);
        ept2 = b2bdypt(b,C_bdy,e2);
        invariant_pts = struct('cpt',cpt,'ept1',ept1,'ept2',ept2);
        
        switch part_coeff.config
            case 'constant' % '1critpt/1endtl'
                % x = parameters, length(x) should = 5 = [cp ep1 ep2 gamma Kp]
                if length(x) ~= 5
                    error('need 5 parameters, [cp ep1 ep2 gamma Kp]');
                end
                % loop over phase coexistence points and calculate spectra at each point
                for c = cin' % 1:ncC
                    cS = S_coex(:,c);
                    cC = C_coex(c,:);

                    [uc,vc] = coexpt2ruledp(x(4),cC,[u directrixO],[v directrixA]); % v = gamma
                    bC = interp1(u,directrixO,uc);
                    aC = interp1(v,directrixA,vc);

                    % interpolate spectra at alpha point on the boundary
                    ab = bdypt2b(aC,b,C_bdy); % bi = bdypt2b(bdypt,b,bdy) 
                %     ap = bdypt2param(C_bdy,aC,'polygon',b);
                    aS = interp1(b,S_bdy',ab);
                    aS = aS';   
                    % normalize spectrum
                    aS = normalize_spectra([spectra.magnetic_fields aS],'d','area',1000);
                    aS = aS(:,2);
                    % interpolate spectra at beta point on the boundary
                    ob = bdypt2b(bC,b,C_bdy);
                    bS = interp1(b,S_bdy',ob);
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
                    kp = x(5);
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
                end
            case 'function' % 1critpt/1endtl
                % x = parameters, length(x) should = 6 = [cp ep1 ep2 gamma c1 c2]
                if length(x) ~= 6
                    error('need 6 parameters, [cp ep1 ep2 gamma c1 c2]');
                end
                % loop over phase coexistence points and calculate spectra at each point
                for c = cin' % 1:ncC
                    cS = S_coex(:,c);
                    cC = C_coex(c,:);

                    [uc,vc] = coexpt2ruledp(x(4),cC,[u directrixO],[v directrixA]); % v = gamma
                    bC = interp1(u,directrixO,uc);
                    aC = interp1(v,directrixA,vc);

                    % interpolate spectra at alpha point on the boundary
                    ab = bdypt2b(aC,b,C_bdy); % bi = bdypt2b(bdypt,b,bdy) 
                %     ap = bdypt2param(C_bdy,aC,'polygon',b);
                    aS = interp1(b,S_bdy',ab);
                    aS = aS';   
                    % normalize spectrum
                    aS = normalize_spectra([spectra.magnetic_fields aS],'d','area',1000);
                    aS = aS(:,2);
                    % interpolate spectra at beta point on the boundary
                    ob = bdypt2b(bC,b,C_bdy);
                    bS = interp1(b,S_bdy',ob);
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
                    kp = Kp_fxn(x(5:6),uc,boundary);
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
                    
                    % another way to calculate covariance above
%                     [u,w,v] = svd(([aS bS]./spectra.cSbline_std(c)^2),0);
%                     w = diag(w);
%                     for i = 1:2
%                         for k = 1:2
%                             covar(i,k) = sum((v(i,:).*v(k,:))'./(w(:).^2));
%                         end
%                     end
                    
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
                    residual(:,c) = ((cS-ecS).*wfxnS)./(stdSc(c)*sqrt(ncS*Nfit)); %(cS-ecS)./(stdSc(c)*sqrt(ncB*ncS));
%                     chisquare(c,1) = ((((cS-ecS).*wfxnS)'*((cS-ecS).*wfxnS))/(varSc(c)*ncB*ncS));
%                     residual(:,c) =
%                     (((cS-ecS).*wfxnS)./(stdSc(c)*sqrt(ncB*ncS)));
                
%                     chisquare(c,1) = (1-lambda)*(((cS-ecS).*wfxn)'*((cS-ecS).*wfxn))/(varSc(c)*ncB*nin) + lambda*((fpa-y(1))^2/(covar(1,1)*nin)) + lambda*((fpb-y(2))^2/(covar(2,2)*nin));% sum
%                     chisquare(c,1) = ((((cS-ecS).*wfxn)'*((cS-ecS).*wfxn))/(varSc(c)*ncB*nin))*(((fpa-y(1))^2/(covar(1,1)*nin)) + ((fpb-y(2))^2/(covar(2,2)*nin)));% product
%                     residual(:,c) = [(1-lambda)*(((cS-ecS).*sqrt(wfxn))./(stdSc(c)*sqrt(ncB*nin)));lambda*((fpa-y(1))/(sqrt(covar(1,1)*nin)));lambda*((fpb-y(2))/(sqrt(covar(2,2)*nin)))];% sum
%                     residual(c,1) = sqrt(((((cS-ecS).*wfxn)'*((cS-ecS).*wfxn))/(varSc(c)*ncB*nin))*(((fpa-y(1))^2/(covar(1,1)*nin)) + ((fpb-y(2))^2/(covar(2,2)*nin))));% product
%                     residual(:,c) = (((cS-ecS).*sqrt(wfxn))./(stdSc(c)*sqrt(ncB*nin)))*sqrt((((fpa-y(1))^2/(covar(1,1)*nin)) + ((fpb-y(2))^2/(covar(2,2)*nin))));% product
%                     residual(:,c) = [((fpa-y(1))/(sqrt(covar(1,1)*nin)));((fpb-y(2))/(sqrt(covar(2,2)*nin)))];
%                     residual(:,c) = ((fpa-y(1))/(sqrt(covar(1,1)*nin)));
%                     residual(:,c) = (((cS-ecS).*sqrt(wfxnS))./(stdSc(c)*sqrt(ncB*nin))).*((fpa-y(1))/(sqrt(covar(1,1)))).*((fpb-y(2))/(sqrt(covar(2,2))));
%                     residual(:,c) = (((cS-ecS))./(stdSc(c))).*((fpa-y(1))/(sqrt(covar(1,1))));
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
                end
            otherwise 
                error('invalid partition coefficient configuration');
        end  
    case '2endtls'
        switch part.coeff.config
            case 'constant'
                [Olength,u] = chord_length(directrixO);
                [Alength,v] = chord_length(directrixA);
                % loop over phase coexistence points and calculate spectra at each point
                for c = 1:ncC
                    cS = S_coex(:,c);
                    cC = C_coex(c,:);

                    [uc,vc] = coexpt2ruledp(x(3),cC,[u directrixO],[v directrixA]); % v = gamma
                    bC = interp1(u,directrixO,uc);
                    aC = interp1(v,directrixA,vc);

                    % interpolate spectra at alpha point on the boundary
                    ab = bdypt2b(aC,b,C_bdy); % bi = bdypt2b(bdypt,b,bdy) 
                %     ap = bdypt2param(C_bdy,aC,'polygon',b);
                    aS = interp1(b,S_bdy',ab);
                    aS = aS';   
                    % normalize spectrum
                    aS = normalize_spectra([spectra.magnetic_fields aS],'d','area',1000);
                    aS = aS(:,2);
                    % interpolate spectra at beta point on the boundary
                    ob = bdypt2b(bC,b,C_bdy);
                    bS = interp1(b,S_bdy',ob);
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
                    kp = x(4);
                    fpb = (kp*fb)/(fa + kp*fb);
                    fpa = fa/(fa + kp*fb);
                    ecS = fpa*aS + fpb*bS;

                    options = optimset('LargeScale','off','Display','off');
                    y = lsqlin([aS bS],cS,[-1 0;0 -1],[0;0],[1 1],1,[],[],[],options);
                    
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
            case 'function'% 2endtls
                [Olength,u] = chord_length(directrixO);
                [Alength,v] = chord_length(directrixA);
                % loop over phase coexistence points and calculate spectra at each point
                for c = 1:ncC
                    cS = S_coex(:,c);
                    cC = C_coex(c,:);

                    [uc,vc] = coexpt2ruledp(x(3),cC,[u directrixO],[v directrixA]); % v = gamma
                    bC = interp1(u,directrixO,uc);
                    aC = interp1(v,directrixA,vc);

                    % interpolate spectra at alpha point on the boundary
                    ab = bdypt2b(aC,b,C_bdy); % bi = bdypt2b(bdypt,b,bdy) 
                %     ap = bdypt2param(C_bdy,aC,'polygon',b);
                    aS = interp1(b,S_bdy',ab);
                    aS = aS';   
                    % normalize spectrum
                    aS = normalize_spectra([spectra.magnetic_fields aS],'d','area',1000);
                    aS = aS(:,2);
                    % interpolate spectra at beta point on the boundary
                    ob = bdypt2b(bC,b,C_bdy);
                    bS = interp1(b,S_bdy',ob);
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
                    kp = Kp_fxn(x(4:end),uc,boundary);
                    fpb = (kp*fb)/(fa + kp*fb);
                    fpa = fa/(fa + kp*fb);
                    ecS = fpa*aS + fpb*bS;
                    
                    options = optimset('LargeScale','off','Display','off');
                    y = lsqlin([aS bS],cS,[-1 0;0 -1],[0;0],[1 1],1,[],[],[],options);

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
X2_TE = nanmean(chisq_TE);
X2_ED = nanmean(chisq_ED);
% Area_DET = nanmean(X2area_triDET);
X2_AD = nanmean(chisq_AD);
X2_BD = nanmean(chisq_BD);
X2_AB = nanmean(chisq_AB);
X2_AE = nanmean(chisq_AE);
X2_BE = nanmean(chisq_BE);
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
kpnorm = norm(Kp-eKp);
corrcoeff = varefp(:,3)./(sqrt(varefp(:,1).*varefp(:,2)));
CC = nanmean(corrcoeff);

% output depends on the number of output arguments
if nargout > 1
    % output fit result data structure     
%     varargout{1} = X2_TD;
%     varargout{1} = X2_TD./abs(CC);
%     varargout{1} = X2_TD*normfp;
%     varargout{1} = X2_TD*kpnorm;
%     varargout{1} = X2_TD*normfp./abs(CC);
%     varargout{1} = nanmean(chisq_ED.*normfp2./abs(corrcoeff));
%     varargout{1} = X2_ED*normfp;
%     varargout{1} = X2_ED.*Normfp2;
    fit_measure = X2_TD/X2_AB;
    
    varargout{1} = fit_measure;
    
    configs = struct('bdyconfig',boundary.config,'tlconfig',{'ruled'},'Kpconfig',part_coeff.config);
    fit_data = struct('x',x','fit_measure',fit_measure,'chisq_TD',chisq_TD,'chisq_TE',chisq_TE,'chisq_ED',chisq_ED,...
        'chisq_AD',chisq_AD,'chisq_BD',chisq_BD,'chisq_AB',chisq_AB,'chisq_AE',chisq_AE,'chisq_BE',chisq_BE,'Nfit',Nfit,...
        'X2_TD',X2_TD,'X2_TE',X2_TE,'X2_ED',X2_ED,'X2_AD',X2_AD,'X2_BD',X2_BD,'X2_AB',X2_AB,'X2_AE',X2_AE,'X2_BE',X2_BE,...
        'projA',projA,'projB',projB,'projT',projT,'Chi2fp',Chi2fp,'chi2fp',chi2fp,'chi2fpa',chi2fpa,'chi2fpb',chi2fpb);
    comps_data = struct('invariant_pts',invariant_pts,'aC',alphaC,'bC',betaC,'cC',coexC,'bdy',C_bdy);
    spectra_data = struct('magnetic_fields',spectra.magnetic_fields,...
        'aS',alphaS,'bS',betaS,'exp_cS',expcS,'obs_cS',obscS,'theory_cS',theorycS,...
        'cSbline_var',spectra.cSbline_var','cSbline_std',spectra.cSbline_std','bSC_std',spectra.bSC_std,'bSC_var',spectra.bSC_var,...
        'weight_fxn_S',spectra.weight_fxn_S,'weight_fxn_fp',spectra.weight_fxn_fp);
    tieline_data = struct('efpafpb',efpafpb,'tfpafpb',tfpafpb,'varefp',varefp,'normfp',normfp,'normfp2',normfp2,'Normfp2',Normfp2,...
        'fafb',fafb,'Kp',Kp,'eKp',eKp,'uC',uC,'vC',vC,'corrcoeff',corrcoeff,'directrixO',directrixO,'directrixA',directrixA);
    varargout{2} = struct('configs',configs,'fit',fit_data,'comps',comps_data,'spectra',spectra_data,'tieline',tieline_data);
    
%     varargout{2} = struct('bdyconfig',boundary.config,'tlconfig',{'ruled'},'Kpconfig',part_coeff.config,'x',x',...
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
            temp = reshape(residual*sqrt(normfp),ncB*ncS,1);
            varargout{1} = temp(~isnan(temp));
        case 'levmar'
        case 'simplex'
%             varargout{1} = X2_TD;
            varargout{1} = X2_TD/X2_AB;
%             varargout{1} = X2_TD./abs(CC);
%             varargout{1} = (X2_TD*X2_ED*X2_TE)^1/3;
%             varargout{1} = X2_TD*normfp;
%             varargout{1} = X2_TD*kpnorm;
%             varargout{1} = X2_TD*normfp./abs(CC);
%             varargout{1} = nanmean(chisq_ED.*normfp2./abs(corrcoeff));
%             varargout{1} = X2_ED*normfp;
%             varargout{1} = X2_ED.*Normfp2;
        otherwise
            error('invalid fit method');
    end
end

return