function fit_out = lincomb_fit(S,C,varargin)
% lincomb_fit performs a fit of spectra in S by linear combination of the
% first and last spectra in S.
% S = spectra, fit spectra = S(2:end-1), basis spectra = S(1,end) 
% C = compositions of the spectra
% varargin = 'method','method_options'

% checking input arguments

disp(sprintf('checking input arguments and spectrum pre-processing:'));

if all(size(S) > 1) & ndims(S) == 2 % if S is a matrix
    [nb,ncol] = size(S);

    if ncol < 2 || rem(ncol,2) ~= 0
        error('each spectrum is two columns: [B-field intensity_values]');
    end
    
    ns = ncol/2;
else
    error('first argument must be a matrix of the spectra');
end

if  any(size(C) == 1) & any(size(C) > 1) % if C is a vector
    nc = length(C); % nc = number of compositions, should equal ns
    nd = 1;
    [nrow,ncol] = size(C);
    
    if nrow == 1
        C = C'; % convert to column vector
    end
    
    if nc ~= ns
        error('number of compositions must equal number of spectra');
    end
elseif all(size(C) > 1) & ndims(C) == 2 % if C is a matrix 
    [nc,nd] = size(C);
    
    if all(sum(C,2) ~= 1)
        error('sum of compositions must equal 1');
    end
    
    if nc ~= ns
        error('number of compositions must equal number of spectra');
    end
end

% check variable arguments

% align spectra
S = align_spectra(S);

% normalize spectra
S = normalize_spectra(S,'d','area',1000);

% reduce spectra to just intensity values
B = S(:,1:2:end);
I = S(:,2:2:end);

% get baseline statistics
for s = 1:ns
    baseline(:,s) = I([1:200 end-200:end]',s);
    avg_bline(s) = mean(baseline(:,s));
    var_bline(s) = var(baseline(:,s));
    std_bline(s) = std(baseline(:,s));
end 

% fS = S(:,2:end-1);
% [nfb,nfs] = size(fS);
% std_fbline = std_bline(2:end-1);
% define basis
Basis = [I(:,1) I(:,end)];
C_alpha = C(1,:);
C_beta = C(end,:);
[nb,ns] = size(I);

% start of fitting process
disp(sprintf('starting fitting process:'));
switch varargin{1}
    case 'fpafpb'
        switch varargin{2} 
            % solve for the fraction of probe in each phase
            case 'lsqlin'
                options = optimset('LargeScale','off','Display','off');
                for s = 1:ns
                    f_beta(s,1) = norm(C(s,:)-C_alpha)/norm(C_beta-C_alpha);
                    fpafpb(:,s) = lsqlin(Basis,I(:,s),[1 0;0 1],[1-(1e-10);1-(1e-10)],[1 1],1,[],[],[],options);
%                     fpafpb(:,s) = lsqlin(Basis,I(:,s),[-1 0;0 -1],[0;0],[1 1],1,[],[],[],options);
                    if s == 1 || s == ns
                        kp(s) = 1;
                    else
                        kp(s) = (fpafpb(2,s)*(1-f_beta(s,1)))/(fpafpb(1,s)*f_beta(s,1));
                    end
                end
                fpafpb = fpafpb';
                covar = inv(Basis'*Basis);
                fpvar = [covar(1,1) covar(2,2) covar(1,2) covar(2,1)];
            case 'svd'
                options = optimset('LargeScale','off','Display','off');
                for s = 1:ns
                    f_beta(s,1) = norm(C(s,:)-C_alpha)/norm(C_beta-C_alpha);
                    [U,W,V] = svd(Basis,0);
                    UB = U(:,1:2); % B = 2 component eigenbasis
		
                    % solve for eigenspectra coefficient matrix M by linear least squares
                    % constraints
                    % Aeq = [sum(U(:,1)) sum(U(:,2))];
                    Aeq = sum(UB);
                    beq1 = sum(Basis(:,1));
                    beq2 = sum(Basis(:,2));
                    %A = -B;
                    %b = zeros(size(U(:,1)));
                    temp1 = lsqlin(UB,Basis(:,1),[],[],Aeq,beq1,[],[],[],options);
                    temp2 = lsqlin(UB,Basis(:,2),[],[],Aeq,beq2,[],[],[],options);
                    M = [temp1 temp2];
                    Aeq = [1 1];
                    beq = 1;
                    ISt = I(:,s);
                    fpafpb(:,s) = lsqlin(M,[ISt'*UB(:,1);ISt'*UB(:,2)],[],[],Aeq,beq,[],[],[],options);
%                     x = inv(M)*[ISt'*B(:,1);ISt'*B(:,2)];
%                     fpalpha(s) = x(1);
%                     fpbeta(s) = x(2);
                    kp(s) = (fpafpb(2,s)*(1-f_beta(s,1)))/(fpafpb(1,s)*f_beta(s,1));
                end
                fpafpb = fpafpb';
                [u,w,v] = svd((Basis./var_bline.^2),0);
                w = diag(w);
                for i = 1:2
                    for k = 1:2
                        covar(i,k) = sum((v(i,:).*v(k,:))'./(w(:).^2));
                    end
                end
                fpvar = [covar(1,1) covar(2,2) covar(1,2) covar(2,1)];
            otherwise
                error('invalid method for fitting with fraction of probe');
        end
        
        % fit experimental coefficients solved above with theoretical
        % prediction from Kp and fraction of phase (ie lever rule)
        if strcmp(varargin{3},'grid')
            kpv = [[0.1:0.01:1] [1.1:0.1:10]]; % kp = [0.1:0.01:10];
            for k = 1:length(kpv)
                E = fpafpb(:,2);
                T = (kpv(k).*f_beta)./((1-f_beta) + kpv(k)*f_beta);      
                chisquare(k) = (E-T)'*(E-T);
            end   
            % find best chisq and corresponding Kp
            [chi2,i] = min(chisquare);
            Kp = kpv(i);
            clear kpv;
        else 
            options = optimset('LargeScale','on','Display','off');
            Kp = lsqcurvefit(@lincomb_fit_fun,1.0,f_beta,fpafpb(:,2),0.1,10,options,'fpafpb');
        end
        efpbKp = feval(@lincomb_fit_fun,Kp,f_beta,'fpafpb');
        efpaKp = 1-efpbKp;
        efpafpbKp = [efpaKp efpbKp];
%         chi2fpKp = (1/2).*sum([((fpafpb(:,1)-efpaKp).^2)./fpvar(1) ((fpafpb(:,2)-efpbKp).^2)./fpvar(2)],2);
% %             Chi2fp = (1/nc)*sum([((fpafpb(:,1)-efpa).^2)./fpvar(1);((fpafpb(:,2)-efpb).^2)./fpvar(2)]);
%         Chi2fpKp = mean(chi2fpKp(2:end-1));
        normfpKp = norm(fpafpb(:)-efpafpbKp(:));
        normfp2Kp = sum((fpafpb-efpafpbKp).^2,2);
        
        avg_kp = mean(kp(2:end-1));
        stdev_kp = std(kp(2:end-1));
        efpbkp = feval(@lincomb_fit_fun,avg_kp,f_beta,'fpafpb');
        efpakp = 1-efpbkp;
        efpafpbkp = [efpakp efpbkp];
%         chi2fpkp = (1/2).*sum([((fpafpb(:,1)-efpakp).^2)./fpvar(1) ((fpafpb(:,2)-efpbkp).^2)./fpvar(2)],2);
% %             Chi2fp = (1/nc)*sum([((fpafpb(:,1)-efpa).^2)./fpvar(1);((fpafpb(:,2)-efpb).^2)./fpvar(2)]);
%         Chi2fpkp = mean(chi2fpkp(2:end-1));
        normfpkp = norm(fpafpb(:)-efpafpbkp(:));
        normfp2kp = sum((fpafpb-efpafpbkp).^2,2);
        
        corrcoeff = fpvar(3)./(sqrt(fpvar(1).*fpvar(2)));
        
        lcSKp(:,1:2:2*ns) = B;
        lcSkp(:,1:2:2*ns) = B;
        lcS(:,1:2:2*ns) = B;
        for s = 1:ns
            E = I(:,s);
%             T = ((1-f_beta(s))/((1-f_beta(s)) + Kp*f_beta(s))).*B_alpha +
%             ((Kp*f_beta(s))/((1-f_beta(s)) + Kp*f_beta(s))).*B_beta;
            TeKp = Basis*efpafpbKp(s,:)';
            lcSKp(:,2*s) = TeKp;
            norm2SeKp(s) = (E-TeKp)'*(E-TeKp);
            
            Tekp = Basis*efpafpbkp(s,:)';
            lcSkp(:,2*s) = Tekp;
            norm2Sekp(s) = (E-Tekp)'*(E-Tekp);
            
            To = Basis*fpafpb(s,:)';
            lcS(:,2*s) = To;
            norm2So(s) = (E-To)'*(E-To);
            
            norm2SoSeKp(s) = (To-TeKp)'*(To-TeKp);
            norm2SoSekp(s) = (To-Tekp)'*(To-Tekp);
        end
        
%         chi2speckp = chi2speckp./(nb.*var_bline);
%         chisqkp = chi2speckp.*chi2fpkp';
%         Chisqkp = mean(chi2speckp(2:end-1))*Chi2fpkp;
%         
%         chi2specKp = chi2specKp./(nb.*var_bline);
%         chisqKp = chi2specKp.*chi2fpKp';
%         ChisqKp = mean(chi2specKp(2:end-1))*Chi2fpKp;
%         
%         chi2spec = chi2spec./(nb.*var_bline);
%         Chisqspec = mean(chi2spec(2:end-1));
        
%         fit_out = struct('Chisqspec',Chisqspec','chi2spec',chi2spec,...
%             'ChisqKp',ChisqKp,'chisqKp',chisqKp,'chi2specKp',chi2specKp,'Chi2fpKp',Chi2fpKp,...
%             'Chisqkp',Chisqkp,'chisqkp',chisqkp,'chi2speckp',chi2speckp,'Chi2fpkp',Chi2fpkp,...
%             'Kp',Kp,'kps',kp,'avg_kp',avg_kp,'stdev_kp',stdev_kp,'fpafpb',fpafpb,'fpvar',fpvar,...
%             'efpafpbKp',efpafpbKp,'efpafpbkp',efpafpbkp,'normfpKp',normfpKp,'normfpkp',normfpkp,...
%             'exp_S',S,'lincomb_S',lcS,'lincomb_SKp',lcSKp,'lincomb_Skp',lcSkp);
        
        X2So = norm2So./(nb.*var_bline);
        X2 = mean(X2So(2:end-1));
        
        X2Sekp = norm2Sekp./(nb.*var_bline);
        X2kp = mean(X2Sekp(2:end-1));
        
        X2SeKp = norm2SeKp./(nb.*var_bline);
        X2Kp = mean(X2SeKp(2:end-1));
        
        X2SoSeKp = norm2SoSeKp./(nb.*var_bline);
        X2oeKp = mean(X2SoSeKp(2:end-1));
        
        X2SoSekp = norm2SoSekp./(nb.*var_bline);
        X2oekp = mean(X2SoSekp(2:end-1));

        fit_out = struct('X2',X2,'X2So',X2So,'norm2So',norm2So,...
            'X2Kp',X2Kp,'X2SeKp',X2SeKp,'norm2SeKp',norm2SeKp,...
            'X2kp',X2kp,'X2Sekp',X2Sekp,'norm2Sekp',norm2Sekp,...
            'X2oeKp',X2oeKp,'X2SoSeKp',X2SoSeKp,'norm2SoSeKp',norm2SoSeKp,...
            'X2oekp',X2oekp,'X2SoSekp',X2SoSekp,'norm2SoSekp',norm2SoSekp,...
            'normfpKp',normfpKp,'normfp2Kp',normfp2Kp,'normfpkp',normfpkp,'normfp2kp',normfp2kp,...
            'Kp',Kp,'kps',kp,'avg_kp',avg_kp,'stdev_kp',stdev_kp,'fpafpb',fpafpb,'fpvar',fpvar,...
            'efpafpbKp',efpafpbKp,'efpafpbkp',efpafpbkp,'corrcoeff',corrcoeff,...
            'exp_S',S,'lincomb_S',lcS,'lincomb_SKp',lcSKp,'lincomb_Skp',lcSkp);
        
%         fit_out = struct('X2',X2,'X2So',X2So,'norm2So',norm2So,...
%             'X2Kp',X2Kp,'X2SeKp',X2SeKp,'norm2SeKp',norm2SeKp,...
%             'X2kp',X2kp,'X2Sekp',X2Sekp,'norm2Sekp',norm2Sekp,...
%             'X2oeKp',X2oeKp,'X2SoSeKp',X2SoSeKp,'norm2SoSeKp',norm2SoSeKp,...
%             'X2oekp',X2oekp,'X2SoSekp',X2SoSekp,'norm2SoSekp',norm2SoSekp,...
%             'Chi2fpKp',Chi2fpKp,'normfpKp',normfpKp,'Chi2fpkp',Chi2fpkp,'normfpkp',normfpkp,...
%             'Kp',Kp,'kps',kp,'avg_kp',avg_kp,'stdev_kp',stdev_kp,'fpafpb',fpafpb,'fpvar',fpvar,...
%             'efpafpbKp',efpafpbKp,'efpafpbkp',efpafpbkp,'corrcoeff',corrcoeff,...
%             'exp_S',S,'lincomb_S',lcS,'lincomb_SKp',lcSKp,'lincomb_Skp',lcSkp);
        
        disp(sprintf('done!'));
    case 'spectrum'
        lcS(:,1:2:2*ns) = B;
        switch varargin{2}
            case 'lsqcurvefit'
                options = optimset('display','off','jacobian','off');
                for s = 1:ns
                    f_beta = norm(C(s,:)-C_alpha)/norm(C_beta-C_alpha);
                    [kp(s),chisq(s)] = lsqcurvefit(@lincomb_fit_fun,1.0,f_beta,I(:,s),0.1,10,options,'spectrum','lsqcurvefit',Basis); 
                    chisq(s) = (1/(nb*var_bline(s)))*chisq(s);
                    fpa = ((1-f_beta)/((1-f_beta) + kp(s)*f_beta));
                    fpb = ((kp(s)*f_beta)/((1-f_beta) + kp(s)*f_beta));
                    fpafpb(s,:) = [fpa fpb];
                    lcS(:,2*s) = Basis*fpafpb(s,:)';
                end
            case 'chisq_min'
                options = optimset('display','off');
                for s = 1:ns
                    f_beta = norm(C(s,:)-C_alpha)/norm(C_beta-C_alpha);
                    [kp(s),chisq(s)] = fminbnd(@lincomb_fit_fun,0.1,10,options,f_beta,'spectrum','chisq_min',I(:,s),Basis,nb,std_bline(s));
                    fpa = ((1-f_beta)/((1-f_beta) + kp(s)*f_beta));
                    fpb = ((kp(s)*f_beta)/((1-f_beta) + kp(s)*f_beta));
                    fpafpb(s,:) = [fpa fpb];
                    lcS(:,2*s) = Basis*fpafpb(s,:)';
                end
            case 'grid'
                kpv = [[0.1:0.01:1] [1.1:0.1:10]];
                for s = 1:ns
                    E = I(:,s);
                    f_beta = norm(C(s,:)-C_alpha)/norm(C_beta-C_alpha);
                    for k = 1:length(kpv)
                        fpa(k) = ((1-f_beta)/((1-f_beta) + kpv(k)*f_beta));
                        fpb(k) = ((kpv(k)*f_beta)/((1-f_beta) + kpv(k)*f_beta));
                        T(:,k) = Basis*[fpa(k);fpb(k)];
%                         T = ((1-f_beta)/((1-f_beta) + kpv(k)*f_beta)).*B_alpha + ((kpv(k)*f_beta)/((1-f_beta) + kpv(k)*f_beta)).*B_beta;
                        chisquare(k) = (1/(nb*var_bline(s)))*((E-T(:,k))'*(E-T(:,k)));  
                    end
                    % find best chisq and corresponding kp for each
                    % spectrum
                    [chisq(s),i] = min(chisquare);
                    kp(s) = kpv(i);
                    fpafpb(s,:) = [fpa(i) fpb(i)];
                    lcS(:,2*s) = T(:,i);
                end   
            otherwise
                error('invalid algorithm for spectrum fitting');
        end 

        Kp = mean(kp(2:end-1));
        stdev_Kp = std(kp(2:end-1));
        Chisq = mean(chisq(2:end-1)); 

        fit_out = struct('Chisq',Chisq,'chisqs',chisq,'Kp',Kp,'kps',kp,'stdevkp',stdev_Kp,'fpafpb',fpafpb,'exp_S',S,'lincomb_S',lcS);
        disp(sprintf('done!'));
    case 'tieline' % just means kp = constant when fitting by linear combination
        lcS(:,1:2:2*ns) = B;
        switch varargin{2}
            case 'lsqcurvefit'
                options = optimset('display','off','jacobian','off');
                for s = 1:ns
                    f_beta(s,1) = norm(C(s,:)-C_alpha)/norm(C_beta-C_alpha);  
                end
                data = reshape(I,nb*ns,1);
                Kp = lsqcurvefit(@lincomb_fit_fun,1.0,f_beta,data,0.1,10,options,'tieline',Basis); 
                for s = 1:ns
                    E = I(:,s);
                    fpa = ((1-f_beta(s))/((1-f_beta(s)) + Kp*f_beta(s)));
                    fpb = 1-fpa;
                    fpafpb(s,:) = [fpa fpb];
                    if s == 1 || s == ns
                        kp(s) = 1;
                    else
                        kp(s) = (fpb*(1-f_beta(s)))/(fpa*f_beta(s));
                    end
                    T = Basis*fpafpb(s,:)';
                    lcS(:,2*s) = T;
                    chisq(s) = (1/(nb*var_bline(s)))*((E-T)'*(E-T));
                end
            case 'grid'
                for k = 1:length(kp)
                    for s = 1:nfs
                        f_beta(s) = norm(C(s,:)-C_alpha)/norm(C_beta-C_alpha);
                        E = S(:,s);
                        T = ((1-f_beta)/((1-f_beta) + kp(k)*f_beta)).*B_alpha + ((kp(k)*f_beta)/((1-f_beta) + kp(k)*f_beta)).*B_beta;
                        chisq(k,s) = (1/((nb-1)*(std_bline(s)^2)))*((E-T)'*(E-T));
                        Kp(k,s) = kp(k);
                    end
                end

                kpv = [[0.1:0.01:1] [1.1:0.1:10]];
                for s = 1:ns
                    E = I(:,s);
                    f_beta = norm(C(s,:)-C_alpha)/norm(C_beta-C_alpha);
                    for k = 1:length(kpv)
                        fpa(k) = ((1-f_beta)/((1-f_beta) + kpv(k)*f_beta));
                        fpb(k) = ((kpv(k)*f_beta)/((1-f_beta) + kpv(k)*f_beta));
                        T(:,k) = Basis*[fpa(k);fpb(k)];
        %                         T = ((1-f_beta)/((1-f_beta) + kpv(k)*f_beta)).*B_alpha + ((kpv(k)*f_beta)/((1-f_beta) + kpv(k)*f_beta)).*B_beta;
                        chisquare(k) = (1/(nb*var_bline(s)))*((E-T)'*(E-T));  
                    end
                    % find best chisq and corresponding kp for each
                    % spectrum
                    [chisq(s),i] = min(chisquare);
                    kp(s) = kpv(i);
                    fpafpb(s,:) = [fpa(i) fpb(i)];
                    lcS(:,2*s) = T(:,i);
                end   

                % find best chisq and corresponding Kp
                avg_chisq = mean(chisq,2);
                [minchisq,imin] = min(avg_chisq);
                minKp = Kp(imin);
            otherwise
                error('invalid algorithm for tieline spectrum fitting');
        end
        
        Chisq = mean(chisq(2:end-1));
        
        fit_out = struct('Chisq',Chisq,'chisqs',chisq,'Kp',Kp,'kps',kp,'fpafpb',fpafpb,'exp_S',S,'lincomb_S',lcS);
        disp(sprintf('done!'));
    otherwise
        error('invalid data to be fit, choices are either "fpafpb", "spectrum", or "tie_line"');
end 

return