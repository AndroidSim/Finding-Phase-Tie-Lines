function output = simulate_tielines(x,C_coex,S_coex,C_bdy,S_bdy,bdyconfig,tlconfig,Kpconfig)

% checking coexistence spectra
if all(size(S_coex) > 1) & ndims(S_coex) == 2 % if spectra is a matrix
    [ncB,ncScol] = size(S_coex);

    if ncScol < 2 || rem(ncScol,2) ~= 0
        error('each coexistence region spectrum is two columns: [B-field intensity_values]');
    end
    
    ncS = ncScol/2;
else
    error('first argument must be a matrix of the coexistence region spectra');
end

% checking coexistence region compositions
if  any(size(C_coex) == 1) & any(size(C_coex) > 1) % if compositions is a vector
    [nrow,ncol] = size(C_coex);
    
    if ncol == 1 % if column vector
        C_coex = C_coex'; % convert to row vector
    end
    
    ncC = 1; 
    cCd = length(C_coex);
    clear nrow ncol;
elseif all(size(C_coex) > 1) & ndims(C_coex) == 2 % if compositions is a matrix 
    [ncC,cCd] = size(C_coex);
end

if ncC ~= ncS
    error('number of coexistence region compositions must equal number of coexistence region spectra');
end

% sum of coexistence region mole fractions must equal one
if all(sum(C_coex')) ~= 1
    error('sum of coexistence region mole fractions must equal one');
end

% checking boundary spectra
if all(size(S_bdy) > 1) & ndims(S_bdy) == 2 % if spectra is a matrix
    [nbB,nbScol] = size(S_bdy);

    if nbScol < 2 || rem(nbScol,2) ~= 0
        error('each boundary spectrum is two columns: [B-field intensity_values]');
    end
    
    nbS = nbScol/2;
else
    error('third argument must be a matrix of the boundary spectra');
end

% checking boundary compositions
if  any(size(C_bdy) == 1) & any(size(C_bdy) > 1) % if compositions is a vector
    [nrow,ncol] = size(C_bdy);
    
    if ncol == 1 % if column vector
        C_bdy = C_bdy'; % convert to row vector
    end
    
    nbC = 1; 
    bCd = length(C_bdy);
    clear nrow ncol; 
elseif all(size(C_bdy) > 1) & ndims(C_bdy) == 2 % if compositions is a matrix 
    [nbC,bCd] = size(C_bdy);
end

if nbC ~= nbS
    error('number of boundary compositions must equal number of boundary spectra');
end
    
% sum of boundary mole fractions must equal one
if all(sum(C_bdy')) ~= 1
    error('sum of boundary mole fractions with spectra must equal one');
end

allS = [S_coex S_bdy];
allS = align_spectra(allS); 
allS = normalize_spectra(allS,'d','area',1000);
S_coex = allS(:,1:ncScol);
S_bdy = allS(:,ncScol+1:ncScol+nbScol);
clear allS;
magnetic_fields = S_coex(:,1); % or S_bdy, magnetic field values the same for all spectra because of alignment
% reduced all spectra to a vector of intensity values (magnetic fields not needed for fitting)
S_coex = S_coex(:,2:2:end);
S_bdy = S_bdy(:,2:2:end);

% linearly transform ternary(i.e. homogeneous barycentric coordinates) to 
% cartesian coordinates using 1st (Sm) and 3rd (chol) mole fractions
C_bdy = tern2cart(C_bdy,1);
C_coex = tern2cart(C_coex,1);

% determine sizes of inputs before fitting
[ncC,ncd] = size(C_coex);
[ncB,ncS] = size(S_coex);
[nbC,nbd] = size(C_bdy);
[nbB,nbS] = size(S_bdy);

% [u,w,v] = svd([S_coex S_bdy],0);
% S = u(:,1:5)*w(1:5,1:5)*v(:,1:5)';
% S_coex = S(:,1:ncS);
% S_bdy = S(:,ncS+1:end);

for s = 1:ncS
    cSbaseline(:,s) = S_coex([1:200 end-200:end]',s);
    avg_cSbline(s) = mean(cSbaseline(:,s));
    var_cSbline(s) = var(cSbaseline(:,s));
    std_cSbline(s) = std(cSbaseline(:,s));
end 
for s = 1:nbS
    bSbaseline(:,s) = S_bdy([1:200 end-200:end]',s);
    avg_bSbline(s) = mean(bSbaseline(:,s));
    var_bSbline(s) = var(bSbaseline(:,s));
    std_bSbline(s) = std(bSbaseline(:,s));
end 
stdbS = std(S_bdy')';
varbS = var(S_bdy')';

% weight_fxn_S = ones(size(stdbS));

weight_fxn_S = zeros(size(stdbS));
weight_fxn_S(300:ncB-200) = ones(length(300:ncB-200),1);

w = zeros(size(stdbS));
w(stdbS >= 10*mean(std_bSbline)) = 1;
weight_fxn_fp = w;
% weight_fxn_S = w;

% w = zeros(size(stdbS));
% w(weight_fxn >= 10*mean(std_bSbline)) = 1;
% w(weight_fxn >= 100*mean(std_bSbline)) = 2;
% w(weight_fxn >= 1000*mean(std_bSbline)) = 3;
% w(weight_fxn >= 10000*mean(std_bSbline)) = 4;
% weight_fxn = w;

spectra = struct('magnetic_fields',magnetic_fields,'cSbline_var',var_cSbline,'cSbline_std',std_cSbline,...
    'bSC_std',stdbS,'bSC_var',varbS,'weight_fxn_S',weight_fxn_S,'weight_fxn_fp',weight_fxn_fp);
boundary.config = bdyconfig;
part_coeff.config = Kpconfig;
[bdy_length,b] = chord_length(C_bdy);
boundary.parameterization = b;

switch tlconfig
    case 'parallel'
        switch bdyconfig
            case '2critpts'
                switch Kpconfig
                    case 'constant'
                        % x = parameters, length(x) should = 2 = [slope Kp]
                        if length(x) ~= 2
                            error('need 2 parameters, [slope Kp]');
                        end
                        [fitmeasure,output] = parallel_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,[]);
                    case 'function'
                        % x = parameters, length(x) should = 2 = [slope c]
                        if length(x) ~= 2
                            error('need 2 parameters, [slope c]');
                        end
                        [fitmeasure,output] = parallel_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,[]);
                    otherwise
                        error('invalid Kp configuration');
                end
            case '1critpt/1endtl' % parallel tielines
                switch Kpconfig
                    case 'constant'
                        % x = parameters, length(x) should = 3 = [ep1 ep2 Kp]
                        if length(x) ~= 3
                            error('need 3 parameters, [ep1 ep2 Kp]');
                        end
                        % analyze boundary before calculating expected coexistence spectra
                        % slope calculated from the end points of the end
                        % tieline
                        [fitmeasure,output] = parallel_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,[]);
                    case 'function'
                        % x = parameters, length(x) should = 4 = [ep1 ep2 c1 c2]
                        if length(x) ~= 4
                            error('need 4 parameters, [ep1 ep2 c1 c2]');
                        end
                        % analyze boundary before calculating expected coexistence spectra
                        % slope calculated from the end points of the end
                        % tieline
                        [fitmeasure,output] = parallel_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,[]);
                    otherwise
                        error('invalid Kp configuration');
                end
            case '2endtls' % parallel tielines
                switch Kpconfig
                    case 'constant'
                    case 'function'
                    otherwise
                        error('invalid Kp configuration');
                end
            otherwise
                error('invalid boundary configuration');
        end      
    case 'tangent'
        switch bdyconfig
            case '2critpts'
                switch Kpconfig
                    case 'constant'
                        % x = parameters, length(x) should = 3 = [cp1 cp2 Kp]
                        if length(x) ~= 3
                            error('need 3 parameters, [cp1 cp2 Kp]');
                        end
                        [fitmeasure,output] = tangent_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,[]);
                    case 'function'
                        % x = parameters, length(x) should = 3 = [cp1 cp2 c]
                        if length(x) ~= 3
                            error('need 3 parameters, [cp1 cp2 c]');
                        end
                        [fitmeasure,output] = tangent_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,[]);
                    otherwise
                        error('invalid Kp configuration');
                end
            case '1critpt/1endtl' % tangent tielines
                switch Kpconfig
                    case 'constant'
                        % x = parameters, length(x) should = 4 = [cp ep1 ep2 Kp]
                        if length(x) ~= 4
                            error('need 4 parameters, [cp ep1 ep2 Kp]');
                        end
                        % analyze boundary before calculating expected coexistence spectra
                        [fitmeasure,output] = tangent_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,[]);
                    case 'function'
                        % x = parameters, length(x) should = 5 = [cp ep1 ep2 c1 c2]
                        if length(x) ~= 5
                            error('need 5 parameters, [cp ep1 ep2 c1 c2]');
                        end
                        % analyze boundary before calculating expected coexistence spectra
                        [fitmeasure,output] = tangent_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,[]);
                    otherwise
                        error('invalid Kp configuration');
                end
            case '2endtls' % tangent tielines
                switch Kpconfig
                    case 'constant'
                    case 'function'
                    otherwise
                        error('invalid Kp configuration');
                end
            otherwise
                error('invalid boundary configuration');
        end
    case 'ruled'
        switch bdyconfig
            case '2critpts'
                switch Kpconfig
                    case 'constant'
                        % x = parameters, length(x) should = 4 = [cp1 cp2 gamma Kp]
                        if length(x) ~= 4
                            error('need 4 parameters, [cp1 cp2 gamma Kp]');
                        end
                        [fitmeasure,output] = ruled_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,[]);
                    case 'function'
                        % x = parameters, length(x) should = 4 = [cp1 cp2 gamma c]
                        if length(x) ~= 4
                            error('need 4 parameters, [cp1 cp2 gamma c]');
                        end
                        [fitmeasure,output] = ruled_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,[]);
                    otherwise
                        error('invalid Kp configuration');
                end
            case '1critpt/1endtl' % ruled tielines
                switch Kpconfig
                    case 'constant'
                        % x = parameters, length(x) should = 5 = [cp ep1 ep2 gamma Kp]
                        if length(x) ~= 5
                            error('need 5 parameters, [cp ep1 ep2 gamma Kp]');
                        end
                        % analyze boundary before calculating expected coexistence spectra
                        [fitmeasure,output] = ruled_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,[]);
                    case 'function'
                        % x = parameters, length(x) should = 6 = [cp ep1 ep2 gamma c1 c2]
                        if length(x) ~= 6
                            error('need 6 parameters, [cp ep1 ep2 gamma c1 c2]');
                        end
                        % analyze boundary before calculating expected coexistence spectra
                        [fitmeasure,output] = ruled_tieline_fit_fxn(x,C_coex,S_coex,C_bdy,S_bdy,boundary,part_coeff,spectra,[]);
                    otherwise
                        error('invalid Kp configuration');
                end
            case '2endtls' % ruled tielines
                switch Kpconfig
                    case 'constant'
                    case 'function'
                    otherwise
                        error('invalid Kp configuration');
                end
            otherwise
                error('invalid boundary configuration');
        end
    otherwise 
        error('invalid tieline configuration');
end

return