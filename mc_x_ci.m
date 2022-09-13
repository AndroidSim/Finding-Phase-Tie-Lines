function output = mc_x_ci(fit_output,C_coex,S_coex,C_bdy,S_bdy,mcmethod,varargin)
% output = mc_x_ci(fit_output,C_coex,S_coex,C_bdy,S_bdy,varargin)

[ncC,ncd] = size(C_coex);
[ncB,ncScol] = size(S_coex);
[nbC,nbd] = size(C_bdy);
[nbB,nbScol] = size(S_bdy);
ncS = ncScol/2;
nbS = nbScol/2;

allS = [S_coex S_bdy];
allS = align_spectra(allS); 
allS = normalize_spectra(allS,'d','area',1000);
Scoex = allS(:,1:ncScol);
Sbdy = allS(:,ncScol+1:ncScol+nbScol);
clear allS;

if ~any(strcmp(mcmethod,{'bootstrap';'sampling'}))
    error('monte carlo method either bootstrap or sampling');
end

switch mcmethod
    case 'bootstrap'
        if isempty(varargin)
            N = 50;
        else
            if isscalar(varargin{1})
                N = varargin{1};
            else
                error('1st variable argument must be a scalar specifying the number of monte carlo iterations');
            end
        end
    case 'sampling'
        if isempty(varargin)
            N = 1000;
            p = 0.95;
        else
            if isscalar(varargin{1})
                N = varargin{1};
            else
                error('1st variable argument must be a scalar specifying the number of monte carlo iterations');
            end
            if isscalar(varargin{2})
                p = varargin{2};
            else
                error('2nd variable argument must be a scalar specifying the confidence level (probability)');
            end
        end
    otherwise
        error('invalid monte carlo method');
end

% initialization
bdyconfig = fit_output.configs.bdyconfig;
tlconfig = fit_output.configs.tlconfig;
Kpconfig = fit_output.configs.Kpconfig;

switch tlconfig
    case 'parallel'
        switch bdyconfig
            case '2critpts'
                switch Kpconfig
                    case 'constant'
                        % x = parameters, length(x) should = 2 = [slope Kp]
                           
                    case 'function'
                         % x = parameters, length(x) should = 2 = [slope c]
                        fit_measure = fit_output.fit.fit_measure;
                        X2TD = fit_output.fit.X2_TD;
                        x = fit_output.fit.x;
                        df = length(x);
                        isaved = 1;
                        for i = 1:N
                            trial_x = x' + [unifrnd(-0.1,0.1);unifrnd(-1,1)];
                            trial_output = simulate_tielines(trial_x,C_coex,S_coex,C_bdy,S_bdy,bdyconfig,tlconfig,Kpconfig);
                            trial_fit_measure = trial_output.fit.fit_measure;
                            if trial_fit_measure - fit_measure < chi2inv(p,df)
                                xsaved(isaved,:) = trial_output.fit.x;
                                fitmsaved(isaved,1) = trial_fit_measure;
                                isaved = isaved+1;
                                pack;
                            end
                        end
                    otherwise
                        error('invalid Kp configuration');
                end
            case '1critpt/1endtl' % parallel tielines
                switch Kpconfig
                    case 'constant'
                        % x = parameters, length(x) should = 3 = [ep1 ep2 Kp]
                        
                    case 'function'
                        % x = parameters, length(x) should = 4 = [ep1 ep2 c1 c2]
                        fit_measure = fit_output.fit.fit_measure;
                        X2TD = fit_output.fit.X2_TD;
                        x = fit_output.fit.x;
                        df = length(x);
                        xsaved = zeros(1,4);
                        fitmsaved = 0;
                        isaved = 1;
                        for i = 1:N
                            trial_x = x' + [unifrnd(-0.05,0.05);unifrnd(-0.05,0.05);unifrnd(-1,1);unifrnd(-1,1)];
                            trial_output = simulate_tielines(trial_x,C_coex,S_coex,C_bdy,S_bdy,bdyconfig,tlconfig,Kpconfig);
                            trial_fit_measure = trial_output.fit.fit_measure;
                            if trial_fit_measure - fit_measure < chi2inv(p,df)
                                xsaved(isaved,:) = trial_output.fit.x;
                                fitmsaved(isaved,1) = trial_fit_measure;
                                isaved = isaved+1;
                                pack;
                            end
                        end
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
                    case 'function'
                        % x = parameters, length(x) should = 3 = [cp1 cp2 c]
                    otherwise
                        error('invalid Kp configuration');
                end
            case '1critpt/1endtl' % tangent tielines
                switch Kpconfig
                    case 'constant'
                        % x = parameters, length(x) should = 4 = [cp ep1 ep2 Kp]
                    case 'function'
                        % x = parameters, length(x) should = 5 = [cp ep1 ep2 c1 c2]
                        fit_measure = fit_output.fit.fit_measure;
                        X2TD = fit_output.fit.X2_TD;
                        x = fit_output.fit.x;
                        df = length(x);
                        xsaved = zeros(1,5);
                        fitmsaved = 0;
                        isaved = 1;
                        for i = 1:N
                            trial_x = x' + [unifrnd(-0.05,0.05);unifrnd(-0.05,0.05);unifrnd(-0.05,0.05);unifrnd(-1,1);unifrnd(-1,1)];
                            trial_output = simulate_tielines(trial_x,C_coex,S_coex,C_bdy,S_bdy,bdyconfig,tlconfig,Kpconfig);
                            trial_fit_measure = trial_output.fit.fit_measure;
                            if trial_fit_measure - fit_measure < chi2inv(p,df)
                                xsaved(isaved,:) = trial_output.fit.x;
                                fitmsaved(isaved,1) = trial_fit_measure;
                                isaved = isaved+1;
                                pack;
                            end
                        end
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
                        
                    case 'function'
                        % x = parameters, length(x) should = 4 = [cp1 cp2 gamma c]
                        
                    otherwise
                        error('invalid Kp configuration');
                end
            case '1critpt/1endtl' % ruled tielines
                switch Kpconfig
                    case 'constant'
                        % x = parameters, length(x) should = 5 = [cp ep1 ep2 gamma Kp]
                        
                    case 'function'
                        % x = parameters, length(x) should = 6 = [cp ep1 ep2 gamma c1 c2]
                        % xmin = [0.3458  0.6974  0.8555  2.3000  4.2630  -3.8680]
%                         ternary_plot(Ccoex,'.k','markersize',20)
%                         ternary_plot(Ccoex+Ccoex.*repmat(unifrnd(-0.02,0.02,51,1),1,3),'.r','markersize',20)
                        switch mcmethod
                            case 'bootstrap'
                                B = repmat(fit_output.spectra.magnetic_fields,1,ncS);
                                ncB = length(fit_output.spectra.magnetic_fields);
                                tScoex = zeros(ncB,ncS*2);
                                tScoex(:,1:2:end) = B;
                                tScoex(:,2:2:end) = fit_output.spectra.theory_cS;
                                oScoex = zeros(ncB,ncS*2);
                                oScoex(:,1:2:end) = B;
                                oScoex(:,2:2:end) = fit_output.spectra.obs_cS;
                                [tr,ts] = find(isnan(tScoex));
                                if ~isempty(ts)
                                    tScoex(tr,ts) = Scoex(tr,ts);
                                end
                                [or,os] = find(isnan(oScoex));
                                if ~isempty(os)
                                    oScoex(or,os) = Scoex(or,os);
                                end
                                if ~all(Scoex(:) == oScoex(:))
                                    error('Scoex not equal to oScoex');
                                end
                                fit_conditions.bdy_config = bdyconfig;
                                fit_conditions.tieline_config = tlconfig;
                                fit_conditions.critptparameter = [0.33 0.38];
                                fit_conditions.endpt1parameter = [0.67 0.72];
                                fit_conditions.endpt2parameter = [0.82 0.86];
                                xsaved = zeros(1,6);
                                fitmsaved = 0;
                                X2TDsaved = 0;
                                for i = 1:N
                                    synScoex = tScoex;
                                    synScoex(:,2:2:end) = tScoex(:,2:2:end)+sqrt(repmat(var(oScoex([1:200 end-200:end]',2:2:end)),ncB,1)).*randn(ncB,ncS);
%                                     synScoex = S_coex;
%                                     synScoex(:,2:2:end) = S_coex(:,2:2:end)+sqrt(repmat(var(S_coex([1:200 end-200:end]',2:2:end)),ncB,1)).*randn(ncB,ncS);
%                                     synSbdy = S_bdy;
%                                     synSbdy(:,2:2:end) = S_bdy(:,2:2:end)+sqrt(repmat(var(S_bdy([1:200 end-200:end]',2:2:end)),nbB,1)).*randn(nbB,nbS);
%                                     trial_output = tieline_field_fitting(synScoex,C_coex,synSbdy,C_bdy,'simplex',fit_conditions);
%                                     synCcoex = C_coex+C_coex.*repmat(unifrnd(-0.02,0.02,ncC,1),1,ncd);
                                    synCcoex = C_coex+C_coex.*[unifrnd(-0.02,0.02,ncC,1) unifrnd(-0.02,0.02,ncC,1) unifrnd(-0.02,0.02,ncC,1)];
%                                     trial_output = tieline_field_fitting(S_coex,synCcoex,S_bdy,C_bdy,'simplex',fit_conditions);
                                    trial_output = tieline_field_fitting(synScoex,synCcoex,Sbdy,C_bdy,'simplex',fit_conditions);
                                    xsaved(i,:) = trial_output.fit.x;
                                    fitmsaved(i,1) = trial_output.fit.fit_measure;
                                    X2TDsaved(i,1) = trial_output.fit.X2_TD;
                                    pack;
                                end
                            case 'sampling'
                                fit_measure = fit_output.fit.fit_measure;
                                X2TD = fit_output.fit.X2_TD;
                                x = fit_output.fit.x;
                                df = length(x);
                                xsaved = zeros(1,6);
                                fitmsaved = 0;
                                isaved = 1;
                                for i = 1:N
                                    trial_x = x' + [unifrnd(-0.1,0.04);unifrnd(-0.03,0.03);unifrnd(-0.04,0.01);unifrnd(-2,2);unifrnd(-2,2);unifrnd(-2,2)];
                                    trial_output = simulate_tielines(trial_x,C_coex,S_coex,C_bdy,S_bdy,bdyconfig,tlconfig,Kpconfig);
                                    trial_X2TD = trial_output.fit.X2_TD;
                                    if trial_X2TD - X2TD < chi2inv(p,df)
                                        xsaved(isaved,:) = trial_output.fit.x;
                                        fitmsaved(isaved,1) = trial_output.fit.fit_measure;;
                                        X2TDsaved = trial_X2TD;
                                        isaved = isaved+1;
                                        pack;
                                    end
                                end
                            otherwise
                                error('invalid monte carlo method');
                        end 
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
output.xsaved = xsaved;
output.fitmsaved = fitmsaved;
output.X2TDsaved = X2TDsaved;
output.avgx = mean(xsaved);
output.stdx = std(xsaved);

return