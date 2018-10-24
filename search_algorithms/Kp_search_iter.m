function [Kps]=Kp_search_iter(exp_surface,boundary_pts,ntiepts,ntielines,xdonr,xaceptr,max_kpd,max_kpa,tie_type,expconsta,expconstb,threshold)
% the function Kp_search performs an interative search to find the best Kp-vs-tie line function to be used to fit
% or simulate experimental FRET data surfaces

options=optimset('MaxFunEvals',1000,'TolFun',1e-10,'TolX',1e-10,'MaxIter',1000,'LargeScale','off','Display','off');
% create lookup table for Kps
disp('calculating lookup table for Kps');
tl_kpd=[1.0:0.2:max_kpd];
tl_kpa=[1.0:0.2:max_kpa];
tl_maxf=zeros(length(tl_kpd),length(tl_kpa));
tl_maxa=zeros(length(tl_kpd),length(tl_kpa));
nelem=1;
for d=1:length(tl_kpd)
    for a=1:length(tl_kpa)
        [maxa,maxf,exitflag]=fmincon(@FRET_eqn,0.5,[],[],[],[],0,1,[],options,-1,-1,xdonr,xaceptr,tl_kpd(d),tl_kpa(a),expconsta,expconstb);
        if (exitflag < 0)
            error('for lookup table, the maximization did not converge')
        end
        tl_maxf(d,a)=-maxf;
        tl_maxa(d,a)=maxa;
        lookup_table(nelem,1)=maxa;
        lookup_table(nelem,2)=-maxf;
        lookup_table(nelem,3)=d;
        lookup_table(nelem,4)=a;
        nelem=nelem+1;
    end
end

% for starting iteration, calculate best fit tie lines and best fit surface with Kp fxn that is 
% parabolic with the maximum kp at middle of tieline for both donor and acceptor
tie_num=[0:ntielines+1];
pd=polyfit([0 median(tie_num) tie_num(end)],[1 max_kpd 1],2); 
pa=polyfit([0 median(tie_num) tie_num(end)],[1 max_kpa 1],2); 
kpd=polyval(pd,tie_num');
kpa=polyval(pa,tie_num');
kpd(1)=[];
kpd(end)=[];
kpa(1)=[];
kpa(end)=[];

% iterate until the total rms between before and after donor and acceptor Kp vectors is below some threshold
niter=0;
rms=Inf;
ifit_value=100;
while (rms >= threshold(1) | dfit_value <= threshold(2))
    disp('using donor and acceptor Kp fxns to find best fit tielines and surface');
    % run fitting routine to get best fit tielines and surface
    [p2_output,data_surface]=fit_surface(exp_surface,boundary_pts,ntiepts,ntielines,xdonr,xaceptr,kpd,kpa,tie_type,0,expconsta,expconstb);
    best_tielines=p2_output(5).data;
    best_surface=p2_output(6).data;
    lowest=p2_output(8).data;
    best_Kps=p2_output(12).data{lowest(1),lowest(2)};
    best_ab_pts=p2_output(4).data{lowest(1),lowest(2)};
    fit_value=p2_output(9).data(lowest(1),lowest(2));
    disp(sprintf('fit value = %f \n',fit_value));
    
    % loop over each tie line in best fit tie lines and determine new Kps from experimental FRET along each
    % tie line
    disp('calculating new donor and acceptor Kp fxns');
    [ntie,ntpts,xy]=size(best_tielines);
    for i=1:ntie
        tie_line=best_tielines(i,:,:);
        tie_line=squeeze(tie_line);
        [alpha,beta]=alphas_betas(tie_line);
        alpha=alpha';
        beta=beta';
        % get FRET data along tie line
        FRET_data=interp2(exp_surface',tie_line(:,1),tie_line(:,2));
        % scale experimental FRET data along each tie line of best fit tie lines
        if (isequal(best_ab_pts(i,1,1),tie_line(1,1)) & isequal(best_ab_pts(i,2,1),tie_line(1,2))) % x and y
            if (isequal(best_ab_pts(i,1,2),tie_line(end,1)) & isequal(best_ab_pts(i,2,2),tie_line(end,2)))
                Fa=FRET_data(1);
                Fb=FRET_data(end);
                % calculte base line FRET
                base_line=feval(@FRET_eqn,alpha,Fa,Fb,xdonr,xaceptr,1,1,expconsta,expconstb);
                FRET_scaled=FRET_scaling(FRET_data,alpha,Fa,Fb,xdonr,xaceptr,best_Kps(i,1),best_Kps(i,2),expconsta,expconstb);
            else
                error('alpha pt == first tie line pt, but beta pt ~= end tie line pt')
            end
        elseif (isequal(best_ab_pts(i,1,1),tie_line(end,1)) & isequal(best_ab_pts(i,2,1),tie_line(end,2)))
            if (isequal(best_ab_pts(i,1,2),tie_line(1,1)) & isequal(best_ab_pts(i,2,2),tie_line(1,2)))
                Fa=FRET_data(end);
                Fb=FRET_data(1);
                % calculte base line FRET
                base_line=feval(@FRET_eqn,beta,Fa,Fb,xdonr,xaceptr,1,1,expconsta,expconstb);
                FRET_scaled=FRET_scaling(FRET_data,beta,Fa,Fb,xdonr,xaceptr,best_Kps(i,1),best_Kps(i,2),expconsta,expconstb);
                FRET_scaled=flipud(FRET_scaled);
            else
                error('beta pt == first tie line pt, but alpha pt ~= end tie line pt')
            end
        end
        FRET_scaled=FRET_scaled-FRET_scaled(1); % or FRET_scaled(end), should be the same value, now FRET_scaled should start and end at 1
        % find maximum value of scaled FRET and the alpha at that maximum value
        % check to see if data has a maximum
        dFRET=gradient(FRET_scaled);
        for j=1:length(dFRET)-1
            data_pt=dFRET(j);
            next_data_pt=dFRET(j+1);
            if (isequal(sign(next_data_pt),-1) & isequal(sign(data_pt),1))
                maxyorn=1;
                break;
            elseif (isequal(sign(next_data_pt),0) & isequal(sign(data_pt),1))
                maxyorn=1;
                break;
            elseif (isequal(sign(data_pt),0) & isequal(sign(next_data_pt),-1))
                maxyorn=1;
                break;
            else
                % no maximum
                maxyorn=0;
            end
        end
        if (maxyorn == 1)
            maxf=max(FRET_scaled);
            maxa=alpha(find(maxf == FRET_scaled));
            if (length(maxa) > 2)
                error('number of alpha at max > 2')
            else
                if (isequal(length(maxa),2))
                    maxa=mean(maxa);
                elseif (isequal(length(maxa),1))
                    maxa=maxa;
                else
                    error('alpha of max of data is empty')
                end
            end
        else
            error('scaled FRET data has no maximum, as it should')
        end
        % lookup the Kp(donor) and Kp(acceptor) that would give maxf and maxa from lookup_table, given the assumption
        % that the Kp(donor) [fast dio] >= Kp(acceptor) [c18:1 dii]
        % maxf_matrix=reshape(lookup_table(:,2),length(tl_kpd),length(tl_kpa));
        % maxa_matrix=reshape(lookup_table(:,1),length(tl_kpd),length(tl_kpa));
        C_maxf=contourc(tl_maxf,[maxf maxf]); % top row ~= kpa, bottom row ~= kpd
        C_maxf(:,1)=[];
        C_maxf=C_maxf';
        C_maxa=contourc(tl_maxa,[maxa maxa]); % top row ~= kpa, bottom row ~= kpd
        C_maxa(:,1)=[];
        C_maxa=C_maxa';
        kps_maxf=[tl_kpd(round(C_maxf(:,2)))' tl_kpa(round(C_maxf(:,1)))'];
        kps_maxa=[tl_kpd(round(C_maxa(:,2)))' tl_kpa(round(C_maxa(:,1)))'];
        common_kps=intersect(kps_maxf,kps_maxa,'rows');
        nint=size(common_kps,1);
        if (isequal(nint,0))
            error('maxf and maxa have no common kps')
        elseif (nint > 1)
            for j=1:nint
                if (common_kps(j,1) >= common_kps(j,2)) % kp donor >= kp acceptor
                    kpd(i)=common_kps(j,1);
                    kpa(i)=common_kps(j,2);
                end
            end
        else
            kpd(i)=common_kps(nint,1);
            kpa(i)=common_kps(nint,2);
        end
        %for j=1:length(C)
        %    if (tl_kpd(round(C(1,j))) >= tl_kpa(round(C(2,j))))
        %        minxC(i)=fminsearch(@FRET_eqn,0.2,options,-1,-1,0.001,0.00025,C(1,i),C(2,i),238,508);
        %    end
        %end   
    end
    
    rmsd=mean((kpd-best_Kps(:,1)).^2);
    rmsa=mean((kpa-best_Kps(:,2)).^2);
    rms(niter)=rmsd+rmsa;
    dfit_value=abs(fit_value-ifit_value);
    ifit_value=fit_value;
    
    niter=niter+1;
    
    disp(sprintf('rms = %f \n',rms));
    disp(sprintf('niter = %f \n',niter));
end

Kps=cat(tie_num,cat(kpd,kpa,2),2);