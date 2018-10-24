function [output,data_surface]=Kp_fxns_search(exp_surf,boundary_pts,ntiepts,ntielines,tie_start,start_inc,tie_stop,top_inc,nbest2save,min_kpd,max_kpd,kpd_inc,min_kpa,max_kpa,kpa_inc,xdonr,xaceptr,expconsta,expconstb,fit_method)
% the function Kp_search performs an interative search to find the best Kp-vs-tie line function to be used to fit
% or simulate experimental FRET data surfaces
initial_time=clock;

Kpd=[min_kpd:kpd_inc:max_kpd];
Kpa=[min_kpa:kpa_inc:max_kpa];
x=[0:ntielines+1];
if (tie_stop <= tie_start)
    error('stop point of Kp fxns vs tie line must be > start point')
else
    if (~isequal(rem((tie_stop-tie_start),start_inc),0))
        error('difference between stop and start pt must be a multiple of start increment')
    end
    if (top_inc > ntielines-tie_stop) % ~isequal(rem((tie_stop-tie_start),top_inc),0)
        error('top increment must be < number tie lines - stop tie line #') % difference between stop and start pt must be a multiple of top increment
    end
end
top_start=[tie_start:start_inc:tie_stop-start_inc]; % or top_start=[tie_start:start_inc:tie_stop-start_inc]
top_stop=cell(length(top_start));
for i=1:length(top_start)
    top_stop{i}=[top_start(i)+top_inc:top_inc:tie_stop]; % or top_stop{i}=[top_start(i)+top_inc:top_inc:tie_stop]
end
nfxns_per_kp=sum(1:length(top_start)); 
nkpd_fxns=nfxns_per_kp*length(Kpd);
nkpa_fxns=nfxns_per_kp*length(Kpa);
% generate all sets of tielines
disp('generating all sets of tielines to be searched');
[tielines,alpha_beta_pts]=gen_tielines(boundary_pts,ntiepts,ntielines);
[nisets,njsets]=size(tielines);

kpd_fxns=zeros(nfxns_per_kp,ntielines+2);
kpa_fxns=zeros(nfxns_per_kp,ntielines+2);
if (nkpd_fxns >= nbest2save)
    nsave=nbest2save;
else
    nsave=nkpd_fxns;
end
best_fit_values=zeros(1,nsave);
best_tielines=cell(1,nsave);
best_abpts=cell(1,nsave);
best_surfaces=cell(1,nsave);
best_Kpd_fxns=zeros(ntielines,nsave);
best_Kpa_fxns=zeros(ntielines,nsave);
nfilled=1;
% loop over donor and acceptor Kps
for d=Kpd 
    disp('creating Kpd trapezoidal fxns');
    n=1;
    for i=top_start
        for j=top_stop{find(i == top_start)} % top_stop(#:end)
            kpd_fxns(n,:)=trapf(x,d,[i j]);
            % kpd_fxns(n,1:i)=linspace(1,d,i+1);
            % kpd_fxns(n,i+1:j)=linspace(d,d,j-i);
            % kpd_fxns(n,j+1:end)=linspace(d,1,ntielines+2-j);
            n=n+1;
        end
    end
    for a=Kpa
        disp('creating Kpa trapezoidal fxns');
        n=1;
        for i=top_start
            for j=top_stop{find(i == top_start)} % top_stop(#:end)
                kpa_fxns(n,:)=trapf(x,a,[i j]);
                n=n+1;
            end
        end
        for ndf=1:nfxns_per_kp
            kpdonr=kpd_fxns(ndf,2:ntielines+1)';
            for naf=1:nfxns_per_kp
                kpaceptr=kpa_fxns(naf,2:ntielines+1)';
                tic
                disp('using donor and acceptor Kp fxns to find best fit tielines and surface');
                % fit each calculated surface with experimental surface by least sum of squares
                % calculate FRET surface with alpha and beta fraction of boundary given by alpha_beta_pts, 
                % the sets of tielines given in tielines, and arbitrary given Kps
                all_calc_surfs=cell(nisets,njsets);
                for i=1:nisets
                    for j=1:njsets
                        tylines=tielines{i,j};
                        [FRET_surface]=gen_surface_dist(exp_surf,tylines,xdonr,xaceptr,kpdonr,kpaceptr,expconsta,expconstb);
                        all_calc_surfs{i,j}=FRET_surface;
                        
                        if (all(all(isnan(FRET_surface))))
                            disp('fret surface with following indices contain NaNs')
                            i,j
                        end
                        calc_surface=FRET_surface;
                        indices=find(~isnan(calc_surface));
                        calc_surface=calc_surface(~isnan(calc_surface));
                        if (isempty(indices))
                            i,j
                            error('all matrix elements of calculated fret surface are NaN')
                        end
                        exp_surface=exp_surf(indices);
                        % fitting calculated surface to experimental surface:  *improved fitting fxn may be needed here*
                        if (strcmp(fit_method,'sse'))
                            Error_Vector=calc_surface-exp_surface;
                            sse(i,j)=sum((Error_Vector.^2)./exp_surface); % length(Error_Vector)
                        end
                        if (strcmp(fit_method,'modal'))
                            [si,sj]=find(~isnan(FRET_surface));
                            e_surface=[si sj exp_surface];
                            c_surface=[si sj calc_surface];
                            e_surface(1:4:length(e_surface),:)=[];
                            c_surface(1:4:length(c_surface),:)=[];
                            sse(i,j)=fit_ec_surfaces(e_surface,c_surface);
                        end
                    end
                end

                % find lowest metric of all calculated surfaces
                degeneracy=nnz(min(min(sse)) == sse);
                [ilowest,jlowest]=find(min(min(sse)) == sse);
                lowest=[ilowest jlowest];
                if (degeneracy ~= 1)
                    best_fit_tielines=cell(degeneracy,1);
                    best_fit_abpts=cell(degeneracy,1);
                    best_fit_surface=cell(degeneracy,1);
                    best_fit_param=sse(ilowest(1),jlowest(1));
                    for k=1:degeneracy
                        best_fit_tielines{k}=tielines{ilowest(k),jlowest(k)};
                        best_fit_abpts{k}=alpha_beta_pts{ilowest(k),jlowest(k)};
                        best_fit_surface{k}=all_calc_surfs{ilowest(k),jlowest(k)};
                    end
                else
                    best_fit_param=sse(ilowest,jlowest);
                    best_fit_tielines=tielines{ilowest,jlowest};
                    best_fit_abpts=alpha_beta_pts{ilowest,jlowest};
                    best_fit_surface=all_calc_surfs{ilowest,jlowest};
                end 
                % update stack of best fit values, tielines, alpha beta points, and surfaces
                if (isempty(find(0 == best_fit_values)))
                    [best_fit_values,index]=sort(best_fit_values);
                    best_tielines={best_tielines{index}};
                    best_abpts={best_abpts{index}};
                    best_surfaces={best_surfaces{index}};
                    best_Kpd_fxns=best_Kpd_fxns(:,index);
                    best_Kpa_fxns=best_Kpa_fxns(:,index);
                    i=find(best_fit_param <= best_fit_values);
                    if (~isempty(i))
                        i=max(i);
                        best_fit_values(i)=best_fit_param;
                        best_tielines{i}=best_fit_tielines;
                        best_abpts{i}=best_fit_abpts;
                        best_surfaces{i}=best_fit_surface;
                        best_Kpd_fxns(:,i)=kpdonr;
                        best_Kpa_fxns(:,i)=kpaceptr;
                    end
                else
                    best_fit_values(nfilled)=best_fit_param;
                    best_tielines{nfilled}=best_fit_tielines;
                    best_abpts{nfilled}=best_fit_abpts;
                    best_surfaces{nfilled}=best_fit_surface;
                    best_Kpd_fxns(:,nfilled)=kpdonr;
                    best_Kpa_fxns(:,nfilled)=kpaceptr;
                    nfilled=nfilled+1;
                end
                clear best_fit_param;
                clear best_fit_tielines;
                clear best_fit_abpts;
                clear best_fit_surface;
                % t=toc
                disp(sprintf('execution time for surface generation and fitting = %f seconds\n',toc));
            end % naf
        end % ndf
    end % kpa
end % kpd

best_Kp_fxns=cat(3,best_Kpd_fxns,best_Kpa_fxns);
[best_fit_values,index]=sort(best_fit_values);
best_tielines={best_tielines{index}};
best_abpts={best_abpts{index}};
best_surfaces={best_surfaces{index}};
best_Kpd_fxns=best_Kpd_fxns(:,index);
best_Kpa_fxns=best_Kpd_fxns(:,index);

data_surface=exp_surf;
data_surface(isnan(best_surfaces{1}))=sqrt(NaN);

output(1).name='best Kps';
output(1).data=best_Kp_fxns;
output(2).name='best tielines';
output(2).data=best_tielines;
output(3).name='best alpha beta pts';
output(3).data=best_abpts;
output(4).name='best surfaces';
output(4).data=best_surfaces;
output(5).name='best fit values';
output(5).data=best_fit_values';

disp(sprintf('execution time = %f seconds\n',etime(clock,initial_time)));