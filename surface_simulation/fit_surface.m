function [p2_output,data_surface,varargout]=fit_surface(exp_surf,boundary_pts,ntiepts,ntielines,xdonr,xaceptr,kpdonr,kpaceptr,tie_type,np32search,expconsta,expconstb)
% the function fit_surface fits a generated FRET surface from a set of non-intersecting tie lines 
% to the experimental FRET surface

% varargout = p3b_output,best2or3b_output
initial_time=clock;
% find sizes of passed arrays
[nx,ny]=size(exp_surf);
[nbpts,bxy]=size(boundary_pts);
if (bxy ~= 2)
	error('num cols of boundary points not equal to 2 (x and y)')
end


disp('generating sets of tielines')
% default tieline generation
if (strcmp(tie_type,''))
    tie_type='parallel';
end
% generate set of non-intersecting trial tie lines
if (strcmp(tie_type,'parallel'))
    degree_inc=2;
    [tielines,alpha_beta_pts]=parallel_tielines(boundary_pts,ntiepts,degree_inc);   
elseif (strcmp(tie_type,'tangent'))
    [tielines,alpha_beta_pts]=tangent_tielines(boundary_pts,ntiepts,ntielines);
elseif (strcmp(tie_type,'distance'))
    [tielines,alpha_beta_pts]=gen_tielines(boundary_pts,ntiepts,ntielines);
end
if (~iscell(tielines))
    error('tielines generated is not a cell array')
else
    [nisets,njsets]=size(tielines);
end
for i=1:nisets
    for j=1:njsets
        tylines=tielines{i,j};
        if (any(any(isnan(tylines))))
            disp('tielines with following indices contain NaNs')
            i,j
        end
    end
end


disp('calculating surfaces for each set of two-phase tielines')
disp('and fitting each calculated surface to the experimental surface')
% fit each calculated surface with experimental surface by least sum of squares
% calculate FRET surface with alpha and beta fraction of boundary given by alpha_beta_pts, 
% the sets of tielines given in tielines, and arbitrary given Kps
all_calc_surfs_alfa=cell(nisets,njsets);
all_calc_surfs_beta=cell(nisets,njsets);
Kps=cell(nisets,njsets);
for i=1:nisets
    for j=1:njsets
        tylines=tielines{i,j};
        warning off;
        [FRET_surf_alpha,FRET_surf_beta,dtie_kpda]=gen_surface(exp_surf,tylines,xdonr,xaceptr,kpdonr,kpaceptr,expconsta,expconstb);
        warning on;
        all_calc_surfs_alfa{i,j}=FRET_surf_alpha;
        all_calc_surfs_beta{i,j}=FRET_surf_beta;
        Kps{i,j}=dtie_kpda(:,2:3);
        a_surface=all_calc_surfs_alfa{i,j};
        if (all(all(isnan(a_surface))))
            disp('alpha surface with following indices contain NaNs')
            i,j
        end
        b_surface=all_calc_surfs_beta{i,j};
        if (all(all(isnan(b_surface))))
            disp('beta surface with following indices contain NaNs')
            i,j
        end
        calc_surface_alfa=all_calc_surfs_alfa{i,j};
        indices_alfa=find(~isnan(calc_surface_alfa));
        calc_surface_alfa=calc_surface_alfa(~isnan(calc_surface_alfa));
        calc_surface_beta=all_calc_surfs_beta{i,j};
        indices_beta=find(~isnan(calc_surface_beta));
        calc_surface_beta=calc_surface_beta(~isnan(calc_surface_beta));
        if (isempty(indices_alfa))
            i,j
            error('all matrix elements of calculated alpha surface are NaN')
        end
        if (isempty(indices_beta))
            i,j
            error('all matrix elements of calculated beta surface are NaN')
        end
        if (indices_alfa == indices_beta)
            exp_surface=exp_surf(indices_alfa);
        else
            error('indices of not NAN in calculated surfaces alpha and beta not equal');
        end   
        Error_Vector_alfa=calc_surface_alfa-exp_surface;
        Error_Vector_beta=calc_surface_beta-exp_surface;
        sse_alfa(i,j)=sum(Error_Vector_alfa.^2)./length(Error_Vector_alfa);
        sse_beta(i,j)=sum(Error_Vector_beta.^2)./length(Error_Vector_beta);
        
        if (sse_alfa(i,j) < sse_beta(i,j))
	        least(i,j)=sse_alfa(i,j);
            aorb(i,j)='a';
        elseif (sse_beta(i,j) < sse_alfa(i,j))
            least(i,j)=sse_beta(i,j);
            aorb(i,j)='b';
        else
            least(i,j)=sse_alfa(i,j);
            aorb(i,j)='s';
        end
    end
end

disp('finding best fit 2-phase tie lines and calculated surface')
% find lowest sum of squares of both alpha and beta calculated surfaces
degeneracy=nnz(min(min(least))==least);
[ilowest,jlowest]=find(min(min(least))==least);
lowest=[ilowest jlowest];
if (degeneracy ~= 1)
    best_fit_param=least(ilowest(1),jlowest(1));
    for k=1:degeneracy
        best_fit_tielines{k}=tielines{ilowest(k),jlowest(k)};
        if (aorb(ilowest(k),jlowest(k))=='a')
            best_fit_surface{k}=all_calc_surfs_alfa{ilowest(k),jlowest(k)};
        elseif (aorb(ilowest(k),jlowest(k))=='b')
            best_fit_surface{k}=all_calc_surfs_beta{ilowest(k),jlowest(k)};
        else
            error('neither alpha nor beta surface best fit!');
        end
    end
else
    best_fit_param=least(ilowest,jlowest);
    best_fit_tielines=tielines{ilowest,jlowest};
    if (aorb(ilowest,jlowest)=='a')
        best_fit_surface=all_calc_surfs_alfa{ilowest,jlowest};
    elseif (aorb(ilowest,jlowest)=='b')
        best_fit_surface=all_calc_surfs_beta{ilowest,jlowest};
    else
        error('neither alpha nor beta surface best fit!');
    end
end
    

if (np32search >= 1) % or ~isequal(np32search,0)
    int_type='int32';
    int_np32search=int32(np32search);
    if (~isa(int_np32search,int_type))
        error('np32search is not an integer')
    end
    disp('generating tie lines with possible three-phase boundary')
    [p3_tielines,p3_ab_pts]=p3_boundary_tielines(tielines,alpha_beta_pts);
    ordered_params=sort(least(:));
    n=1;
    while (n <= np32search)
        [i,j]=find(ordered_params(n)==least);
        if (length(i) > 1)
            for m=1:length(i)
                i_p3(n+(m-1))=i(m);
                j_p3(n+(m-1))=j(m);
            end
        else
            i_p3(n)=i;
            j_p3(n)=j;
        end
        n=n+1;
    end
    if (length(i_p3) ~= length(j_p3))
        error('number indices of tielines to search for three-phase boundary not equal')
    else
        nsets=length(i_p3);
    end
    p3_indices=zeros(2,nsets);
    p3_indices(1,:)=i_p3;
    p3_indices(2,:)=j_p3;

    disp('calculating surfaces for each set of three-phase boundary tielines')
    disp('and fitting each calculated surface to the experimental surface')
    % fit each calculated surface with experimental surface by least sum of squares
    % calculate FRET surface with alpha and beta fraction of boundary given by alpha_beta_pts, 
    % the sets of tielines given in tielines, and arbitrary given Kps
    all_p3_calc_surfs_a=cell(1,nsets);
    all_p3_calc_surfs_b=cell(1,nsets);
    p3_sse_alfa=cell(1,nsets);
    p3_sse_beta=cell(1,nsets);
    p3_least=cell(1,nsets);
    for n=1:nsets
        all_p3_calc_surfs_a{n}=cell(2,length(p3_tielines{i_p3(n),j_p3(n)}));
        all_p3_calc_surfs_b{n}=cell(2,length(p3_tielines{i_p3(n),j_p3(n)}));
        p3_sse_alfa{n}=zeros(2,length(p3_tielines{i_p3(n),j_p3(n)}));
        p3_sse_beta{n}=zeros(2,length(p3_tielines{i_p3(n),j_p3(n)}));
        p3_least{n}=zeros(2,length(p3_tielines{i_p3(n),j_p3(n)}));
        p3_aorb{n}=zeros(2,length(p3_tielines{i_p3(n),j_p3(n)}));
        for k=1:length(p3_tielines{i_p3(n),j_p3(n)})
            for m=1:2
                tylines=p3_tielines{i_p3(n),j_p3(n)}{m,k};
                warning off;
                [FRET_surf_alpha,FRET_surf_beta,dtie_kpda]=gen_surface(exp_surf,tylines,xdonr,xaceptr,kpdonr,kpaceptr,expconsta,expconstb);
                warning on;
                all_p3_calc_surfs_a{n}{m,k}=FRET_surf_alpha;
                all_p3_calc_surfs_b{n}{m,k}=FRET_surf_beta;
                p3_a_surface=all_p3_calc_surfs_a{n}{m,k};
                if (all(all(isnan(p3_a_surface))))
                    disp('three phase boundary alpha surface with following indices contain NaNs')
                    n,k,m
                end
                p3_b_surface=all_p3_calc_surfs_b{n}{m,k};
                if (all(all(isnan(p3_b_surface))))
                    disp('three phase boundary beta surface with following indices contain NaNs')
                    n,k,m
                end
                calc_surface_alfa=all_p3_calc_surfs_a{n}{m,k};
                indices_alfa=find(~isnan(calc_surface_alfa));
                calc_surface_alfa=calc_surface_alfa(~isnan(calc_surface_alfa));
                calc_surface_beta=all_p3_calc_surfs_b{n}{m,k};
                indices_beta=find(~isnan(calc_surface_beta));
                calc_surface_beta=calc_surface_beta(~isnan(calc_surface_beta));
                if (isempty(indices_alfa))
                    n,k,m
                    error('all matrix elements of calculated three phase boundary alpha surface are NaN')
                end
                if (isempty(indices_beta))
                    n,k,m
                    error('all matrix elements of calculated three phase boundary beta surface are NaN')
                end
                if (indices_alfa == indices_beta)
                    exp_surface=exp_surf(indices_alfa);
                else
                    error('indices of not NAN in calculated surfaces three phase boundary alpha and beta not equal');
                end    
                Error_Vector_alfa=calc_surface_alfa-exp_surface;
                Error_Vector_beta=calc_surface_beta-exp_surface;
                p3_sse_alfa{n}(m,k)=sum(Error_Vector_alfa.^2)./length(Error_Vector_alfa);
                p3_sse_beta{n}(m,k)=sum(Error_Vector_beta.^2)./length(Error_Vector_beta);
                
                if (p3_sse_alfa{n}(m,k) < p3_sse_beta{n}(m,k))
	                p3_least{n}(m,k)=p3_sse_alfa{n}(m,k);
                    p3_aorb{n}(m,k)='a';
                elseif (p3_sse_beta{n}(m,k) < p3_sse_alfa{n}(m,k))
                    p3_least{n}(m,k)=p3_sse_beta{n}(m,k);
                    p3_aorb{n}(m,k)='b';
                else
                    p3_least{n}(m,k)=p3_sse_alfa{n}(m,k);
                    p3_aorb{n}(m,k)='s';
                end
            end
        end
    end
    
    disp('finding best fit 3-phase boundary tie lines and calculated surface')
    % find lowest sum of squares of both alpha and beta calculated surfaces
    best_fit_p3_tielines=cell(1,nsets);
    best_fit_p3_surface=cell(1,nsets);
    p3_lowest=cell(1,nsets);
    best_fit_p3_param=zeros(1,nsets);
    for n=1:nsets
        p3_degeneracy(n)=nnz(min(min(p3_least{n}))==p3_least{n});
        [mlowest,klowest]=find(min(min(p3_least{n}))==p3_least{n});
        p3_lowest{n}=[mlowest klowest];
        if (p3_degeneracy(n) ~= 1)
            best_fit_p3_param(n)=p3_least{n}(mlowest(1),klowest(1));
            for r=1:p3_degeneracy
                best_fit_p3_tielines{n}{r}=p3_tielines{i_p3(n),j_p3(n)}{mlowest(r),klowest(r)};
                if (p3_aorb{n}(mlowest(r),klowest(r))=='a')
                    best_fit_p3_surface{n}{r}=all_p3_calc_surfs_a{n}{mlowest(r),klowest(r)};
                elseif (p3_aorb{n}(mlowest(r),klowest(r))=='b')
                    best_fit_p3_surface{n}{r}=all_p3_calc_surfs_b{n}{mlowest(r),klowest(r)};
                else
                    error('neither three phase boundary alpha nor beta surface best fit!');
                end
            end
        else
            best_fit_p3_param(n)=p3_least{n}(mlowest,klowest);
            best_fit_p3_tielines{n}=p3_tielines{i_p3(n),j_p3(n)}{mlowest,klowest};
            if (p3_aorb{n}(mlowest,klowest)=='a')
                best_fit_p3_surface{n}=all_p3_calc_surfs_a{n}{mlowest,klowest};
            elseif (p3_aorb{n}(mlowest,klowest)=='b')
                best_fit_p3_surface{n}=all_p3_calc_surfs_b{n}{mlowest,klowest};
            else
                error('neither three phase boundary alpha nor beta surface best fit!');
            end
        end
    end
    best_p3_param_overall=min(best_fit_p3_param);
    p3_lowest_overall=find(min(best_fit_p3_param)==best_fit_p3_param);
end

data_surface=exp_surf;
data_surface(isnan(best_fit_surface))=sqrt(NaN);

if (np32search >= 1)
    disp('comparing tie lines with 3-phase boundary and tie lines without 3-phase boundary')
    for n=1:nsets
        if (best_fit_p3_param(n) < least(i_p3(n),j_p3(n)))
            best_2or3(n)=3;
            bestfit_param_2or3(n)=best_fit_p3_param(n);
        elseif (least(i_p3(n),j_p3(n)) < best_fit_p3_param(n))
            best_2or3(n)=2;
            bestfit_param_2or3(n)=least(i_p3(n),j_p3(n));
        else
            best_2or3(n)=0;
            bestfit_param_2or3(n)=0;
        end
    end
end

%p2_output={all_calc_surfs_alfa,all_calc_surfs_beta,tielines,alpha_beta_pts,best_fit_tielines,best_fit_surface,aorb,lowest,least,degeneracy};
%p3_output={all_p3_calc_surfs_a,all_p3_calc_surfs_b,p3_tielines,p3_ab_pts,best_fit_p3_tielines,best_fit_p3_surface,p3_aorb,p3_lowest,p3_least,p3_degeneracy,best_p3_param_overall};
%best2or3_output={best_2or3,bestfit_param_2or3};

p2_output(1).name='all_calc_surfs_alfa';
p2_output(1).data=all_calc_surfs_alfa;
p2_output(2).name='all_calc_surfs_beta';
p2_output(2).data=all_calc_surfs_beta;
p2_output(3).name='tielines';
p2_output(3).data=tielines;
p2_output(4).name='alpha_beta_pts';
p2_output(4).data=alpha_beta_pts;
p2_output(5).name='best_fit_tielines';
p2_output(5).data=best_fit_tielines;
p2_output(6).name='best_fit_surface';
p2_output(6).data=best_fit_surface;
p2_output(7).name='aorb';
p2_output(7).data=aorb;
p2_output(8).name='lowest';
p2_output(8).data=lowest;
p2_output(9).name='least';
p2_output(9).data=least;
p2_output(10).name='degeneracy';
p2_output(10).data=degeneracy;
p2_output(11).name='top_4p3b_indices';
p2_output(11).data=p3_indices;
p2_output(12).name='Kps_donr_acceptr';
p2_output(12).data=Kps;

if (np32search >= 1)
    p3b_output(1).name='all_p3_calc_surfs_a';
    p3b_output(1).data=all_p3_calc_surfs_a;
    p3b_output(2).name='all_p3_calc_surfs_b';
    p3b_output(2).data=all_p3_calc_surfs_b;
    p3b_output(3).name='p3_tielines';
    p3b_output(3).data=p3_tielines;
    p3b_output(4).name='p3_ab_pts';
    p3b_output(4).data=p3_ab_pts;
    p3b_output(5).name='best_fit_p3_tielines';
    p3b_output(5).data=best_fit_p3_tielines;
    p3b_output(6).name='best_fit_p3_surface';
    p3b_output(6).data=best_fit_p3_surface;
    p3b_output(7).name='p3_aorb';
    p3b_output(7).data=p3_aorb;
    p3b_output(8).name='p3_lowest';
    p3b_output(8).data=p3_lowest;
    p3b_output(9).name='p3_least';
    p3b_output(9).data=p3_least;
    p3b_output(10).name='p3_degeneracy';
    p3b_output(10).data=p3_degeneracy;
    p3b_output(11).name='best_p3_param_overall';
    p3b_output(11).data=best_p3_param_overall;
    p3b_output(12).name='p3_lowest_overall';
    p3b_output(12).data=p3_lowest_overall;

    best2or3b_output(1).name='best_2or3';
    best2or3b_output(1).data=best_2or3;
    best2or3b_output(2).name='bestfit_param_2or3';
    best2or3b_output(2).data=bestfit_param_2or3;
    
    for i=1:nargout-2
        varargout{1}=p3b_output;
        varargout{2}=best2or3b_output;
    end
end

disp(sprintf('execution time = %f seconds\n',etime(clock,initial_time)));