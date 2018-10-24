function FRET_surface = create_surface(exp_surface,bdy_pts,tielines,xdonr,xaceptr,kpdonr,kpaceptr,expconsta,expconstb)
% the function gen_surface generates a calculated FRET surface for the given set of tie lines
% in the 3D matrix tielines

[nx,ny]=size(exp_surface);
[ntie,ntpts,txy]=size(tielines);
calc_FRET=zeros(ntie,ntpts);

if (~isequal(length(kpdonr),length(kpaceptr)))
    error('lengths of donor Kp and acceptor Kp arrays not equal')
end
if (~isequal(length(kpdonr),1)) 
    if (~isequal(length(kpdonr),ntie))
        error('length of donor Kp array not equal to number of tie lines')
    end
end
if (~isequal(length(kpaceptr),1)) 
    if (~isequal(length(kpaceptr),ntie))
        error('length of acceptor Kp array not equal to number of tie lines')
    end
end

for j=1:ntie
    % calculate fraction of phases alpha and beta
    tie_line=tielines(j,:,:);
    tie_line=squeeze(tie_line);
    [alpha,beta]=alphas_betas(tie_line);
    % Fa = [Fap1 Fap2], Fax = [Fap1x Fap2x], Fay = [Fap1y Fap2y]
    Fax=[tie_line(1,1) tie_line(end,1)];
    Fay=[tie_line(1,2) tie_line(end,2)];
    Fa=interp2(exp_surface',Fax,Fay);
    
    % determine Kps by fitting the the experimental data on a tie line, interpolated from the experimental surface,
    % to the FRET equation where the Kps are adjustable parameters
    %FRET_data=interp2(exp_surface',tie_line(:,1),tie_line(:,2));
    %[x,resnorm,residual,exitflag]=lsqcurvefit(@FRET_eqn_kp_fit,[kpdonr kpaceptr],alpha',FRET_data,[1.01 1.01],[kpdonr kpaceptr],options,xdonr,xaceptr,expconsta,expconstb,Fa(1),Fa(2));
    %if (exitflag < 0)
    %    error('kp fitting did not converge to a solution')
    %end
    %kpd=x(1);
    %kpa=x(2);
    
    kpd=kpdonr(j);
    kpa=kpaceptr(j);
    
    calc_FRET(j,:)=feval(@FRET_eqn,alpha,Fa(1),Fa(2),xdonr,xaceptr,kpd,kpa,expconsta,expconstb);
end
xtie=tielines(:,:,1)';
nxcoord=ntie*ntpts;
xcoord=reshape(xtie,nxcoord,1);
ytie=tielines(:,:,2)';
nycoord=ntie*ntpts;
ycoord=reshape(ytie,nycoord,1);
nzcoord=ntie*ntpts;
zfret=reshape(calc_FRET',nzcoord,1);
xcoord=cat(1,xcoord,bdy_pts(:,1));
ycoord=cat(1,ycoord,bdy_pts(:,2));
zcoord=cat(1,zfret,interp2(exp_surface',bdy_pts(:,1),bdy_pts(:,2)));
% find points outside of boundary points
minx=min(bdy_pts(:,1));
maxx=max(bdy_pts(:,1));
miny=min(bdy_pts(:,2));
maxy=max(bdy_pts(:,2));
x=[round(minx):round(maxx)]';
y=[round(miny):round(maxy)]';
y_within=cell(length(x),1);
y_without=cell(length(x),1);
b_pts=bdy_pts;
%b_pts(end,:)=[];
for i=1:length(x)
    xx=linspace(x(i),x(i),length(y));
    [xi,yi]=polyxpoly(xx,y,b_pts(:,1),b_pts(:,2),'unique');
    if (length(yi) == 2)
        y_within{i}=y(yi(1) < y & y < yi(2));
        y_without{i}=y(yi(1) > y | y > yi(2));
    elseif (length(yi) == 1)
        y_within{i}=yi;
        y_without{i}=y(yi ~= y);
    else
        continue
        %error('either 0 or > 2 intersection points')
    end
    xv=linspace(x(i),x(i),length(y_without{i}));
    xcoord=[xcoord;xv'];
    ycoord=[ycoord;y_without{i}];
    zcoord=[zcoord;interp2(exp_surface',xv',y_without{i})];
end
coords=[xcoord ycoord zcoord];
coords=unique(coords,'rows');
[X,Y]=meshgrid(1:nx,1:ny);
FRET_surface=griddata(coords(:,1),coords(:,2),coords(:,3),X,Y)';