function [FRET_surf_alpha,FRET_surf_beta,dtie_kpda]=gen_surface(exp_surface,tielines,xdonr,xaceptr,kpdonr,kpaceptr,expconsta,expconstb)
% the function gen_surface generates a calculated FRET surface for the given set of tie lines
% in the 3D matrix tielines

[nx,ny]=size(exp_surface);
[ntie,ntpts,txy]=size(tielines);
calc_FRET_alpha=zeros(ntie,ntpts);
calc_FRET_beta=zeros(ntie,ntpts);
kpda=zeros(ntie,2);
%options=optimset('MaxFunEvals',1000,'TolFun',1e-10,'TolX',1e-10,'MaxIter',1000);

% calculate distances of all tie lines
for i=1:ntie
    tie_dist(i)=sqrt(((tielines(i,1,1)-tielines(i,end,1)).^2)+((tielines(i,1,2)-tielines(i,end,2)).^2));
end
kpd_slope=(kpdonr-1)./(max(tie_dist)-0);
kpa_slope=(kpaceptr-1)./(max(tie_dist)-0);

for j=1:ntie
    % calculate fraction of phases alpha and beta
    tie_line=tielines(j,:,:);
    tie_line=squeeze(tie_line);
    [alpha,beta]=alphas_betas(tie_line);
    % Fa = [Fap1 Fap2], Fax = [Fap1x Fap2x], Fay = [Fap1y Fap2y]
    Fax=[tie_line(1,1) tie_line(end,1)];
    Fay=[tie_line(1,2) tie_line(end,2)];
    Fa=interp2(exp_surface',Fax,Fay);
    
    % calculate FRET along tie_line using Kps that are a linear function of tie line length
    % at critical point all kps == 1
    % [length kp]: from [0 1] to [max_length kpdonr] and from [0 1] to [max_length kpaceptr]
    kpd=kpd_slope.*(tie_dist(j)-0)+1;
    kpa=kpa_slope.*(tie_dist(j)-0)+1;
    
    % determine Kps by fitting the the experimental data on a tie line, interpolated from the experimental surface,
    % to the FRET equation where the Kps are adjustable parameters
    %FRET_data=interp2(exp_surface',tie_line(:,1),tie_line(:,2));
    %[x,resnorm,residual,exitflag]=lsqcurvefit(@FRET_eqn_kp_fit,[kpdonr kpaceptr],alpha',FRET_data,[1.01 1.01],[kpdonr kpaceptr],options,xdonr,xaceptr,expconsta,expconstb,Fa(1),Fa(2));
    %if (exitflag < 0)
    %    error('kp fitting did not converge to a solution')
    %end
    %kpd=x(1);
    %kpa=x(2);
    
    %kps_str = mat2str(dtie_kpda(:,2:3))
    %kps_str=[1.73691462212056 1.76458854169224;2.15397509342345 2.24085707133219;2.54836037815066 2.71388599797179;2.84548903993077 3.08643839004484;3.04993442732429 3.3538072844766;3.19140191845434 3.54512196054547;3.31354355675135 3.71407387897298;3.46886063380908 3.93480597882162;3.74616485844583 4.33490914880498;4.19606946606177 5.04383697105172;4.74541587702006 6.03511835326094;5.25905273803078 7.19390028549455;5.79576492275115 9.08049244017748;6.91924830242094 9.99999999999997;7.49205211749208 9.99999999999938;6.52206418774663 9.99999999999998;9.99999999999998 5.14119965817009;9.99999999999998 4.37772214200639;9.99999999999998 3.80658413440884;9.99999999999998 3.49169673932302;9.99999999999996 3.15342725328385;9.99999999999907 2.85902400090665;9.99999999999998 2.68637719017987;9.99999999999998 2.59848945025592;9.99999999999997 2.60284278331668;9.99999999999998 2.56344757962748;6.77651580099789 2.6329771181408;3.09873172812351 3.5695609334546;2.42749761123667 2.58776594196619;1.56801264059107 1.56911061257134];
    %eval(kps_str);
    %kpd=kps_str(j,1);
    %kpa=kps_str(j,2);
    
    dtie_kpda(j,:)=[tie_dist(j) kpd kpa];
    
    calc_FRET_alpha(j,:)=feval(@FRET_eqn,alpha,Fa(1),Fa(2),xdonr,xaceptr,kpd,kpa,expconsta,expconstb);
    calc_FRET_beta(j,:)=feval(@FRET_eqn,beta,Fa(1),Fa(2),xdonr,xaceptr,kpd,kpa,expconsta,expconstb);
end
xtie=tielines(:,:,1)';
nxcoord=ntie*ntpts;
xcoord=reshape(xtie,nxcoord,1);
ytie=tielines(:,:,2)';
nycoord=ntie*ntpts;
ycoord=reshape(ytie,nycoord,1);
nzcoord=ntie*ntpts;
zfret_alfa=reshape(calc_FRET_alpha',nzcoord,1);
zfret_beta=reshape(calc_FRET_beta',nzcoord,1);
[X,Y]=meshgrid(1:nx,1:ny);
FRET_surf_alpha=griddata(xcoord,ycoord,zfret_alfa,X,Y)';
FRET_surf_beta=griddata(xcoord,ycoord,zfret_beta,X,Y)';