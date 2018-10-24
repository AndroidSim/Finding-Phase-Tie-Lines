function [data]=convert_fret(raw_data,x_inc,xmax,y_inc,ymax)
%fret_analysis is a function that analyzes FRET data to pick out boundaries and %tie lines

% [X,Y]=meshgrid(0:0.02:1,0:0.02:0.66);
% Pauls_data=griddata(Pauls_data_to_fit(:,1),Pauls_data_to_fit(:,2),Pauls_data_to_fit(:,3),X,Y)';

% boundary_pts=[0.3 0.125;0.6 0.15;0.925 0.25;0.8 0.45;0.5 0.45;0.3 .45;0.2 0.35;0.1 0.225;0.15 0.175;0.2 0.15;0.8 0.2;0.6 0.45;0.4 0.45;0.4 0.135;0.88 0.35;0.84 0.4;0.24 0.4;0.24 0.137;0.7 0.175;0.5 0.13;0.15 0.3];
% bs_pts(:,1)=boundary_pts(:,1)*50+1;
% bs_pts(:,2)=boundary_pts(:,2)*50+1;

[ndata,ncol]=size(raw_data);
n=(ymax/y_inc)+1;
m=(xmax/x_inc)+1;
data=reshape(raw_data(:,ncol),m,n);
return