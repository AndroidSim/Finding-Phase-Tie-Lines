function [data]=convert_fret(raw_data,x_inc,xmax,y_inc,ymax)
%fret_analysis is a function that analyzes FRET data to pick out boundaries and %tie lines

[ndata,ncol]=size(raw_data);
n=(ymax/y_inc)+1;
m=(xmax/x_inc)+1;
data=reshape(raw_data(:,ncol),m,n);
return