
%fret_analysis is a function that analyzes FRET data to pick out boundaries and 
%tie lines

contour(data,30);
[Fx,Fy]=gradient(data);
hold on;
quiver(Fx,Fy);
hold off;
return