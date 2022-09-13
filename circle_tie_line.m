function v = circle_tie_line(x,m,x1,y1)

% v = sin(x)-y1-m*cos(x)+m*x1;

v = [x(2)-y1-m*x(1)+m*x1;
    (x(1)^2)+(x(2)^2)-1]; 