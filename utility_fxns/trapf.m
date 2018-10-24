function y = trapf(x,kp,params)

if nargin ~= 3
    error('three arguments are required by the trap fxn');
elseif length(params) < 2
    error('the trap fxn needs at least four parameters.');
end

start = params(1); 
stop = params(2); 

if start > stop,
    error('Illegal parameter condition: start > stop');
end

y1 = zeros(size(x));
y2 = zeros(size(x));

% Compute y1
i = find(x > start);
y1(i) = ones(size(i))*kp;
i = find(x >= 0 & x <= start);
y1(i) = linspace(1,kp,length(i));

% Compute y2
i = find(x < stop);
y2(i) = ones(size(i))*kp;
i = find(x >= stop & x <= x(end));
y2(i) = linspace(kp,1,length(i));

% Compute y
if (kp >= 1)
    y = min(y1, y2);
else
    y = max(y1,y2);
end