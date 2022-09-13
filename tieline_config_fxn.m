function v = tieline_config_fxn(a,u)
% gamma = tieline_config_fxn(x,u)
% input: x = parameters to vary
%        u = position parameter on Lo boundary (directrix)
% output: v = (gamma) position parameter on La boundary (directrix)

% given as a fourier series
% w = pi/2;
% v = a(1)*sin(w*u) + a(2)*cos(w*u) + a(3)*sin(2*w*u) + a(4)*cos(2*w*u);
% v = (sin(a)^2)*sin(w*u) + (cos(a)^2)*sin(5*w*u);
% v = (sin(a)^2)*sin(w*u) + (cos(a)^2)*(1-cos(w*u));
% v = (sin(a(1))^2)*sin(w*u) + (cos(a(1))^2)*(1-cos(w*u)) + (sin(a(2))^2)*sin(2*w*u) + (cos(a(2))^2)*(1-cos(2*w*u));
% v = (sin(u).^2).*sin(w*u) + (cos(u).^2).*(1-cos(w*u));

% given as sum of exponentials
% v = (exp(a(1)*u) - exp(a(2)*u))./(exp(a(1)) - exp(a(2)));

% given in a form similar to fpa and fpb
v = u./(u + a.*(1-u));
% v = (sin(a(3))^2)*(u./(u + a(1).*(1-u))) + (cos(a(3))^2).*(u./(u + a(2).*(1-u)));
% v = (sin(u).^2).*(u./(u + a(1).*(1-u))) + (cos(u).^2).*(u./(u + a(2).*(1-u)));
% v = (u./(u + a(3).*(1-u))).*(u./(u + a(1).*(1-u))) + (1 - u./(u + a(3).*(1-u))).*(u./(u + a(2).*(1-u)));
% v = a(3).*(u./(u + a(1).*(1-u))) + (1 - a(3)).*(u./(u + a(2).*(1-u)));
% v = u.*(u./(u + a(1).*(1-u))) + (1 - u).*(u./(u + a(2).*(1-u)));
% v = 0.5.*((u./(u + a(1).*(1-u))) + (u./(u + a(2).*(1-u))));
% v = (1/3).*((u./(u + a(1).*(1-u))) + (u./(u + a(2).*(1-u))) + (u./(u + a(3).*(1-u))));

return