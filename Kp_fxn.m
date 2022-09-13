function k = Kp_fxn(c,u,boundary)
% k = Kp_fxn(u)
% input: u = position parameter on Lo boundary
% output: k = Kp for tieline with u parameter

% given as a taylor (power) series
switch boundary.config
    case '2critpts'
        k = 1 - (1/2).*(c.*u).*(1-u);
    case '1critpt/1endtl'
%         k = 1 + c(1).*u + (c(2)/2).*u.^2;
        k = 1 + c(1).*u + c(2).*u.^2;
    case '2endtls'
%         k = c(1) + c(2).*u + (c(3)/2).*u.^2;
        k = c(1) + c(2).*u + c(3).*u.^2;
    otherwise
        error('invalid boundary configuration');
end

return