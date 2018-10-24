function [Kpd,Kpa] = Kp_fxns(x,fxn_type,gamma3,gamma4)
% Kp_fxns calculates the donor and acceptor Kps (Kpd and Kpa respectively) for tie line x with fxn given by fxn_type

switch fxn_type
    case 'para'
        % x = 0 -> 1
        pd=polyfit([0 0.5 1],[1 gamma3 1],2); 
        pa=polyfit([0 0.5 1],[1 gamma4 1],2); 
        Kpd=polyval(pd,x);
        Kpa=polyval(pa,x);
    case 'trap'
%        Kpd_fxn=trapf(x,gamma3,[gamma4 gamma5]);
%        Kpa_fxn=trapf(x,gamma6,[gamma7 gamma8]);
    otherwise
        [Kpd,Kpa]=Kp_fxns(x,'para');
end