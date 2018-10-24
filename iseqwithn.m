function yesorno = iseqwithn(lhs,rhs,acc)
% iseqwithn determines whether left-hand side (lhs) and the right-hand side (rhs) of the equality relational
% operator is equal within a certain given accuracy (acc).  the reason i 
% have this function is because the matlab operator and functions check for absolute equality to machine precision
%,which is way more than i need, and also those matlab functions are giving me a headache, those sobs
%[ll,wl]=size(lhs);
%[lr,wr]=size(rhs);
%if (ll ~= 1 & wl ~= 1 & lr ~= 1 & wr ~= 1)
%    error('lhs and rhs of equality must be scalars')
%end
%acc=10^-2;
if (lhs == 0 & rhs == 0) % damn 0
    yesorno=logical(1);
elseif (lhs <= 0 & rhs >= 0) 
    if (rhs-acc > 0 & 0 < lhs+acc)
        yesorno=logical(1);
    else
        yesorno=logical(0);
    end
elseif (lhs > 0 & rhs < 0)
    if (lhs-acc > 0 & 0 < rhs+acc)
        yesorno=logical(1);
    else
        yesorno=logical(0);
    end
else
    if (rhs-acc < lhs & lhs < rhs+acc)
        yesorno=logical(1);
    else
        yesorno=logical(0);
    end
end
return