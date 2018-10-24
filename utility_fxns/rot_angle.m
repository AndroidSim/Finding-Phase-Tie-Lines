function [theta]=rot_angle(vec1,vec2,acos_or_asin)
% the function rot_angle returns the angle from which vec1 is rotated to get vec2,
% given the rotation matrix for counterclockwise rotation R=[cos(theta) -sin(theta);sin(theta) cos(theta)]
% and that vec1 and vec2 are unit vectors (normalized)
topacos=(vec1(1).*vec2(1))+(vec1(2).*vec2(2));
bottom=((vec1(1).^2)+(vec1(2).^2));
r_angle_acos=acos(topacos./bottom);
topasin=((vec2(2).*vec1(1))-(vec1(2).*vec2(1)));
r_angle_asin=asin(topasin./bottom);
%keyboard
%if (abs(r_angle_acos) ~= abs(r_angle_asin))
%    error('rotation angle from asin and acos not equal')
%else
%    theta=r_angle_asin;
%end
if (strcmp(acos_or_asin,'acos'))
    theta=r_angle_acos;
elseif (strcmp(acos_or_asin,'asin'))
    theta=r_angle_asin;
else
    error('angle returned from rot_angle must be acos or asin')
end