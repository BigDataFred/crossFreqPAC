function [dle] = dipoleLocaliziationError(c1,c2)

[dle] = sqrt((c1(1)-c2(1))^2+(c1(2)-c2(2))^2+(c1(3)-c2(3))^2);

