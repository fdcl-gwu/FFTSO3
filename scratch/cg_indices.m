function [row , col] = cg_indices (l, m, l1 , m1 , l2 , m2)

row = (l1+m1 )*(2* l2 +1)+ l2+m2 +1;
col = l^2 -(l2 -l1 )^2+ l+m+1;

end