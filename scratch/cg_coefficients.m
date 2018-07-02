function [rows , cols , cg] = cg_coefficients (l, m, l1 , l2)

% Calculate all CG - coefficients with m=m1+m2;
% adapted from Straub (2014).
mm = (m-l1 -l2+abs (l1 -l2+m ))/2;
n = (m+l1+l2 - abs (l1 -l2 -m))/2 - mm +1;
B = zeros (2*n ,1); C = zeros (n+1 ,1);

count = 0;
C(n) = 1;
for x = n -1: -1:1
    B(2*x) = l1 *( l1 +1)+ l2 *( l2 +1)+2*( mm+x)*(m-mm -x)-l*(l +1);
    B(2*x -1) = sqrt (( l1 *( l1 +1) -( mm+x)*( mm+x -1))*( l2 *( l2 +1) ...
        -(m-mm -x)*(m-mm -x +1)));
    C(x) = -(B(2*x)*C(x +1)+ B(2*x +1)* C(x +2))/ B(2*x -1);
    count = count + C(x ).^2;
end

C(n) = sqrt (1/( count +1));
for x = n -1: -1:1
    C(x) = -(B(2*x)*C(x +1)+ B(2*x +1)* C(x +2))/ B(2*x -1);
    count = count + C(x ).^2;
end

% Store the calculated CG - coefficients in a vector .
% The entries are ordered by ascending values of m1.
cg = C(1: n);

% Find all coefficients with m=m1+m2 ...
rows = zeros (n ,1); cols = zeros (n ,1);
count = 1;

for m1=-l1:l1
    if abs (m-m1) <= l2
        % ... and calculate their row and column indices .
        [ rows(count), cols(count)] = cg_indices (l,m,l1 ,m1 ,l2 ,m-m1 );
        count = count +1;
    end
end

end