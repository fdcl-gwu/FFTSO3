function tmp_Euler323
e1=[1 0 0]';
e2=[0 1 0]';
e3=[0 0 1]';

% expmso3(e1)-exp1(1)
% expmso3(e2)-exp2(1)
% expmso3(e3)-exp3(1)

syms a b g;
R=exp3(a)*exp2(b)*exp3(g)

[a b g]=R2Euler323(Euler3232R(-0.5,0,+1))



end


function R=exp1(c)
R=[1 0 0;
    0 cos(c) -sin(c);
    0 sin(c) cos(c)];
end
function R=exp3(c)
R=[cos(c) -sin(c) 0;
    sin(c) cos(c) 0;
    0 0 1];
end
function R=exp2(c)
R=[cos(c) 0 sin(c);
    0 1 0
    -sin(c) 0 cos(c)];
end

function R=Euler3232R(a,b,g)
R=exp3(a)*exp2(b)*exp3(g);
end

function [a,b,g]=R2Euler323(R)
b=acos(R(3,3));

if b~=0
    g=atan2(R(3,2),-R(3,1));
    a=atan2(R(2,3),R(1,3));
else
    a=acos(R(1,1))/2;
    g=a;
end

end