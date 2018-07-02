function d=fdcl_wigner_d(l,b)

d=wignerD(l,0,b,0);

for m=-l:l
    for n=-l:l
        d(m+l+1,n+l+1)=(-1)^(m+n)*d(m+l+1,n+l+1);
    end
end



end