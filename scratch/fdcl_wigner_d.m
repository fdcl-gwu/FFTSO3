function D=fdcl_wigner_D(a,b,g)


D=[(1+cos(b))/2*exp(i*(a+g)) 1/sqrt(2)*sin(b)*exp(i*a) (1-cos(b))/2*exp(i*(a-g));
    -1/sqrt(2)*sin(b)*exp(i*g) cos(b) 1/sqrt(2)*sin(b)*exp(-i*g);
    (1-cos(b))/2*exp(-i*(a-g)) -1/sqrt(2)*sin(b)*exp(-i*a) (1+cos(b))/2*exp(-i*(a+g))];

end