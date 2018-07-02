function test_Phi_sum
clear all;
close all;

B=10;
l=5;

for m=-l:l
    for n=-l:l
        Phi1_sum=0;
        Phi2_sum=0;
        for j1=0:2*B-1
            alpha_j1=j1*pi/B;
            for j2=0:2*B-1
                gamma_j2=j2*pi/B;
                Phi1_sum=Phi1_sum+Phi1(m,n,alpha_j1,gamma_j2);
                Phi2_sum=Phi2_sum+Phi2(m,n,alpha_j1,gamma_j2);
            end
        end
        if abs(Phi1_sum) > 1e-5
            disp('Phi1');
            disp([m, n, Phi1_sum])
        end
        if abs(Phi2_sum) > 1e-5
            disp('Phi2');
            disp([Phi2_sum])
        end
        
    end
end

    
end

function y=Phi1(m,n,alpha,gamma)

if m*n > 0
    y=(-1)^(m-n)*cos(m*alpha+n*gamma);
elseif m*n < 0
    y=-(-1)^(m-n)*sin(m*alpha-n*gamma);
elseif ( m>0 & n==0) | ( m==0 & n >0)
    y=(-1)^(m-n)*sqrt(2)*cos(m*alpha+n*gamma);
elseif ( m<0 & n==0) | ( m==0 & n <0)
    y=-(-1)^(m-n)*sqrt(2)*sin(m*alpha-n*gamma);
elseif m==0 & n==0;
    y=1;
end

end

function y=Phi2(m,n,alpha,gamma)

if m*n > 0
    y=(-1)^(m)*sign(m)*cos(m*alpha-n*gamma);
elseif m*n < 0
    y=(-1)^(m)*sign(m)*sin(m*alpha+n*gamma);
else
    y=0;
end

end

     