clear all;
close all;

B=10;
m=0;

for j1=0:2*B-1
    alpha_j1(j1+1)=j1*pi/B;
    y(j1+1)=exp(-j*m*alpha_j1(j1+1));
end

figure;
plot(0:2*B-1,alpha_j1/pi,'b.');

figure;
plot(y,'b.');
axis equal;
sum(y)