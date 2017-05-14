function [ betar,betak,width,betatheta ] = rpriorfinder( A,plots )
%uses fft to find an estimate for the size of banding in image A and
%converts this to the beta value of an exponential distribution.
%   Detailed explanation goes here
n = length(A);
A = double(A);
%K = fftshift(abs( fft2(A) ).^2) ;
K1 = mean(abs( fft(A) ).^2) ;
K2 = mean(abs( fft(A') ).^2) ;
K1 = n*K1/sum(K1(:));
K2=n*K2/sum(K2(:));
% k = linspace(-n,n,n);
% [yy,xx] = meshgrid(k,k);
% r = sqrt(xx.^2+yy.^2);
% l = 2*n/n;
% D = 1:n;
% powspec=zeros(n,1);
% for d = D
%     ker = double((d<=r).*(r<=(d+l)));
%     powspec(d)=sum(sum(ker.*K));
% end

[C1,I1]=max(K1(2:floor(end/2)));
[C2,I2] = max(K2(2:floor(end/2)));
if plots==1
    close all;
    figure(1)
    plot(K1,'b--');
    hold on;
    plot(K2,'r--');
    plot(I1+1,C1,'b*');
    plot(I2+1,C2,'r*');
end
I = max(I1,I2);
betar=I+1;
betak = 1/max(C1,C2);

%calculate theta using the idea that if banding is in y direction theta 0
%and if in x direction theta is likely to be pi/2. If it seems equal it's
%likely to be pi/4
betatheta = pi/4*tanh(C1-C2)+pi/4;
%%now calculate widths.
K1 = mean(A,1);
K2 = mean(A,2);
K1 = K1(:);
K2 = K2(:);
modelfun = @(b,x)(b(1)+b(4)*sin(2*pi*(b(2)*x+b(3))/75));
beta0 = 1+rand(4,1);
t = (1:n)';
coef = nlinfit(t,K1,modelfun,beta0);
coef2 = nlinfit(t,K2,modelfun,beta0);
if plots==1
    figure(2)
    plot(K1,'b-');
    hold on;
    plot(t,coef(1)+coef(4)*sin(2*pi*(coef(2)*t+coef(3))/75),'r--');
    plot(K2,'k-');
    plot(t,coef2(1)+coef2(4)*sin(2*pi*(coef2(2)*t+coef2(3))/75),'r--');
end
 width1 = 0.5*75/coef(2);
 width2 = 0.5*75/coef2(2);
 width = sqrt(sqrt(max(width1,width2)/2)); %the division by two here is due to the fact standard deviation is typically half the width
% plot(ix1,mx1,'r*');
% plot(in1,mn1,'k*');
end

