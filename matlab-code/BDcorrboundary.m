function [ BiB,BiD,DeB,DeD,N00,N01,N11 ] = BDcorrboundary( A,m )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% m parameter defines sizes of box over which to take correlations. Used
% for real data where need to seriously worry about boundary conditions
n=length(A);
A = double(A);

k = linspace(-n,n,n);
[yy,xx] = meshgrid(k,k);
r = sqrt(xx.^2+yy.^2);
l = 2*n/n;
mid = floor(n/2);
m2 = floor(m/2);
rho1 = sum(A(:));
rho0 = sum(1-A(:));
%declare correlations memory
BiB = zeros(n,n,m); %correlations in births locations to births
DeB = zeros(n,n,m); %correlations in death loacations to births
BiD = zeros(n,n,m); %correlations in births locations to deaths
DeD = zeros(n,n,m); %correlations in deaths locations to dearth
D = 1:m;
%normalisation constants for w_xx(d)

   N00=zeros(m,1);
   N01=zeros(m,1);
   N11=zeros(m,1);


for d = D
    K = double((d<=r).*(r<=(d+l)));
     N = sum(sum(imfilter(ones(n,n),K,0,'conv'))); %sum(K(:))
    convd = imfilter(A,K,'circular','conv');
    bib = convd.*A;
    deb = convd.*(1-A);
    N11(d) = sum(bib(:));%N*(rho1+rho1);
    N01(d) = sum(deb(:)); %N*(rho0+rho1);%
    BiB(:,:,d) = bib;
    DeB(:,:,d) = deb;
    
    
    convd = imfilter(1-A,K,0,'conv');
    bid = convd.*A;
    ded = convd.*(1-A);
    N00(d) = sum(ded(:));%N*(rho0+rho0);
    BiD(:,:,d) = bid;
    DeD(:,:,d) = ded;
    
    
end
% N00 = permute(repmat(N00,[1 n n]),[2 3 1]);
% N01 = permute(repmat(N01,[1 n n]),[2 3 1]);
% N11 = permute(repmat(N11,[1 n n]),[2 3 1]);
% nn = N01;
mid = floor(n/2);
m2 = floor(m/2);
BiB=BiB(mid-m2:mid+m2,mid-m2:mid+m2,:);
DeB = DeB(mid-m2:mid+m2,mid-m2:mid+m2,:);
DeD = DeD(mid-m2:mid+m2,mid-m2:mid+m2,:);
BiD = BiD(mid-m2:mid+m2,mid-m2:mid+m2,:);
% N01 = N01(mid-m2:mid+m2,mid-m2:mid+m2);
% N00 = N00(mid-m2:mid+m2,mid-m2:mid+m2);
% N11 = N11(mid-m2:mid+m2,mid-m2:mid+m2);
% plot(squeeze(sum(sum(BiB./nn,1),2)),'rs-');
% hold on;
% plot(0,(0.2)^2,'rp');
% plot(squeeze(sum(sum(DeB./nn,1),2)),'kp--');
% plot(squeeze(sum(sum(DeB./nn,1),2)),'m^-.');
% plot(0,(0.2)*(1-0.2),'kp');
% plot(squeeze(sum(sum(DeD./nn,1),2)),'b.-');
% plot(0,(1-0.2).^2,'bp');
% legend('B to B','D to B','D to B','D to D');


end

