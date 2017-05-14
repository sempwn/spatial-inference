N = 100;
M=100;
 A = rand(N,N)<0.1;
% A(floor(N/2),floor(N/2))=1;
A = double(A);
kx = ((-(N-1)/2):((N-1)/2));
ky = ((-(M-1)/2):((M-1)/2));
[xx,yy] = meshgrid(ky,kx);
eps = 0.01;
l2 = 1; %1
l1 = 0.5; %0.5
k = 0.1;%0.2;
theta =pi/4; %1.5
r = 20;
%eta = 0.1;
lambda = 1;
k1 = floor(sqrt(2*log(1/0.025))*l1^2)+1; %this selects the size of the kernel such that the mass of a unit outside of it is less than 0.025
k1= k1 + mod(k1,2); %increses by one if k1 odd in order to do next operation
ker1 = exp(-(.5/l1^2)*(xx.^2+yy.^2));
ker1 = ker1/sum(ker1(:));
ker1 = ker1((N-k1)/2:(N+k1)/2,(N-k1)/2:(N+k1)/2); 

k2 = floor(sqrt(2*log(1/0.025))*l2^2)+1; %this selects the size of the kernel such that the mass of a unit outside of it is less than 0.025
k2= k2 + mod(k2,2); %increses by one if k1 odd in order to do next operation
ker2 = exp(-(.5/l2^2)*((xx).^2+(yy).^2)); %exp(-(.5/alpha)*(xx.^2+yy.^2));%
ker2 = ker2/sum(ker2(:));
ker2 = ker2((N-k2)/2:(N+k2)/2,(N-k2)/2:(N+k2)/2); 
ox = round(r*cos(theta));
oy = round(r*sin(theta));

for t = 1:2000
P1 = eps*lambda*imfilter(A,ker1,'circular','conv');%P1 = lambda*K1'*(S_vec.*(1+omega*S_vec));
        
P2 = max(0,imfilter(A,ker2,'circular','conv')); %P2 = 1 - K2'*S_vec/K;
P2 = circshift(eps*P2.^k,[ox,oy]);
%P4 = imfilter(S,ker2,'circular')-a;
%P3 = double(eps*eta*(A==1));%eta*Sp_vec;

B = rand(N,N)<double(P1); %prob of birth event
        
D = rand(N,N)<P2; %prob of death event

A = A.*(1-B-D) + 1*B;

if mod(t,50) == 0
        disp(['t = ' num2str(t)]);
        imagesc(A);
        
        colormap([0.95 0.95 0.5; 0 0.6 0]);
        title(['rho at ' num2str(mean(A(:)))]);
        axis off
        drawnow
end
end

%title('')
%set(gcf,'PaperPositionMode','auto');
%file = strcat('recovery_plot','.eps');
%print('banding_model_example','-depsc','-r0');
