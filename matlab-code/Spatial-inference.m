%%EDD inference
%using a metropolis-hastings style algorithm on abs(EDD), we infer
%parameters for the banding PCA model and use it to make science fun. 
 %%%functions
 close all
 frac = @(x) x-floor(x);
 
% [ A ] = bdcreator( 100);
%% 
% A = imread('banding3.png','png');
% A = sum(A,3)<765;
% A = (imresize(A*265,0.25,'bicubic'));%treats as an image so this multiplying factor smooths out image when it is resized
% A = double(A>200); %threshold above some value to give what (for banding site 2 and 1 this threshold is 1). 
% [N,M] = size(A);
% N=min(N,M);
% if mod(N,2)==1
%     N=N-1; %if odd then make even
% end
% M=N;
% A = A(1:100,50:149);%A=A(41:140,41:140);
% A=double(BWn);
[N,M] = size(A);
kx = ((-(N-1)/2):((N-1)/2)); 
ky = ((-(M-1)/2):((M-1)/2));
[xx,yy] = meshgrid(ky,kx);
N=max(N,M);
% define MCMC parameters

burnt = 0;
inft = 1.5E6;

%shpp = 1; %likelihood shape parameter
%declare memory for parameters and EDD records
rl2 = zeros(inft,1);
rl1 = zeros(inft,1);
rk = zeros(inft,1);
rr = zeros(inft,1);
rtheta = zeros(inft,1);
rc = zeros(inft,1);
rlambda = zeros(inft,1);
as = zeros(inft,1);

EDD = zeros(inft,1);
prior = zeros(inft,1);

bbox = 80; %size of interior box over which zeta is calculated
[ BiB,BiD,DeB,DeD,N00,N01,N11 ] = BDcorrboundary( A,bbox );
% [ betar,betak,width,betatheta ] = rpriorfinder( A,0 );
[ betar,betak,width,betatheta ] = rpriorfinder( A,0 );
rho = mean(A(:));
betar=10;
betamu = 1/rho - 1;
betamu = 1-exp(-betamu);
mu = 0.3;
betac=10;
width=1;
%  betar = 1;
%  betatheta=1.5;
  betak=10;
%initialise parameters and EDD
%c=1;
%  r = exprnd(betar);
%  theta = mod(normrnd(betatheta,pi/16),pi);
%  l2 = exprnd(width);%gamrnd(width^2/width,width/width); %calculate from distance between bands??
%  l1 = exprnd(width);%gamrnd(width^2/width,width/width); %these need to be widths, calculate one from width of band
%  k = exprnd(betak); %hardest one to form a prior on. Maybe relate it to the strength of the signal in the power spectrum. i.e. if it is weak then competition probably is as well.
%  mu = exprnd(betamu);
 figure(1)
imagesc(A)
waitfor(msgbox(['r = ' num2str(betar) ' theta = ' num2str(betatheta) 'width = ' num2str(width) ' k = ' num2str(k)]));
close(1)

%[ EDD(1) ] = EDDtempcalc( A,l1,l2,k,r,theta,DeD,BiB,DeB,BiD,N01,N,xx,yy,1 );

eps=0.01;%eps = 0.01;
            k1 = floor(sqrt(2*log(1/0.025))*l1^2)+1; %this selects the size of the kernel such that the mass of a unit outside of it is less than 0.025
            k1= k1 + mod(k1,2); %increses by one if k1 odd in order to do next operation
            if k1 > N-2
                k1 = N-2;
            end
            ker1 = exp(-(.5/l1^2)*(xx.^2+yy.^2));
            ker1 = ker1/sum(ker1(:));
            ker1 = ker1((N-k1)/2:(N+k1)/2,(N-k1)/2:(N+k1)/2); 

            k2 = floor(sqrt(2*log(1/0.025))*l2^2)+1; %this selects the size of the kernel such that the mass of a unit outside of it is less than 0.025
            k2= k2 + mod(k2,2); %increses by one if k1 odd in order to do next operation
            if k2 > N-2
                k2 = N-2;
            end
            ker2 = exp(-(.5/l2^2)*((xx).^2+(yy).^2)); %exp(-(.5/alpha)*(xx.^2+yy.^2));%
            ker2 = ker2/sum(ker2(:));
            ker2 = ker2((N-k2)/2:(N+k2)/2,(N-k2)/2:(N+k2)/2); 


            P1 = imfilter(A,ker1,'circular','conv');%P1 = lambda*K1'*(S_vec.*(1+omega*S_vec));

            P2 = imfilter(A,ker2,'circular','conv'); %P2 = 1 - K2'*S_vec/K;
            
            
            oxs = round(r*cos(theta));
            oys = round(r*sin(theta));
            P2 = circshift(P2,[oxs,oys]);
            
            B = 1-exp(-eps*P1.*(1-A));
            D = 1-exp(-k*eps.*A.*P2.^c);
            
%             B = double(P1.*(A==0)); %prob of birth event
%            
%            
%             D = double(eps.*P2.*(A==1)); %prob of death event
            mid = floor(N/2);
            m2 = floor(bbox/2);
            B=B(mid-m2:mid+m2,mid-m2:mid+m2);
            D = D(mid-m2:mid+m2,mid-m2:mid+m2);
            totr = 1;%sum(B(:))+sum(D(:));
            BB = repmat(B,[1 1 bbox]);
            DB = repmat(D,[1 1 bbox]);
            nn= N01;
            H1 = (-BB.*DeD+DB.*BiD)/totr; %P00
            H2 = (-DB.*BiB+BB.*DeB)/totr; %P11
            H3 = (DB.*BiB-DB.*BiD-BB.*DeB+BB.*DeD)/totr; %H3 = (2*D.*BiB-D.*BiD-B.*DeB+2*B.*DeD); %P01
            %H = -B.*(DeD)+D.*(BiD)-D.*(BiB)+B.*(DeB)-D.*(BiD)+B.*(DeD);
            H=squeeze(sum(sum(H1,1),2).^2)./N00.^2+squeeze(sum(sum(H3,1),2).^2)./N11.^2;
            H = squeeze(H);
            %pEDD = exp(-sum(H(2:end))/temp);
            %EDD(1) = gampdf(sum(H(2:end)),1/2,(2/1)*expl); use if
            %likelihood not logged
            %EDD(1) = (shpp(1)-1)*log(sum(H(2:end))) - sum(H(2:end))/(shpp(2));
            EDD(1) = -sum(H(2:end))/expl;
            %prior(1)=log(exppdf(l1,width)*exppdf(l2,width)*exppdf(k,betak)*exppdf(r,betar)*(r<25));
            prior(1) = -(l1/width+l2/width+k/betak+r/betar)+ log(double(r<25)); 
            
            
            %ind=2;
            ind=1;
            h = waitbar(0,'Please wait...');
ra = 0;
    for t = 1:(burnt+inft)
            
                
                    switch mod(t,4)
                        case 0
                            pl1 = max(0,l1+ normrnd(0,0.1)); %0.01
                            %pl1 = l1;
                            pl2=l2;
                            pk=k;
                            ptheta=theta;
                            pr=r;
                            pc = c;
                        case 1
                            pl2 = max(0,l2+ normrnd(0,0.1)); %0.01
                            %pl2=l2;
                            pl1=l1;
                            ptheta=theta;
                            pr=r;
                            pk=k;
                            pc = c;
                        case 2
                            pk = max(0,k+normrnd(0,0.1)); %0.01 %need to change this, the boundary is not implemented correctly
                            %pk=10;
                            pl1=l1;
                            pl2=l2;
                            ptheta=theta;
                            pr=r;
                            pc = c;
                        case 3
                            pr = max(0,r+ normrnd(0,0.5)); %0.05 %.5*frac(o+normrnd(0,0.05));
                            
                            ptheta =  mod(theta + normrnd(0,pi/16),2*pi); %pi/32 %wrap around co-ordinates since director
                            pl1=l1;
                            pl2=l2;
                            pk=k;
                            pc = c;
%                         case 4
%                             pr = r; %0.05 %.5*frac(o+normrnd(0,0.05));
%                             
%                             ptheta =  theta; %pi/32 %wrap around co-ordinates since director
%                             pl1=l1;
%                             pl2=l2;
%                             pk=k;
%                             pc = max(0,k+normrnd(0,0.01));
                            
                    end
                    %             pl1=1;
%             pl2=3;
%             pk=2;
%             pr=10;
%             ptheta=1.5;
                    
            
            k1 = floor(sqrt(2*log(1/0.025))*pl1^2)+1; %this selects the size of the kernel such that the mass of a unit outside of it is less than 0.025
            k1= k1 + mod(k1,2); %increses by one if k1 odd in order to do next operation
            if k1 > N-2
                k1 = N-2;
            end
            ker1 = exp(-(.5/pl1^2)*(xx.^2+yy.^2));
            ker1 = ker1/sum(ker1(:));
            ker1 = ker1((N-k1)/2:(N+k1)/2,(N-k1)/2:(N+k1)/2); 

            k2 = floor(sqrt(2*log(1/0.025))*pl2^2)+1; %this selects the size of the kernel such that the mass of a unit outside of it is less than 0.025
            k2= k2 + mod(k2,2); %increses by one if k1 odd in order to do next operation
            if k2 > N-2
                k2 = N-2;
            end
            ker2 = exp(-(.5/pl2^2)*((xx).^2+(yy).^2)); %exp(-(.5/alpha)*(xx.^2+yy.^2));%
            ker2 = ker2/sum(ker2(:));
            ker2 = ker2((N-k2)/2:(N+k2)/2,(N-k2)/2:(N+k2)/2); 


            P1 = imfilter(A,ker1,'circular','conv');%P1 = lambda*K1'*(S_vec.*(1+omega*S_vec));

            P2 = imfilter(A,ker2,'circular','conv'); %P2 = 1 - K2'*S_vec/K;
%             Dus = P2.^pk;
            
            oxs = round(pr*cos(ptheta));
            oys = round(pr*sin(ptheta));
            %D = circshift(eps*Dus,[oxs,oys]);
            
            P2 = circshift(P2,[oxs,oys]);
%             B = double(min(1,P1.*(A==0))); %prob of birth event
%         
%             D = double(min(1,D.*(A==1))); %prob of death event
            B = 1-exp(-eps*P1.*(1-A));
            D = 1-exp(-pk*eps.*A.*P2.^c);
            %B = double(P1.*(A==0)); %prob of birth event
            
        
            %D = double(eps.*P2.*(A==1)); %prob of death event
            
            B=B(mid-m2:mid+m2,mid-m2:mid+m2);
            D = D(mid-m2:mid+m2,mid-m2:mid+m2);
            totr = 1;%sum(B(:))+sum(D(:));
            BB = repmat(B,[1 1 bbox]);
            DB = repmat(D,[1 1 bbox]);
            nn= N01;
            H1 = (-BB.*DeD+DB.*BiD)/totr; %P00
            H2 = (-DB.*BiB+BB.*DeB)/totr; %P11
            H3 = (DB.*BiB-DB.*BiD-BB.*DeB+BB.*DeD)/totr; %H3 = (2*D.*BiB-D.*BiD-B.*DeB+2*B.*DeD); %P01
            %H = -B.*(DeD)+D.*(BiD)-D.*(BiB)+B.*(DeB)-D.*(BiD)+B.*(DeD);
            H=squeeze(sum(sum(H1,1),2).^2)./N00.^2+squeeze(sum(sum(H3,1),2).^2)./N11.^2;
            H = squeeze(H);
            
            %pEDD = gampdf(sum(H(2:end)),1/2,(2/1)*expl); use if not
            %likelihood not log
            %pEDD = (shpp(1)-1)*log(sum(H(2:end))) - sum(H(2:end))/(shpp(2));
            pEDD = -sum(H(2:end))/expl;
            %calculate prior
             if mod(t,100)==0
                 disp(['acceptance rate =' num2str((100-ra)/100)]);
                 
                disp(['L = ' num2str(EDD(ind)) ', pL = ' num2str(pEDD) ', Error= ' num2str(sum(H(2:end))) 't = ' num2str(t)]);
                disp(['r = ' num2str(r) ' k = ' num2str(k) ' l1 = ' num2str(l1) ' l2 = ' num2str(l2) ' theta = ' num2str(theta) ' c = ' num2str(c)]); 
                ra = 0;
             end
            %pprior=log(exppdf(pl1,width)*exppdf(pl2,width)*exppdf(pk,betak)*exppdf(pr,betar)*(pr<25));
            pprior = -(pl1/width+pl2/width+pk/betak+pr/betar)+ log(double(pr<25)); 
            a = pEDD+ pprior - EDD(ind)-prior(ind); %this is the inverse of the standard EDD(t)/pEDD;
            %a = a*pprior/prior(ind);
            as(ind) = a;
            %metropolis-hastings update procdure, look it up on the wikipedia
            %article
            %
            
            if (a>=  0 && ~isinf(pprior)) % a is now a log
                EDD(ind+1) = pEDD; %update new EDD
                prior(ind+1) = pprior;
                %update to proposal parameters
                 l1 = pl1;
                 l2 = pl2; 
                
                 r = pr;
                 theta = ptheta;
                
                 k = pk;
                 c=pc;
            elseif (log(rand)<a && ~isinf(pprior)) %  a is now a log
                EDD(ind+1) = pEDD; %update new EDD
                prior(ind+1) = pprior;
                %update to proposal parameters
                 l1 = pl1;
                 l2 = pl2; 
                
                 r = pr;
                 theta = ptheta;
                
                 k = pk;
                 c=pc;
            else
                prior(ind+1) = prior(ind);
                EDD(ind+1) = EDD(ind);
                ra = ra+1;
                %don't update proposal parameters
            end
            if (r>25)
                pause
            end
            if t> 0 %if greater than burn time record parameters.
                rl2(ind) = l2;
                rl1(ind) = l1;
                rk(ind) = k;
                rr(ind) = r;
                rtheta(ind)=theta;
                %rc(ind)=c;
            end
            %random walk on parameter space
    
            
    
            waitbar(ind/(burnt+inft));
            ind = ind+1;
            
    end

close(h);
%MBMCMCsite.mat - mussels for large image.
save('seagrasssite2_postphd.mat','rr','rtheta','rl1','rl2','rk','rc','EDD','A');
%% plots
% close all;
% set(0,'defaultaxesfontsize',16,'defaultaxeslinewidth',1.2,...
%     'defaultlinelinewidth',1.8,'defaultpatchlinewidth',.7);
% set(0,'defaultAxesFontName', 'Arial')
% set(0,'defaultTextFontName', 'Arial')
% load('MusselsMCMCsite1.mat');
% betak = 1;
% betar=5;
% width = 1;
% 
% ts = 6E5:length(rr);
% figure(1)
% %%%scatter plots
% for i =1:5
%     for j=1:5
%         switch i
%             case 1
%                 x = rr;
%             case 2
%                 x = rtheta;
%             case 3
%                 x = rl1;
%             case 4
%                 x = rl2;
%             case 5
%                 x = rk;
%         end
%         switch j
%             case 1
%                 y = rr;
%             case 2
%                 y = rtheta;
%             case 3
%                 y = rl1;
%             case 4
%                 y = rl2;
%             case 5
%                 y = rk;
%         end
%         I = i+5*(j-1);
% subplot(5,5,I);
% X = [x(ts)';y(ts)'];
% X=X';
% smoothhist2D(X,5,[100, 100],0);
% if j ==1
%     switch i
%             case 1
%                 title('r');
%             case 2
%                 title('\theta');
%             case 3
%                 title('l_1');
%             case 4
%                 title('l_2');
%             case 5
%                 title('k');
%     end
% end
% 
% if i ==1
%     switch j
%             case 1
%                 ylabel('r','rot',0);
%             case 2
%                 ylabel('\theta','rot',0);
%             case 3
%                 ylabel('l_1','rot',0);
%             case 4
%                 ylabel('l_2','rot',0);
%             case 5
%                 ylabel('k','rot',0);
%     end
% end
% 
%     end
% end
% %%diagonal plots
% subplot(5,5,1);
% [n,xout]=hist(rr(ts),300);
% bar(xout,n/sum(n),1);
% hold on; plot(xout,exppdf(xout,betar)/sum(exppdf(xout,betar)),'r-','LineWidth',3)
% ly = ylim;
% % line([20 20],[ly(1,1) ly(1,2)],'Color','g','LineWidth',3);
% predv=mean(rr(ts));
% line([predv predv],[ly(1,1) ly(1,2)],'Color','b','LineWidth',3);
%     title('r');
%     ylabel('r');
% 
% axis([0 100 0 0.04]);
% 
% subplot(5,5,7);
% [n,xout]=hist(rtheta(ts),200);
% bar(xout,n/sum(n),1);
% hold on; plot(xout,ones(length(xout),1)/sum(xout),'r-','LineWidth',3);
% ly = ylim;
% % line([pi/2 pi/2],[ly(1,1) ly(1,2)],'Color','g','LineWidth',3);
% predvx = mean(cos(rtheta(ts)));
% predvy = mean(sin(rtheta(ts)));
% predv=abs(atan(predvy/predvx));
% line([predv predv],[ly(1,1) ly(1,2)],'Color','b','LineWidth',3);
% %axis([0 0.5 0 4000]);
% 
% subplot(5,5,13);
% [n,xout]=hist(rl1(ts),200);
% bar(xout,n/sum(n),1);
% hold on; plot(xout,exppdf(xout,width)/sum(exppdf(xout,width)),'r-','LineWidth',3);
% ly = ylim;
% % line([1 1],[ly(1,1) ly(1,2)],'Color','g','LineWidth',3);
% predv=mean(rl1(ts));
% line([predv predv],[ly(1,1) ly(1,2)],'Color','b','LineWidth',3);
% 
% subplot(5,5,19);
% [n,xout]=hist(rl2(ts),200);
% bar(xout,n/sum(n),1);
% hold on; plot(xout,exppdf(xout,width)/sum(exppdf(xout,width)),'r-','LineWidth',3);
% ly = ylim;
% % line([3 3],[ly(1,1) ly(1,2)],'Color','g','LineWidth',3);
% predv=mean(rl2(ts));
% line([predv predv],[ly(1,1) ly(1,2)],'Color','b','LineWidth',3);
% 
% 
% subplot(5,5,25);
% hist(rk(ts),30);
% [n,xout]=hist(rk(ts),200);
% bar(xout,n/sum(n),1);
% hold on; plot(xout,exppdf(xout,betak)/sum(exppdf(xout,betak)),'r-','LineWidth',3);
% ly = ylim;
% %line([0.2 0.2],[ly(1,1) ly(1,2)],'Color','g','LineWidth',3);
% predv=mean(rk(ts));
% line([predv predv],[ly(1,1) ly(1,2)],'Color','b','LineWidth',3);
% % 
% % set(gcf,'paperunits','centimeters')
% % set(gcf,'papersize',[40,40]) % Desired outer dimensions of figure
% % 
% % set(gcf,'paperposition',[0,0,40,40]) % Place plot on figure
% % file = strcat(pwd,'\\','SAMHbandingsite2','.pdf');
% % print('-dpdf',file);
% %%
% % subplot(3,2,4);
% % hist(ro,30);
% % title('o');
% % subplot(3,2,5);
% % hist(reta,30);
% % title('\eta');
% 
% figure(2); 
% %hax = axes;
% semilogy(EDD);
% hold on
% 
% %line([burnt burnt],get(hax,'YLim'),'Color',[1 0 0]);
% xlabel('MCMC time');
% ylabel('Pseudo Likelihood');
%  %axis([burnt burnt+inft 0 0.5]);
%  axis tight
%  
% %   set(gcf,'paperunits','centimeters')
% % set(gcf,'papersize',[15,15]) % Desired outer dimensions of figure
% % 
% % set(gcf,'paperposition',[0,0,15,15]) % Place plot on figure
% % file = strcat(pwd,'\\','MHex1Likelihood','.pdf');
% % print('-dpdf',file);
%  
%  figure(3)
% 
% 
%  plot(ts,rr(ts),'b--',ts,rtheta(ts),'k.-');
%  xlabel('MCMC Time (t)');
%  ylabel('Offset');
%  hold on;
%  plot(ts(end),0,'bp','MarkerSize',15);
%  plot(ts(end),20,'kp','MarkerSize',15);
%  
% %  set(gcf,'paperunits','centimeters')
% % set(gcf,'papersize',[15,15]) % Desired outer dimensions of figure
% % 
% % set(gcf,'paperposition',[0,0,15,15]) % Place plot on figure
% % file = strcat(pwd,'\\','MHex1oxy','.pdf');
% % print('-dpdf',file);
%  figure(4)
%  plot(ts,rl1(ts),'b--',ts,rl2(ts),'k.-');
%   xlabel('MCMC Time (t)');
%  ylabel('\tau_i');
%  hold on;
%  plot(ts(end),1,'bp','MarkerSize',15);
%  plot(ts(end),3,'kp','MarkerSize',15);
% %  
% %   set(gcf,'paperunits','centimeters')
% % set(gcf,'papersize',[10,10]) % Desired outer dimensions of figure
% % 
% % set(gcf,'paperposition',[0,0,10,10]) % Place plot on figure
% % file = strcat(pwd,'\\','KorcakExample','.pdf');
% % print('-dpdf',file);
% 
%   
% figure(5)
% imagesc(A)
% colormap([0.95 0.95 0.5; 0 0.6 0]);
% mytheta = mode(round(rtheta*100)/100);%-pi;
% myr = mode(round(rr));
% [xx,yy] = meshgrid(1:10:100,1:10:100);
% hold on;
% % line([50, 50+myr*cos(mytheta)],[50,50+myr*sin(mytheta)],'Color','w');
% quiver(xx,yy,myr*cos(mytheta)*sign(xx),myr*sin(mytheta)*sign(yy),'LineWidth',3,'Color','k');
