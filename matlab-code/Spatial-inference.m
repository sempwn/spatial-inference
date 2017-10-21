%%EDD inference
%using a metropolis-hastings style algorithm on abs(EDD), we infer
%parameters for the banding PCA model.
%%%functions
 close all
 frac = @(x) x-floor(x);

%If fitting from simulated data, uncomment the below line
% [ A ] = bdcreator( 100);

% If fitting from binary image data, uncomment the below line and define
% the file name
% A = imread(FILENAME,'png');

[N,M] = size(A);
kx = ((-(N-1)/2):((N-1)/2));
ky = ((-(M-1)/2):((M-1)/2));
[xx,yy] = meshgrid(ky,kx);
N=max(N,M);

% define MCMC parameters

burnt = 0;
inft = 1.5E6;


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

[ betar,betak,width,betatheta ] = rpriorfinder( A,0 );

rho = mean(A(:));
betar=10;
betamu = 1/rho - 1;
betamu = 1-exp(-betamu);
mu = 0.3;
betac=10;
width=1;
betak=10;

%initialise parameters and EDD
c=1;
r = exprnd(betar);
theta = mod(normrnd(betatheta,pi/16),pi);
l2 = exprnd(width);%l2 prior
l1 = exprnd(width);%l1 prior
k = exprnd(betak); %k prior

 figure(1)
imagesc(A)
waitfor(msgbox(['r = ' num2str(betar) ' theta = ' num2str(betatheta) 'width = ' num2str(width) ' k = ' num2str(k)]));
close(1)



eps=0.01;
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

            H=squeeze(sum(sum(H1,1),2).^2)./N00.^2+squeeze(sum(sum(H3,1),2).^2)./N11.^2;
            H = squeeze(H);

            EDD(1) = -sum(H(2:end))/expl;

            prior(1) = -(l1/width+l2/width+k/betak+r/betar)+ log(double(r<25));



            ind=1;
            h = waitbar(0,'Please wait...');
ra = 0;
    for t = 1:(burnt+inft)


                    switch mod(t,4)
                        case 0
                            pl1 = max(0,l1+ normrnd(0,0.1));

                            pl2=l2;
                            pk=k;
                            ptheta=theta;
                            pr=r;
                            pc = c;
                        case 1
                            pl2 = max(0,l2+ normrnd(0,0.1));

                            pl1=l1;
                            ptheta=theta;
                            pr=r;
                            pk=k;
                            pc = c;
                        case 2
                            pk = max(0,k+normrnd(0,0.1));

                            pl1=l1;
                            pl2=l2;
                            ptheta=theta;
                            pr=r;
                            pc = c;
                        case 3
                            pr = max(0,r+ normrnd(0,0.5));

                            ptheta =  mod(theta + normrnd(0,pi/16),2*pi);
                            pl1=l1;
                            pl2=l2;
                            pk=k;
                            pc = c;
                    end



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


            oxs = round(pr*cos(ptheta));
            oys = round(pr*sin(ptheta));


            P2 = circshift(P2,[oxs,oys]);

            B = 1-exp(-eps*P1.*(1-A));
            D = 1-exp(-pk*eps.*A.*P2.^c);


            B=B(mid-m2:mid+m2,mid-m2:mid+m2);
            D = D(mid-m2:mid+m2,mid-m2:mid+m2);
            totr = 1;%sum(B(:))+sum(D(:));
            BB = repmat(B,[1 1 bbox]);
            DB = repmat(D,[1 1 bbox]);
            nn= N01;
            H1 = (-BB.*DeD+DB.*BiD)/totr; %P00
            H2 = (-DB.*BiB+BB.*DeB)/totr; %P11
            H3 = (DB.*BiB-DB.*BiD-BB.*DeB+BB.*DeD)/totr; %H3 = (2*D.*BiB-D.*BiD-B.*DeB+2*B.*DeD); %P01

            H=squeeze(sum(sum(H1,1),2).^2)./N00.^2+squeeze(sum(sum(H3,1),2).^2)./N11.^2;
            H = squeeze(H);
            pEDD = -sum(H(2:end))/expl;

            %print summary statistics
             if mod(t,100)==0
                 disp(['acceptance rate =' num2str((100-ra)/100)]);

                disp(['L = ' num2str(EDD(ind)) ', pL = ' num2str(pEDD) ', Error= ' num2str(sum(H(2:end))) 't = ' num2str(t)]);
                disp(['r = ' num2str(r) ' k = ' num2str(k) ' l1 = ' num2str(l1) ' l2 = ' num2str(l2) ' theta = ' num2str(theta) ' c = ' num2str(c)]);
                ra = 0;
             end

            pprior = -(pl1/width+pl2/width+pk/betak+pr/betar)+ log(double(pr<25));
            a = pEDD+ pprior - EDD(ind)-prior(ind); %this is the inverse of the standard EDD(t)/pEDD;

            as(ind) = a;
            %metropolis-hastings update procdure

            if (a>=  0 && ~isinf(pprior))
                EDD(ind+1) = pEDD; %update new EDD
                prior(ind+1) = pprior;
                %update to proposal parameters
                 l1 = pl1;
                 l2 = pl2;

                 r = pr;
                 theta = ptheta;

                 k = pk;
                 c=pc;
            elseif (log(rand)<a && ~isinf(pprior))
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




            waitbar(ind/(burnt+inft));
            ind = ind+1;

    end

close(h);
%save results
save('results.mat','rr','rtheta','rl1','rl2','rk','rc','EDD','A');
