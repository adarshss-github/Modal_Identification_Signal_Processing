function [wnat,Phi,zeta] = efdd(Y,N,fs,nsing,opt)

%===========================================================================================================================================================
%***********************************************************************************************************************************************************
%
%                                      --------------------------------------------------
%                                                          PREAMBLE
%                                      --------------------------------------------------
%
%[wnat,Phi,zeta] = efdd(Y,N,fs,nsing)
%
%By Adarsh S, Ph.D. Candidate IIT Kanpur
%
%What this function does:
%------------------------
%This function calculates natural frequencies, modeshapes and damping by Frequency Domain Decomposition (FDD)
% 
%Input Arguments:
%---------------
% Y is the different channels of data arranged in rows
%Ig opt = 1, N is the number of points which are used for Fourier Transform while calculating spectral densities (should be a power of 2)  
%and if opt = 2 N is the number of lags required in the calculation of the cross-correlations (should be even)
%fs is the sampling frequencies
% nsing is the number of singular values to be plotted (1 to 4)
%If opt = 1, the spectral density matrix is calculated by FFT , and if opt = 2, the spectral density matrix is calculated from correlation matrix after applying an exponential window
%
%Output Arguments:
%-----------------
%wnat is the identified natural frequencies (rad/sec)
%Phi is the identified mode shapes
%zeta is the identified damping values
%
%Note
%-----
%nsing should not be greater than 4
%
%Ex: [wnat,Phi,zeta] = efdd(Y,1024,200,4,1) ;
%
%
%                                      -----------------------------------------------
%                            *********|| All rights reserved; Adarsh S; October, 2019 || *********
%                                      -----------------------------------------------
%
%
%***********************************************************************************************************************************************************
%===========================================================================================================================================================


if nsing>4
    error('Number of singular values to be displayed should be less than or equal to 4') ;
end

if opt == 1
    [Gyy,f] = specmatfft(Y,N,fs) ;
elseif opt == 2
    [Gyy, f] = specmatcormat(Y,N,1/fs) ;
else
    error('The value of opt should be 1 or 2') ;
end

clear Y ;
lf = length(f) ;
[mGyy,nGyy,oGyy] = size(Gyy) ;

U1col = zeros(mGyy,oGyy) ;
Singv = zeros(nsing,oGyy) ;

for i = 1:1:oGyy
    
    [U,S,~] = svd(Gyy(:,:,i)) ;
    U1col(:,i) = U(:,1) ;
    
    for j = 1:1:nsing
        Singv(j,i) = S(j,j) ;
    end
    
end

clear U S

plotCol = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.4940 0.1840 0.5560; 0.6350 0.0780 0.1840] ;

for i = 1:1:nsing
    plot(f,20*log10(Singv(i,:)),'Color',plotCol(i,:),'LineWidth',2);
    hold on
end
title('Singular Values vs Frequencies','fontsize',14)
xlabel('Frequency Hz','fontsize',14,'FontWeight','bold')
ylabel('Singular Values in dB','fontsize',14,'FontWeight','bold')

if nsing == 1
    legend('First Singular Value')
elseif nsing == 2
    legend('First Singular Value','Second Singular Value')
elseif nsing == 3
    legend('First Singular Value','Second Singular Value','Third Singular Value')
elseif nsing == 4
    legend('First Singular Value','Second Singular Value','Third Singular Value','Fourth Singular Value')
end

fprintf('************************************************************************************************************************\n\n')
fprintf('Please note down the number of modes visibe and their frequency ranges; after that please press any key to continue :-)\n\n')
fprintf('************************************************************************************************************************\n\n')
pause
nummod = input('Please enter the number of modes/significant peaks you see in the singular value plot ') ;

moddata = zeros(nummod,3) ;
filtf = zeros(2*nummod,lf) ;

for i = 1:1:nummod
    disp(['For Mode ',num2str(i)]) ;
    moddata(i,2) = input('Enter the lower bound frequency for mode ') ;
    moddata(i,3) = input('Enter the upper bound frequency for mode ') ;
    k = 0 ;
    for j = 1:1:lf
        if f(j)>=moddata(i,2) && f(j)<=moddata(i,3)
            k = k + 1 ;
            filtf(2*i-1,k) = j ;
            filtf(2*i,k) = f(j) ;
        end
    end
end


for i = 1:1:lf
    a = nnz(filtf(:,i)) ;
    if a == 0
        break ;
    end
end

filtf(:,i-1:lf) = [] ;
[mfiltf nfiltf] = size(filtf) ;
filtfsize = zeros(nummod,1) ;
maxfreq = zeros(nummod,1) ;
Phi = zeros(mGyy,nummod) ;

for i = 1:1:nummod
    
    filtfsize(i,1) = nnz(filtf(2*i,:)) ;
    filtfInd = filtf(2*i-1,1:filtfsize(i,1)) ;
    Singdum = Singv(1,filtfInd) ;
    [~,SingMaxI] = max(Singdum) ;
    maxfreq(i,1) = filtf(2*i,SingMaxI) ;
    Phi(:,i) = U1col(:,filtf(2*i-1,SingMaxI));
    
    
end

clear U1col

Gqq = zeros(nummod,oGyy) ;
PinvPhi = pinv(Phi) ;

for i = 1:1:nummod
    
    
    for j = filtf(2*i-1,1):1:filtf(2*i-1,filtfsize(i))
        
        Gqq(i,j) =  real(PinvPhi(i,:)*Gyy(:,:,j)*PinvPhi(i,:)') ;
        
    end
    
end

clear PinvPhi Gyy
GqqD = [Gqq,(fliplr(Gqq(:,2:oGyy-1)))] ;
clear Gqq
Rqq1 = ifft(GqqD.') ;
Rqq = Rqq1(1:oGyy,:).';
clear Rqq1
df = f(2) - f(1) ;
dt = 1/(df*N) ;
t = (0:1:length(f)-1)*dt ;
wnat = zeros(nummod,1) ;
zeta = zeros(nummod,1) ;

for i = 1:1:nummod
    [wnat(i,1),zeta(i,1),~,~] = compexft(Rqq(i,:),t,dt,1) ;
end
    
end