function [Singv,f] = SingVal(Y,N,fs,nsing,opt)

%===========================================================================================================================================================
%***********************************************************************************************************************************************************
%
%                                      --------------------------------------------------
%                                                          PREAMBLE
%                                      --------------------------------------------------
%
%[Singv,f] = SingVal(Y,N,fs,nsing,opt)
%
%By Adarsh S, Ph.D. Candidate IIT Kanpur
%
%What this function does:
%------------------------
%This function calculates and plots the singular values of the response spectral density matrix
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
%----------------
%Singv is the matrix of singular values arranged in rows 
%f is the frequency axis
%
%Note:
%-----
%nsing should not be greater than 4
%
%Ex: [Singv,f] = SingVal(Yacc,1024*16,2000,4,2) ;
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
    error('Number of singular values to be displayed should be less than 4')
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
xlabel('Frequency [Hz]','fontsize',14,'FontWeight','bold')
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

end
