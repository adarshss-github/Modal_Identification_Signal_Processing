function [Y12,f] =  specden(y1,y2,N,fs)

%===========================================================================================================================================================
%***********************************************************************************************************************************************************
%
%                                      --------------------------------------------------
%                                                          PREAMBLE
%                                      --------------------------------------------------
%
% By Adarsh S, Ph.D. Candidate IIT Kanpur
%[Y12,f] =  specden(y1,y2,N,fs)
%
%What this function does
%-----------------------
%Estimates the one sided cross spectral densities of two signals y1 and y2 by Welch method (uses 50 % overlapping Hann windows)
%
%Input arguments
%----------------
%y1 and y2 are the two signals
%N is the number of points in the window used for fft (should be a power of 2) 
%fs is the sampling rate of the measurements
%N is the number of points for the window (Should be a power of 2)
%fs is the sampling frequency of the signals
%
%Ouput arguments
%----------------
%Y12 is the half cross spectral density and f is the frequency axis 
%f is the frequency axis (Hz)
%
%Note
%------
%y1 and y2 should be row vectors
%If y1 and y2 are same signal, it estimates the power spectral density
%%To get Y12 in (unit)^2/Hz, divide Y12 by fs ; ||***||
%
%Example
%-------
%Ex: [Y12,f] =  specden(y1,y2,1024,200) ;
%
%                                      -----------------------------------------------
%                            *********|| All rights reserved; Adarsh S; October, 2019 || *********
%                                      -----------------------------------------------
%
%
%***********************************************************************************************************************************************************
%===========================================================================================================================================================


dt = 1/fs ;
df = 1/(N*dt) ;
[m1 n1] = size(y1) ;
[m2 n2] = size(y2) ;

r = 2^(nextpow2(N)) ;

%============================================================================================================================

if r~=N
    
    error('N should be a power of 2');
    
end

if m1~=1 || m2~=1
    
    error('Signals not row vectors') ;
    
end

if n1~=n2
    
    error('Signals not of same size') ;
    
end

%============================================================================================================================

n = n1 ;

if (n - fix(n/N)*N + 1) >= 0.5*N
    
    nseg = fix(n/N) + fix(n/N) ;
    seg = N*0.5:N*0.5:(N*fix(n/N)) ;
    
else
    
    nseg = fix(n/N) + fix(n/N) - 1  ;
    seg = N*0.5:N*0.5:(N*fix(n/N)-N*0.5) ;
    
end

Y12dum = zeros(1,N) ;

for i = 1:1:nseg
    
    [Y1,~] = fftwin(y1(1,(seg(i)-N*0.5 + 1):(seg(i)+ N*0.5 ))',fs,1) ;
    Y1 = Y1.'*N ;
    [Y2,~] = fftwin(y2(1,(seg(i)-N*0.5 + 1 ):(seg(i)+ N*0.5))',fs,1) ;
    Y2 = Y2.'*N ;
    Y12dum = Y12dum + conj(Y1).*Y2 ;
    
end

winfn = tukeywin(N,1) ;

 [~,f] = fftwin(y1(seg(1)-N*0.5+1:seg(1)+ N*0.5)',fs,1) ;
 f = f(1:N/2+1) ;

Y12 = Y12dum(1:N/2+1)/(nseg*(winfn')*(winfn)) ; %Window correction

%============================================================================================================================

end
    
    

