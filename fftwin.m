function [Y,f] = fftwin(y,fs,r,win)

%===========================================================================================================================================================
%***********************************************************************************************************************************************************
%
%                                      --------------------------------------------------
%                                                          PREAMBLE
%                                      --------------------------------------------------
%[Y,f] = fftwin(y,fs,r,win)
%
%By Adarsh S, Ph.D. Candidate IIT Kanpur
%This function calculates the double sided Discrete Fourier Tranform of a signal after windowing 
% and zero padding, to make the number of signal bits a power of 2
% y is the 1-D signal (column vector), fs is the sampling frequency of the signal, r is the
% windowing option (by default is r =1 i.e., Hann/Cosine Tapering  ), if win = 1 (by default applies window and zero-padding) else if win = 0 (applies only zero-padding and no window)
% Y is the DFT values, f is the frequency axis in Hz
%Ex: [Y,f] = fftwin(y,100) ;
%Note: The term 1/N is included while calculating DFT (by default not there in matlab) i.e., the FFT values are divided by number of fft points
%
%                                      -----------------------------------------------
%                            *********|| All rights reserved; Adarsh S; October, 2019 || *********
%                                      -----------------------------------------------
%
%
%***********************************************************************************************************************************************************
%===========================================================================================================================================================

if nargin == 2
    
    r = 1 ;
    win = 1 ;
    
end

if nargin == 3
    
    win = 1 ;
    
end

[m n] = size(y) ;
[mr] = length(r) ;
dt = 1/fs ;

if n~=1
    
    error('Signal should be a column vector') ;
    
end

if mr == 0
    
    r = 1 ;
    
end

if win~=0 && win~=1
    
    error('Value of win should be either 0 or 1') ;
    
end

if r~=0 && r~=1
    
    error('Value of r should be either 0 or 1') ;
    
end



if win == 1
    
    winfn = tukeywin(m,r) ;
    y = y.*winfn ;
    
end

numfft = 2^(nextpow2(m)) ;
df = 1/(dt*numfft) ;
f = df*(0:1:numfft-1) ;
Y = fft(y,numfft)/numfft ;

end