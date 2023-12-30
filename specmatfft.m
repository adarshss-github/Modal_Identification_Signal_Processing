function [Gyy f] = specmatfft(Y,N,fs)

%===========================================================================================================================================================
%***********************************************************************************************************************************************************
%
%                                      --------------------------------------------------
%                                                          PREAMBLE
%                                      --------------------------------------------------
%[Gyy f] = specmatfft(Y,N,fs)
%
% By Adarsh S, Ph.D. Candidate IIT Kanpur
%
%What this function does
%----------------------
%Calculates the spectral density matrix (one sided) from the measurements
%
%Input arguments
%----------------
%Y is the measurements arranged in rows
%N is the number of points in the window used for fft (should be a power of 2) 
%fs is the sampling rate of the measurements
%
%Ouput arguments
%----------------
%Gyy is the one sided cross specttal density matrix (3D matrix)
%f is the frequency axis (Hz)
%
%
%Ex:[Gyy f] = specmatfft(Y,4*1024,100) ;
%                                      -----------------------------------------------
%                            *********|| All rights reserved; Adarsh S; October, 2019 || *********
%                                      -----------------------------------------------


%***********************************************************************************************************************************************************
%===========================================================================================================================================================



[m n] = size(Y) ;

[~,f] = specden(Y(1,:),Y(1,:),N,fs) ; ;
nf = length(f) ;
Gyy = zeros(m,m,nf) ;


for i = 1:1:m
    for j = 1:1:m
        [Gyy(i,j,:),~] = specden(Y(i,:),Y(j,:),N,fs) ; 
    end
end

      
end
