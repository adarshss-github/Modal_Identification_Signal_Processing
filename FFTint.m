function [Yfilt,t,Hf,f] = FFTint(Y,fs,N,f1,f2,Intn)

%===========================================================================================================================================================
%***********************************************************************************************************************************************************
%
%                                      --------------------------------------------------
%                                                          PREAMBLE
%                                      --------------------------------------------------
%
%By Adarsh S, Ph.D. Candidate IIT Kanpur
%[Yfilt,t,Hf,f] = FFTint(Y,fs,N,f1,f2,Intn)
%
%Function description:
%-----------------------
%Can integrate and differentiate signals using FFT filtering
%
%
%Input arguments:
%----------------
%1) Y: The measured output response channels arranged in rows
%2) fs: Sampling frequency
%3) N: the size of the segments (should be a power of 2)
%4) f1 is the frequency at which filter transition for cut-off begins (Hz)
%5) f2 is the frequency at which filter transition ends (Hz)
%6) Intn: Number for integration 
%
%Ouput arguments:
%----------------
%1) Yfilt: The integrated or diffrentiated data
%2) t: The time axis
%3)Hf: The FRF of the filter
%4) f: The frequency axis (Hz)
%
%Example:
%-------
%Ex: [Yfilt,t,Hf,f] = FFTint(Y,200,2048,0,0.8,1) ;
%
%Note:
%-----
%1) If n = 1, a single integration is performed; if n = 2, a double
%integration is performed; if n = -1, a single diffrentiation is performed
%2) For differntiation, if the filtered signal becomes complex, consider
%only the real part
%3) Segmentation is useful for very long signals in where the frequency
%content varies with time
%
%
%
%                                      -----------------------------------------------
%                            *********|| Copyright; Adarsh S; October, 2019 || *********
%                                      -----------------------------------------------
%
%
%***********************************************************************************************************************************************************
%===========================================================================================================================================================

%===========================================================================================================================================================



[Yseg,f,seg] = fftsegdat(Y,N,fs) ;
lf = length(f) ;
nY = size(Y,2) ;
[m,n,p] = size(Yseg) ;
fdum = zeros(2,lf) ;
k = 1 ;

for j = 1:1:lf 
    if f(j)>=f1 & f(j)<=f2
        fdum(1,k) = f(j) ;
        fdum(2,k) = j ;
        k = k + 1 ;
    end
end



fdum(:,k:end) = [] ;
lfdum = length(fdum(1,:)) ;
win = tukeywin(2*lfdum,1) ;

if Intn>=0

Hf1 = [zeros(1, ( fdum(2,1)-1 ) ) win(1:lfdum)' ones(1, lf - (fdum(2,end)+1) + 1 )] ;
Hf2 = [0 (i*2*pi*f(2:end)).^(-Intn)] ;
Hf = Hf1.*Hf2 ;

else

Hf1 = [zeros(1, ( fdum(2,1)-1 ) ) win(1:lfdum)' ones(1, lf - (fdum(2,end)+1) + 1 )] ;
Hf2 = [(i*2*pi*f(1:end)).^(-Intn)] ;
Hf = Hf1.*Hf2 ;

end

for h = 1:1:m
    for j = 1:1:p
      Yseg(h,:,j) =   Yseg(h,:,j).*Hf ;
    end
end

[Yfilt,t] = fftreconsdat(Yseg,nY,seg,fs) ;

end