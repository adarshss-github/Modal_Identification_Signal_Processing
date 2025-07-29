function [Yfilt,t,Hf,f] = hpfilt(Y,fs,N,f1,f2)

%===========================================================================================================================================================
%***********************************************************************************************************************************************************
%
%                                      --------------------------------------------------
%                                                          PREAMBLE
%                                      --------------------------------------------------
%
%
%By Adarsh S, Ph.D. Candidate IIT Kanpur
%Performs high pass filtering on data
%Y is the measurements in row wise manner
%N is the size of the segments (should be a power of 2)
%fs is the sampling frequency (Hz)
%f1 is the frequency at which filter transition for cut-off begins (Hz)
%f2 is the frequency at which filter transition ends (Hz)
%Yfilt is the filtered data
%t is the time axis
%Hf is the FRF of the filter
% f is the frequency axis (Hz)
%Ex: [Yfilt,t] = hpfilt(Y,200,8192,2,3.5) ; Here the flat band will be
%greater than 3.5 Hz
%Note: Make sure that the transition is smooth and N is high enough

%                                               -----------------------------------------------
%                                     *********|| All rights reserved; Adarsh S; October, 2019 || *********
%                                               -----------------------------------------------


%***********************************************************************************************************************************************************
%===========================================================================================================================================================



[Yseg,f,seg] = fftsegdat(Y,N,fs) ;
lf = length(f) ;
nY = size(Y,2) ;
[m,n,p] = size(Yseg) ;
fdum = zeros(2,lf) ;
k = 1 ;

for i = 1:1:lf 
    if f(i)>=f1 & f(i)<=f2
        fdum(1,k) = f(i) ;
        fdum(2,k) = i ;
        k = k + 1 ;
    end
end



fdum(:,k:end) = [] ;
lfdum = length(fdum(1,:)) ;
win = tukeywin(2*lfdum,1) ;

Hf = [zeros(1, ( fdum(2,1)-1 ) ) win(1:lfdum)' ones(1, lf - (fdum(2,end)+1) + 1 )] ;

for i = 1:1:m
    for j = 1:1:p
      Yseg(i,:,j) =   Yseg(i,:,j).*Hf ;
    end
end

[Yfilt,t] = fftreconsdat(Yseg,nY,seg,fs) ;

end