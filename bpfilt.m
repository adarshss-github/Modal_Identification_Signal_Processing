function [Yfilt,t,Hf,f] = bpfilt(Y,fs,N,f1,f2,f3,f4)

%===========================================================================================================================================================
%***********************************************************************************************************************************************************
%
%                                      --------------------------------------------------
%                                                          PREAMBLE
%                                      --------------------------------------------------
%
%
%By Adarsh S, Ph.D. Candidate IIT Kanpur
%Performs band pass filtering on data
%Y is the measurements arranged in row wise manner
%N is the size of the segments (should be a power of 2)
%fs is the sampling frequency of the signal in Hz
%f1 is the frequency at which filter transition for high pass cut-off begins
%(Hz)
%f2 is the frequency at which filter transition for high pass ends (Hz)
%f3 is the frequency at which filter transition for low pass cut-off
%begins (Hz)
%f2 is the frequency at which filter transition for low pass ends (Hz)
%Yfilt is the filtered data
%t is the time axis
%Hf is the FRF of the filter
% f is the frequency axis (Hz)
%Ex: [Yfilt,t] = bpfilt(Y,200,8192,0.8,1,1.3,1.5) ; Here the flat band will
%be between 1 and 1.3 Hz
%Note: Make sure that the transition is smooth and N is high enough
%
%                                                 ------------------------------
%                                        *********|| Copyright 2019 Adarsh S  || *********
%                                                -------------------------------
%
%
%***********************************************************************************************************************************************************
%===========================================================================================================================================================

%===========================================================================================================================================================


[Yseg,f,seg] = fftsegdat(Y,N,fs) ;
lf = length(f) ;
nY = size(Y,2) ;
[m,n,p] = size(Yseg) ;
fdum1 = zeros(2,lf) ;
fdum2 = zeros(2,lf) ; 

k = 1 ;
q = 1 ;

for i = 1:1:lf 
    if f(i)>=f1 & f(i)<=f2
        fdum1(1,k) = f(i) ;
        fdum1(2,k) = i ;
        k = k + 1 ;
    end
end



fdum1(:,k:end) = [] ;
lfdum1 = length(fdum1(1,:)) ;
win1 = tukeywin(2*lfdum1,1) ;

for i = 1:1:lf 
    if f(i)>=f3 & f(i)<=f4
        fdum2(1,q) = f(i) ;
        fdum2(2,q) = i ;
        q = q + 1 ;
    end
end



fdum2(:,q:end) = [] ;
lfdum2 = length(fdum2(1,:)) ;
win2 = tukeywin(2*lfdum2,1) ;

Hf = [zeros(1, ( fdum1(2,1)-1 ) ) win1(1:lfdum1)' ones( 1,( (fdum2(2,1)-1) - (fdum1(2,end)+1) + 1 ) ) win2(lfdum2+1:end)'  zeros(1, (lf - (fdum2(2,end)+1) + 1 )     ) ] ;

for i = 1:1:m
    for j = 1:1:p
      Yseg(i,:,j) =   Yseg(i,:,j).*Hf ;
    end
end

[Yfilt,t] = fftreconsdat(Yseg,nY,seg,fs) ;

end

%===========================================================================================================================================================

