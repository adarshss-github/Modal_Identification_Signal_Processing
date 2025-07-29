function [Yseg,f,seg] = fftsegdat(Y,N,fs)

%===========================================================================================================================================================
%***********************************************************************************************************************************************************
%
%                                      --------------------------------------------------
%                                                          PREAMBLE
%                                      --------------------------------------------------
%
%
%By Adarsh S, Ph.D. Candidate IIT Kanpur
%Segments the data and calculates the FFT (half spectrum) of each segment after applying
%Hann window; the segments are overlapping ones with 50% overlap
%Y contains the measured data in row wise form
%N is the size of the segments (sholud be a power of 2)
%fs is the sampling frequency
%Yseg is  3-D matrix containing the FFT values of the segments
%f is the frequency axis
%seg is the midpoint array values of the segments
%Ex: [Yseg,f,seg] = fftsegdat(Y,8192,200)

%                                      -----------------------------------------------
%                            *********|| All rights reserved; Adarsh S; October, 2019 || *********
%                                      -----------------------------------------------


%***********************************************************************************************************************************************************
%===========================================================================================================================================================



dt = 1/fs ;
df = 1/(N*dt) ;
[m n] = size(Y) ;


r = rem(N,2) ;

%============================================================================================================================

if r~0
    
    error('N should be a power of 2');
    
end

%============================================================================================================================


if (n - fix(n/N)*N + 1) >= 0.5*N
    
    nseg = fix(n/N) + fix(n/N) ;
    seg = N*0.5:N*0.5:(N*fix(n/N)) ;
    
else
    
    nseg = fix(n/N) + fix(n/N) - 1  ;
    seg = N*0.5:N*0.5:(N*fix(n/N)-N*0.5) ;
    
end

    
    

Yseg = zeros(m,N/2+1,nseg) ;
Ydum = zeros(N,1) ;
dum = 0 ;

for i = 1:1:m
    
 for j = 1:1:nseg
    
     if dum == 0
       [Ydum,f] = fftwin(Y(i,(seg(j)-N*0.5 + 1):(seg(j)+ N*0.5 ))',fs,1) ;
         dum = dum + 1 ;
        end
         
         
    [Ydum,~] = fftwin(Y(i,(seg(j)-N*0.5 + 1):(seg(j)+ N*0.5 ))',fs,1) ;
    Ydum = Ydum.'*N ;
    Yseg(i,:,j) = Ydum(1:N/2+1) ;
     
     end
 
end
 
 
f = f(1:N/2+1) ;

end
 
 



