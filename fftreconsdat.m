function [Ytime,t] = fftreconsdat(Yseg,nY,seg,fs)

%===========================================================================================================================================================
%***********************************************************************************************************************************************************
%
%                                      --------------------------------------------------
%                                                          PREAMBLE
%                                      --------------------------------------------------
%
%
%By Adarsh S, Ph.D. Candidate IIT Kanpur
%Constructs back the time domain signal from the segmented fft values (see fftsegdat.m)
%Yseg contains the segmented data in frequency domain
% seg is the segment midpoints from fftsegdat
%fs is the sampling frequency
%nY is the total data points in the measurements
%Ytime is the reconstructed signal
%t is the time axis


%                                      -----------------------------------------------
%                            *********|| All rights reserved; Adarsh S; October, 2019 || *********
%                                      -----------------------------------------------


%***********************************************************************************************************************************************************
%===========================================================================================================================================================

[m n p] = size(Yseg) ;
nfft = 2*(n-1) ;
dt = 1/fs ;

Ytime = zeros(m,nY) ;
Ydum = zeros(1,n) ;
t = (0:1:nY-1)*dt ; 


for i = 1:1:m
    for j = 1:1:p
        
        Ydum = ifft( [Yseg(i,:,j) (flipud( (conj(Yseg(i,2:n-1,j) )).' )).' ] ) ;
       Ytime(i,(seg(j)-(n-1) + 1):(seg(j)+ n-1 )) = Ytime(i,(seg(j)-(n-1) + 1):(seg(j)+ n-1 )) + Ydum ;
      
    end
    Ydum = zeros(1,n) ;
end

%===========================================================================================================================================================

end
        
    