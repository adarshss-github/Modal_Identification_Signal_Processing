function [Yds,tds] = dwnsa(Y,fs,N,Nseg,f1,f2)

%===========================================================================================================================================================
%***********************************************************************************************************************************************************
%
%                                      --------------------------------------------------
%                                                          PREAMBLE
%                                      --------------------------------------------------
%
%[Yds,tds] = dwnsa(Y,fs,N,Nseg,f1,f2)
%By Adarsh S, Ph.D. Candidate IIT Kanpur
%
%This function downsamples the measurements to fd = fs/N times (N should be an integer); the signals are filtered using a low-pass filter first to prevent aliasing 
% The signal is low-pass filtered with a cut-off of approximately fd/2
%Y contains the measurements in rows
%fs is the original sampling frequency
%N is the downsmapling factor
%Nseg is the number of segments used for filering (2^N)
%f1 is the frequency at which filter transition for cut-off begins (Hz)
%f2 is the frequency at which filter transition ends (Hz)
%Yds is the downsampled measurements
%tds is the downsampled time axis
%Ex:   [Yds,tds] = dwnsa(Yacc,2000,10,1024*16,95,100)  ;
%
%                                      -----------------------------------------------
%                            *********|| All rights reserved; Adarsh S; October, 2019 || *********
%                                      -----------------------------------------------
%
%
%***********************************************************************************************************************************************************
%===========================================================================================================================================================

[Yfilt,~,~,~] = lpfilt(Y,fs,Nseg,f1,f2) ;
clear Y
[mYfilt nYfilt] = size(Yfilt) ;
NYds = floor(nYfilt/N) ;
Yds = zeros(mYfilt,NYds) ;

s = 1:N:nYfilt ;
tds = (s-1)*(1/fs) ;
Yds = Yfilt(:,s) ;
clear s Yfilt

end
