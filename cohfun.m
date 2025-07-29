function [coh,f] = cohfun(op,inp,N,fs)

%===========================================================================================================================================================
%***********************************************************************************************************************************************************
%
%                                      --------------------------------------------------
%                                                          PREAMBLE
%                                      --------------------------------------------------
%
%[coh,f] = cohfun(op,inp,N,fs)
%By Adarsh S, Ph.D. Candidate IIT Kanpur
%
%What this function does
%-----------------------
%This function function calculates the coherence function from given output and input
%
%Input Arguments
%----------------
%op is the output time series
%inp is the input time series
%N is the number of points in the window used for fft (should be a power of 2) 
%fs is the sampling frequency of the signals
%
%Output Arguments
%----------------
%coh is the coherence function
%f is the frequency axis
%
%--------------------------------------------------------------------------------------------------------------------------
%Note: 
%1) op  and inp should be row vectors
%--------------------------------------------------------------------------------------------------------------------------
%
%Example: [coh,f] = cohfun(Yacc(1,:),inp,1024*16,2000) ;
%
%
%                                      -----------------------------------------------
%                            *********|| All rights reserved; Adarsh S; October, 2019 || *********
%                                      -----------------------------------------------
%
%
%***********************************************************************************************************************************************************
%===========================================================================================================================================================

[Sinpop,f] = specden(inp,op,N,fs) ;
[Sinin,~] = specden(inp,inp,N,fs) ;
[Sopop,~] = specden(op,op,N,fs) ;

coh = (abs(Sinpop).^2) ./ (Sinin.*Sopop) ;

clear Sinpop Sinin Sopop ;

end