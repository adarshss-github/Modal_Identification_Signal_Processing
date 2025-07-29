function [Ynoise] = AddNois(Y,noiseqty,opt,fs,ti,tf)

%===========================================================================================================================================================
%***********************************************************************************************************************************************************
%
%                                      --------------------------------------------------
%                                                          PREAMBLE
%                                      --------------------------------------------------
%
%[[Ynoise] = AddNois(Y,noiseqty,opt,fs,ti,tf)
%
%By Adarsh S, Ph.D. Candidate IIT Kanpur
%
%Function description:
%------------------------
%Adds Gaussian white noise to respones in Y 
%
%%Input Arguments:
%----------------
%1) Y: Responses arranged in rows
%2) noiseqty : Column vector containing the power sectral density values or rms of noise to be added to each channel
%3) opt : If opt = 1, noiseqty should contain power sectral density values, and if opt = 2,
%noiseqty should contain rms
%4) fs: Sampling frequency used (see function, gausswn)
%5) ti: Initial time (see function, gausswn)
%6) tf: Final time (see function, gausswn)
%
%Output Arguments:
%-----------------
%1) Ynoise: TResponse with noise added
%
%Example: [Ynoise] = errAddRmsPer(Y,[10 10],100,0,200) ; ;
%
%                                      -----------------------------------------------
%                            *********|| All rights reserved; Adarsh S; October, 2019 || *********
%                                      -----------------------------------------------
%
%

%===========================================================================================================================================================
%***********************************************************************************************************************************************************

[mY,nY] = size(Y) ;
dumerr = zeros(1,nY) ;
Ynoise = zeros(mY,nY) ;

if opt == 1
    
    for i = 1:1:mY
        
        [dumerr,~] = gausswn(noiseqty(i),1,fs,ti,tf,1) ;
        Ynoise(i,:) = Y(i,:) + dumerr ;
        
    end
    
elseif opt == 2
    
    for i = 1:1:mY
        
        [dumerr,~] = gausswn(noiseqty(i),2,fs,ti,tf,1) ;
        Ynoise(i,:) = Y(i,:) + dumerr ;
        
    end
    
end
    

clear dumerr


end