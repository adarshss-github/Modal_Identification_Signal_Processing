function [R,Tau] = corrmat(Y,N,side,dt)

%===========================================================================================================================================================
%***********************************************************************************************************************************************************
%
%                                      --------------------------------------------------
%                                                          PREAMBLE
%                                      --------------------------------------------------
%
%[R,Tau] = corrmat(Y,N,side,dt)
%
%By Adarsh S, Ph.D. Candidate IIT Kanpur
%
%%What this function does:
%------------------------
%Calculates the unbiased correlation matrix of responses by temporal averaging
%
%%Input Arguments:
%----------------
%Y is the matrix containing row wise arranged data
%N is the total number of lags in correlation matrix, including positive and negative time lags, N should be even
%Give side = 1 for 1 sided correlation matrix and side = 2 for double sided correlation matrix 
%dt is the sampling time of the measurements
%
%Output Arguments:
%-----------------
%R is the correlation matrix and Tau is the lag axis
%
%Example [corrmatunb,Tau] = corrmat(Y,16384,2,1/200) ;
%
%                                      -----------------------------------------------
%                            *********|| All rights reserved; Adarsh S; October, 2019 || *********
%                                      -----------------------------------------------
%
%
%***********************************************************************************************************************************************************
%===========================================================================================================================================================


if mod(N,2) ~= 0
  error('N is not even') ;
end

[m,n] = size(Y) ;
% m is the number of responses  amd n is the total discrete time points

if side == 1
   R = zeros(m,m,N*0.5 + 1) ;
   for i = 1:1:N*0.5+1
      Y1 = Y(:,1:n-(i-1)) ;
      Y2 = Y(:,1+(i-1):n) ;
      R(:,:,i) = (1/(n-i))*Y1*Y2' ;
      clear Y1 Y2 ;
   end
   Tau = (0:1:N*0.5)*dt ;
elseif side == 2
    R = zeros(m,m,N + 1) ;
    k = N*0.5 + 1 ;
    for i = 1:1:N*0.5+1
       Y1 = Y(:,1:n-(i-1)) ;
       Y2 = Y(:,1+(i-1):n) ;
       R(:,:,k) = (1/(n-i))*Y1*Y2' ;
       clear Y1 Y2 ;
       if k~= N*0.5 + 1
          l = N + 2 -k ;
          R(:,:,l) = R(:,:,k) ;
       end
       k = k + 1 ;
    end
    Tau = (-N*0.5:1:N*0.5)*dt ;
else
    error('The value of side should be either 1 or 2') ;
end

end









