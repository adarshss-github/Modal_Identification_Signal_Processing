function [Gyy f] = specmatcormat(Y,N,dt)

%===========================================================================================================================================================
%***********************************************************************************************************************************************************
%
%                                      --------------------------------------------------
%                                                          PREAMBLE
%                                      --------------------------------------------------
%[Gyy f] = specmatcormat(Y,N,dt)
%
% By Adarsh S, Ph.D. Candidate IIT Kanpur
%
%What this function does
%----------------------
%Computes the one sided cross spectral density matrix from the one sided cross correlation matrix, R, (after the application of exponential windows)
%
%Input arguments
%----------------
%1)Y is the matrix containing the recorded ouputs arranged in rows
%N is the number of points in the window used for fft (should be a power of 2) 
%2)N is the number of lags required in the calculation of the cross (should be even)
%correlation matrix
%3)dt is the sampling time
%
%Ouput arguments
%----------------
%1)Gyy is the computed cross spectral density matrix, 
%2)f is the corresponding frequency axis
%
%Ex:[Gyy f] = specmatcormat(Y,8194,1/200) ;
%
%                                      -----------------------------------------------
%                            *********|| All rights reserved; Adarsh S; October, 2019 || *********
%                                      -----------------------------------------------


%***********************************************************************************************************************************************************
%===========================================================================================================================================================



[R,Tau] = corrmat(Y,N,1,dt) ;
[W] = expwin(Tau,0.01,1) ;
Rw = corrmatexp(R,Tau,0.01,1) ;
[Gyy f] = corrmat2specmat(Rw,Tau) ;


end
