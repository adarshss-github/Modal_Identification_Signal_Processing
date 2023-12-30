function [Psiact Mu1] = era(R,s,N,dt)

%===========================================================================================================================================================
%***********************************************************************************************************************************************************
%
%                                      --------------------------------------------------
%                                                          PREAMBLE
%                                      --------------------------------------------------
%[Psiact Mu] = era(R,s,N,dt)
%
%By Adarsh S, Ph.D. Candidate IIT Kanpur
%
%What this function does
%-----------------------
%Identifies the system modal parameters by Eigen Realization Algorithm
%
%Input Arguments
%---------------
% R is the 3-D correlation matrix of the measured responses
%s is the number of block rows required for the construction of the Hankel
%matrix
%N is the number of modes to be identified
%dt is the sampling time of the measurements
%
%Output Arguments
%-----------------
%Psiact is the complex mode shapes corresponding to the measured locations
%Mu is the discrete poles 
%
%Ex:[Psiact Mu] = era(R,200,2,1/200) ; 
%
%                                      -----------------------------------------------
%                            *********|| All rights reserved; Adarsh S; October, 2019 || *********
%                                      -----------------------------------------------
%
%
%***********************************************************************************************************************************************************
%===========================================================================================================================================================



%=========================================================================================================================

[m,n,o] = size(R) ; %Size of 3-D correlation matrix 
N = 2*N ; %Size of state vector

%=========================================================================================================================

H = bHankel(R(:,:,1:o-1),s) ; % Block Hankel matrix
Hp = bHankel(R(:,:,2:o),s) ;  % Shifted block Hankel matrix
clear R
[U,Sig,V] = svd(H) ;
Ob = U(:,1:N)*sqrt(Sig(1:N,1:N)) ;  %Extended Observability matrix
Co =  sqrt(Sig(1:N,1:N))*(V(:,1:N))';

%=========================================================================================================================

Ad = pinv(Ob)*Hp*pinv(Co) ; %Discrete state matrix
Cd = Ob(1:m,:) ;            %Discrete observation matrix

%=========================================================================================================================

[Psi,Mu] = eig(Ad) ; %Mu contains the discrete time poles 
[Psi1,Mu1] = stablepoles(Psi,Mu,dt) ; %Only stable poles and corresponding eigen vectors are kept
Psiact = Cd*Psi1 ; %Complex mode shapes corresponding to the measured locations

%=========================================================================================================================

end
