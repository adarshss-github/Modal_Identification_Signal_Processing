function [wnat,zeta,res,hrecov] = compexft(h,t,dt,Nmod)

%===========================================================================================================================================================
%***********************************************************************************************************************************************************
%
%                                      --------------------------------------------------
%                                                          PREAMBLE
%                                      --------------------------------------------------
%
%
%By Adarsh S, Ph.D. Candidate IIT Kanpur
%
%This function calculates the modal parameters from the given IRF by Complex Exponential Fitting 
%h is the Impulse Response Function (row vector) 
% t is the time corresponding to h in sec
% dt is the sampling time
%wnat is the identified natural frequencies (rad/sec)
%zeta is the identified damping ratios
%res are the values of the residues
%hrecov is the fitted IRF
%Ex:  [wnat,zeta,res,hrecov] = compexft(h,t,0.01,1) ;
%
%                                      -----------------------------------------------
%                            *********|| All rights reserved; Adarsh S; October, 2019 || *********
%                                      -----------------------------------------------
%
%
%***********************************************************************************************************************************************************
%===========================================================================================================================================================

[mh,nh] = size(h) ;

if mh ~= 1 
    
    error('h should be a row vector') ;

end

hHank = zeros(nh-2*Nmod,2*Nmod) ;
hr = zeros(nh-2*Nmod,1) ;

%===========================================================================================================================================================

for i = 1:1:(nh-2*Nmod)
     init = i ;
     for j = 1:1:2*Nmod
         hHank(i,j) = h(init) ;
         init = init + 1 ;
     end
     hr(i) = -h(init) ;
end

%===========================================================================================================================================================

beta = hHank\hr ;
beta = [1; beta(end:-1:1)] ;
Vr = roots(beta) ;
poldum = log(Vr)./dt ;
pol = zeros(Nmod,1) ;
j = 1 ;

for i = 1:2:2*Nmod
    pol(j) = poldum(i) ;
    j = j + 1 ;
end

wnat = abs(pol) ;
zeta = -real(pol)./(wnat) ;

resmat = zeros(2*Nmod,2*Nmod) ;
resh = zeros(2*Nmod,1) ;

for i = 1:1:2*Nmod
      
      %resmat(i,:) = (Vr.^(i-1))' ;
      resmat(i,:) = (Vr.^(t(i)/dt))' ;
      resh(i) = h(i) ;
end

res = resmat\resh ;

hrecov = zeros(2*Nmod,1) ;

for i = 1:1:nh
    
    hrecov(i) = ((Vr.^(t(i)/dt))')*res ;
    
end

hrecov = real(hrecov) ;

%===========================================================================================================================================================

end