function [fud,zeta,muc] = dispol2contmod(Mu,dt)

%===========================================================================================================================================================
%***********************************************************************************************************************************************************
%
%                                      --------------------------------------------------
%                                                          PREAMBLE
%                                      --------------------------------------------------
%
%By Adarsh S, Ph.D. Candidate IIT Kanpur
%[fud,zeta,muc] = dispol2contmod(Mu,dt)
%
%What this function does:
%-----------------------
%%Converts discrete time poles to continuous time modal parameters
%
%Input arguments:
%----------------
%Mu is the diagonal eigen value matrix of the discrete time state matrix
%dt is the sampling time
%
%Ouput arguments:
%----------------
%fud is the continuous time undmaped natural frequency in Hz (will be repeated due to the complex conjugate pole pairs in Mu)
%zeta is the continuous time damping ratio (will be repeated due to the complex conjugate pole pairs in Mu)
%
%
%Example:
%-------
%Ex: [fud,zeta] = dispol2contmod(Mu,1/200) ;
%
%
%                                      -----------------------------------------------
%                            *********|| All rights reserved; Adarsh S; October, 2019 || *********
%                                      -----------------------------------------------
%
%
%***********************************************************************************************************************************************************
%===========================================================================================================================================================

mud = diag(Mu) ;
n = length(mud) ;
muc = zeros(n,1) ;
fud = zeros(n,1) ;
zeta =  zeros(n,1) ;
wud = 0 ;

for i = 1:1:n
    
    muc(i) = log(mud(i))/dt ;
     wud = abs(muc(i)) ;
     fud(i) = wud/(2*pi) ;
     zeta(i) = -real(muc(i))/(wud) ;
    
end

clear mud ;

end

    
