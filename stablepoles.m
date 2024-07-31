function [PsiSt,MuSt] = stablepoles(Psi,Mu,dt)

%===========================================================================================================================================================
%***********************************************************************************************************************************************************
%
%                                      --------------------------------------------------
%                                                          PREAMBLE
%                                      --------------------------------------------------
%
%
%
%By Adarsh S, Ph.D. Candidate IIT Kanpur
% This function returns the stable poles and the corresponding eigen vectors
%Psi is the matrix containing conjugate eigen vector pairs of the discrete time state space matrix
% Mu is the matrix containing the conjugate eigen value paisrs of the discrete  time state space matrix
% dt is the sampling time
%PsiSt is the matrix containing conjuagte eigen vector pairs corresponding to stable poles
%MuSt is the matrix containing the stable poles


%                                      -----------------------------------------------
%                            *********|| All rights reserved; Adarsh S; October, 2019 || *********
%                                      -----------------------------------------------


%***********************************************************************************************************************************************************
%===========================================================================================================================================================


zetamax = 1 ;
zetamin = 0 ;

[~,zeta] = dispol2contmod(Mu,dt) ;
n = length(zeta) ;
PsiSt = zeros(n,n) ;
MuSt = zeros(n,n) ;
k = 1 ;

for i = 1:1:n
    
    if zeta(i)>zetamin && zeta(i)<zetamax
        MuSt(k,k) = Mu(i,i) ;
        PsiSt(:,k) = Psi(:,i) ;
        k = k + 1 ;
        end
end

k = 0 ;

for i = 1:1:n
    
    if MuSt(i,i) == 0
        k = i-1 ;
        break 
    end
    
end

if i == n
    k = n ;
end

MuSt = MuSt(1:k,1:k) ;
PsiSt = PsiSt(:,1:k) ;

end

%===========================================================================================================================================================