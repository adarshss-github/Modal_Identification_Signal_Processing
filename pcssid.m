function [Psiact,Mu1,Ad,Cd,Res] = pcssid(Y,s,N,dt)

%===========================================================================================================================================================
%***********************************************************************************************************************************************************
%
%                                      --------------------------------------------------
%                                                          PREAMBLE
%                                      --------------------------------------------------
%
%[Psi1,Mu1,Ad,Cd,Res] = pcssid(Y,s,N,dt)
%
%By Adarsh S, Ph.D. Candidate IIT Kanpur
%
%Description:
%-----------
% This is a program which performs the Principal Component Data-Driven Stochastic Subspace Identification algorithm
% 
%Input Arguments:
%---------------
% Y is the matrix containing row wise arranged data
% s is the number of block rows of future or past observed outputs
% N is the number of modes to be found out
% dt is the sampling time
%
%Output Arguments:
%-----------------
%Psiact consists of the complex mode shapes found corresponding to the measured quantities
% Mu1 consists of the discrete poles found (stable ones)
%Ad is the identified discrete state space matrix
%Cd is the identified observation matrix
%Res is the calculated residue matrix which can be used for characterizing noise
%
%Note
%-----
%
%
%Ex: [Psi1,Mu1,Ad,Cd,Res] = pcssid(Y,200,4,1/2000) ;
%
%
%                                      -----------------------------------------------
%                            *********|| All rights reserved; Adarsh S; October, 2019 || *********
%                                      -----------------------------------------------
%
%
%***********************************************************************************************************************************************************
%===========================================================================================================================================================

%=================================================================================================================================
                                          %Computation of Hankel matrix of past and future observations
                                          %------------------------------------------------------------                                          
[m n] = size(Y) ; % m is the number of channels of data and n is the total number of data points
N = 2*N ; %Total number of poles
H = bHankel(Y,2*s) ; %Hankel matrix (top half is past and bottom half is future observations)
clear Y ;
[~,nc] = size(H) ; %Number of columns of the Hankel matrix

%=================================================================================================================================
                                               % LQ factorization of the Hankel matrix  
                                               %--------------------------------------
L = triu(qr(H'))' ;% Start of LQ factorization of Hankel matrix (Reduces storage space and the number of flops while calculating the projection) 
clear H ;
L = L((1:2*m*s),(1:m*s+m)) ;
L11 = L((1:m*s),(1:m*s)) ;
L21 = L((m*s+1:m*s+ m),(1:m*s)) ;
L22 = L((m*s+1:m*s+m),(m*s+1:end)) ;
L31 = L((m*s+m+1:end),(1:m*s)) ;
L32 = L((m*s+m+1:2*m*s),(m*s+1:end)) ;
clear L ;

%=================================================================================================================================
                                                   %Projection step
                                                   %--------------
                                                   
W1OW2 = [L21;L31]*L11' ; %Projection of row-space of future observation matrix to the row space of past observation matrix (Main Theorem of SSID), i.e, E(Yf/Yp) multiplied by weighing matrices
[U,Sig,V] = svd(W1OW2,'econ') ;  

%=================================================================================================================================
                         
                                           %Calculation of System Matrices
                                           %------------------------------
                                           
Ob = U(:,1:N)*sqrt(Sig(1:N,1:N)) ; %Extended observability matrix
X = pinv(Ob)*[L21;L31] ; %Kalman state sequence 
Obp = U(1:m*s-m,1:N)*sqrt(Sig(1:N,1:N)) ;%Extended observability matrix without last m rows
Xp = pinv(Obp)*[L31 L32] ; %Kalman state sequence (plus one) 
Ys = [L21 L22] ; %Output sequence at s+1 th bock row 
clear W1OW2 U Sig V L11 L21 L22 L31 L32 Ob Obp ;
AC = [Xp;Ys]/[X zeros(N,m)] ; %Least squares solution
Ad = AC(1:N,1:end) ; %Discrete state space matrix
Cd = AC(N+1:end,1:end) ; %Discrete observation matrix 
Res = [Xp;Ys] - AC*[X zeros(N,m)] ; %Residues
clear AC ;
[Psi,Mu] = eig(Ad) ;
[Psi1,Mu1] = stablepoles(Psi,Mu,dt) ; % Only stable poles and Eigen vectors kept
Psiact = Cd*Psi1 ;


%=================================================================================================================================


end


function bH = bHankel(Y,r)

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
%Forms the block Hankel natrix with r block rows from response/correlation matrix Y (uses all data in Y)
% Y is the matrix containing row wise arranged data or Y can be the 3-D correlation matrix of the measured responses
% r is the number of block rows of the Hankel matrix
%bH is the block Hankel matrix formed with 'r' block rows (the entire data has been used)
%Example  bH = bHankel(Y,200) ;
%
%
%
%
%                                      -----------------------------------------------
%                            *********|| All rights reserved; Adarsh S; October, 2019 || *********
%                                      -----------------------------------------------
%
%
%***********************************************************************************************************************************************************
%===========================================================================================================================================================


[m,n,o] = size(Y) ;


if o == 1

bH = zeros(r*m,n-r+1) ;

for i = 1:1:r
    
   bH(1+m*(i-1):m*i,:) =  Y(:,i:n-r+i) ;
   
end

else

bH = zeros(r*m,(o-r+1)*n) ;

for i = 1:1:r
    
    for j = 1:1:o-r+1
        
        bH((i-1)*m+1:i*m,(j-1)*n+1:j*n) = Y(:,:,i+j-1) ;
    end
    
end

end

end

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



