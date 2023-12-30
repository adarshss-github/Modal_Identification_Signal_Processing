function [fud,zeta,Phi,Ad,Bd,C,D] = N4SID(Y,U,s,N,fs)

%===========================================================================================================================================================
%***********************************************************************************************************************************************************
%
%                                      --------------------------------------------------
%                                                          PREAMBLE
%                                      --------------------------------------------------
%
%By Adarsh S, Ph.D. Candidate IIT Kanpur
%[fud,zeta,Phi,Ad,Bd,C,D] = N4SID(Y,U,s,N,fs) ;
%
%Function description:
%-----------------------
%Does a modal and discrete state-space input-ouput identification by using
%N4SID subsapce identification algorithm [Combined deterministic-stochastic identification ]
%
%Input arguments:
%----------------
%1)Y: The measured output response channels arranged in rows
%2)U: Measured input channels arranged in rows
%3)s: Number of block rows of future or past quantities
%4)N: Number of significant modes
%5)fs: Sampling frequency
%
%Ouput arguments:
%----------------
%1)fud: The continuous time undmaped natural frequency in Hz 
%2)zeta: The continuous time damping ratio 
%3)Phi: The real mode shape matrix (undamped /classical damping assumption)
%4)Ad: Discrete state matrix
%5)Bd: Disccrete input matrix
%6)C: Output matrix
%7)D: Feedback matrix
%
% x[k+1] = Ad*x[k] + B*u[k] + v[k]
% y[k] = C*x[k] + D*u[k] + u[k]
%
%Example:
%-------
%Ex: [fud,zeta,Phi,Ad,Bd,C,D] = N4SID(Y,U,50,18,200) ;
%
%Notes:
%-----
%1)ZOH is used for discretization
%2)Approximates OX_f as the projection of the rowspace of Y_f along the
%rowspace of U_f (oblique projection) on to the rowspace of W_p, W_p = [U_p;y_p] 
%3)LQ tranformation is used for calculating the oblique projections
%4) The identified modal parameters appear in pairs due to complex conjugate
%poles
%5)Don't confuse with n4sid (lower case) MATLAB built-in function
%6)N4SID stands for 'Numerical algorithms for subspace state space sysytem identification; this acronym should be read like
% the Californian license plate: 'enforce it'. :-) 
%
%Reference:
%----------
% Subspace Identification for Linear Systems by Peter VAN
%OVERSCHEE and Bart De MOOR. Page:124
%
%
%                                      --------------------------------------
%                            *********|| Copyright; Adarsh S; October, 2019 || *********
%                                      --------------------------------------
%
%
%***********************************************************************************************************************************************************
%===========================================================================================================================================================

[mY,~] = size(Y) ; %mY is the number of output channels
[mU,~] = size(U) ; %mU is the number of input channels
N = 2*N ; % 2*N is the number of singular values of the projection matrix to be retained down the line

%===========================================================================================================================================================

                                          %Computation of P_i, Gamm_i and X_i
                                         %------------------------------------

                                       

HU = bHankel(U,2*s) ; %Block Hankel of inputs
HY = bHankel(Y,2*s) ;  %Block Hankel of outputs
Ui = HU(s*mU+1:s*mU+mU,:) ; %Used later for estimationg C and D matrices
Yi = HY(s*mY+1:s*mY+mY,:) ; %Used later for estimationg C and D matrices
H1 = [HU(s*mU+1:2*s*mU,:); HU(1:s*mU,:);HY ;] ; % H1 = [U_f;U_p;Y_p;Y_f]
Up = HU(1:s*mU,:); %Used for oblique projection later
clear HU ;
Yp = HY(1:s*mY,:) ; %Used for oblique projection later
clear HY ;

[~,L1] = qr(H1',0) ; % Lower triangular matrix
L1 = L1' ;
clear H1
L122 = L1(s*mU+1:s*(2*mU+mY),s*mU+1:s*(2*mU+mY)) ;
L132 = L1(s*(2*mU+mY)+1:2*s*(mU+mY),s*mU+1:s*(2*mU+mY)) ;
clear L1
Pi = L132*(L122\[Up;Yp]) ; % First oblique projection
clear L132 L122 Up Yp

[Usvd,Sig,~] = svd(Pi,'econ') ; 
Gami = Usvd(:,1:N)*sqrt(Sig(1:N,1:N)) ;
clear Usvd Sig 
Xi = Gami\Pi ;
clear  Pi

%===========================================================================================================================================================
                                         
                                               %Computation of P_ip, Gamm_ip and X_ip
                                              %---------------------------------------


HU = bHankel(U,2*s) ;  %Block Hankel of inputs
HY = bHankel(Y,2*s) ;  %Block Hankel of outputs
clear Y U
Upp = HU(1:s*mU+mU,:);  %Used for oblique projection later
Ufm = HU(s*mU+mU+1:2*s*mU,:);
clear HU ;
Ypp = HY(1:s*mY+mY,:) ;  %Used for oblique projection later
Yfm = HY(s*mY+mY+1:2*s*mY,:) ;
clear HY
H2 = [Ufm;Upp;Ypp;Yfm] ;  % H2 = [U_fm;U_pp;Y_pp;Y_fm]
clear Ufm Yfm

[~,L2] = qr(H2',0) ;  % Lower triangular matrix
L2 = L2' ;
clear H2
L222 = L2(mU*(s-1)+1:s*(2*mU+mY)+mY,mU*(s-1)+1:s*(2*mU+mY)+mY) ;
L232 = L2(s*(2*mU+mY)+mY+1:2*s*(mU+mY),mU*(s-1)+1:s*(2*mU+mY)+mY) ;
clear L2
Pip =  L232*(L222\[Upp;Ypp]) ; % Second oblique projection
clear L232 L222 Upp Ypp

Gamip = Gami(1:mY*(s-1),:) ;
clear Usvd Sig 
Xip = Gamip\Pip ;
clear Gamip Gami Pip

%===========================================================================================================================================================
 
                                          %Computation of discrete system matrices (A_d, B_d, C, D)
                                        %----------------------------------------------------------

ABCD = [Xip;Yi]/[Xi;Ui] ; %[A_d B_d ; C D]
clear Xip Yi X1 Ui
Ad = ABCD(1:N,1:N) ;
Bd = ABCD(1:N,N+1:end) ;
C = ABCD(N+1:end,1:N) ;
D = ABCD(N+1:end,N+1:end) ;
clear ABCD

%--------------
%| Easter egg |: Finding my wife, Vidya, was the greatest thing that ever happened in my life.
%--------------

%===========================================================================================================================================================

                                             %Computation of modal parameters
                                            %---------------------------------

[Psi,Mu] = eig(Ad) ;
[Psi,Mu] = stablepoles(Psi,Mu,1/fs) ; % Only stable poles and eigen vectors kept
[fud,zeta,~] = dispol2contmod(Mu,1/fs) ;
clear Mu

[Phicomp] = C*Psi ;
Phi = comp2realmodmat(Phicomp) ;
clear Psi Phicomp

[~,I] = sort(fud) ;
fud = fud(I) ;
[fud,fudUI] = unique(fud) ;
zeta = zeta(I) ;
zeta = zeta(fudUI) ;
Phi = Phi(:,I) ;
Phi = Phi(:,fudUI) ;
clear I fudUI

end

%===========================================================================================================================================================

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

%===========================================================================================================================================================

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

%===========================================================================================================================================================
    
function realmat = comp2realmodmat(compmat)

%===========================================================================================================================================================
%***********************************************************************************************************************************************************
%
%                                      --------------------------------------------------
%                                                          PREAMBLE
%                                      --------------------------------------------------
%
%
%realm = comp2realmodmat(compm)
%
%By Adarsh S, Ph.D. Candidate IIT Kanpur
%
%Description
%-----------
% Converts the identified complex mode shape matrix to real mode shape
% matrix
%
%Input arguments
%---------------
%1)compmat is the matrix representing the complex mode shapes in columns
%
%Output arguments
%----------------
%1)realmat is the matrix representing the converted real mode shapes in
%columns
%
%Ex: realmat = comp2realmodmat(compmat) ;
%
%                                      -----------------------------------------------
%                            *********|| All rights reserved; Adarsh S; October, 2019 || *********
%                                      -----------------------------------------------
%
%
%***********************************************************************************************************************************************************
%===========================================================================================================================================================


[m,n] = size(compmat) ;

realmat = zeros(m,n) ;

for i = 1:1:n
   realmat(:,i) = comp2realmod(compmat(:,i)) ;
end

end

%===========================================================================================================================================================

function realm = comp2realmod(compm)

%===========================================================================================================================================================
%***********************************************************************************************************************************************************
%
%                                      --------------------------------------------------
%                                                          PREAMBLE
%                                      --------------------------------------------------
%
%
%realm = comp2realmod(compm)
%
%By Adarsh S, Ph.D. Candidate IIT Kanpur
%
%Description
%-----------
% Converts the identified complex modes to real modes
%
%Input arguments
%---------------
%1)compm is the column vector representing the complex mode shape
%
%Output arguments
%----------------
%1)realm is the column vector representing the converted real mode shape
%
%Ex: realm = comp2realmod(compm) ;
%
%                                      -----------------------------------------------
%                            *********|| All rights reserved; Adarsh S; October, 2019 || *********
%                                      -----------------------------------------------
%
%
%***********************************************************************************************************************************************************
%===========================================================================================================================================================



[m,n] = size(compm) ;

if n~=1
    error('Complex mode should be a column vector') ;
end

realm = complex(zeros(m,1)) ;
imagP = zeros(m,1) ;
realP = zeros(m,1) ;
mag = zeros(m,1) ;
rotcompm = zeros(m,1) ;
dum = zeros(2,1) ;
phang = zeros(m,1) ;

realP = real(compm) ;
imagP = imag(compm) ;
p = polyfit(realP,imagP,1) ;

rotang = -atan(p(1)) ;

RotTrMat = [cos(rotang) -sin(rotang) ; sin(rotang)  cos(rotang)]; 

for i = 1:1:m
    dum = RotTrMat*[realP(i,1);imagP(i,1)] ;
    rotcompm(i,1) = dum(1) + dum(2)*1i ;
end
 
clear realP imagP compm ;

mag = abs(rotcompm) ;
mag = mag/max(mag) ;

phang = angle(rotcompm) ;

s = 0 ;

for i = 1:1:m
    
    if mag(i,1)~= 0
        if phang(i,1)>=0 && phang(i,1)<pi/2
            s = 1 ;
        elseif phang(i,1)<0 && phang(i,1)>= -pi/2
            s = 1 ;
        elseif phang(i,1)>=pi/2 && phang(i,1)<=pi
            s = -1 ;
        elseif  phang(i,1)>=-pi && phang(i,1)<-pi/2
            s = -1 ;
        end
    else
        
        s = 0 ;
        
    end
  
  realm(i,1) = s*mag(i,1) ;
  
  
    
end
clear phang mag ;
end

%===========================================================================================================================================================
