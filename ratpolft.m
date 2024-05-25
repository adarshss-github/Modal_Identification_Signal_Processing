function [wnat,zeta,peaks,Hrecov] = ratpolft(w,H,Nmod,Nnum,Nden)

%===========================================================================================================================================================
%***********************************************************************************************************************************************************
%
%                                      --------------------------------------------------
%                                                          PREAMBLE
%                                      --------------------------------------------------
%
%[wnat,zeta,peaks,Hrecov] = ratpolft(w,H,Nmod,Nnum,Nden)
%
%By Adarsh S, Ph.D. Candidate IIT Kanpur
%
%What this function does:
%-----------------------
%This function calculates the modal parameters from a noisy FRF by using Rational Polynomial Curve Fitting (SIMO version)  
% Forsythe algorithm was used to evaluate the values of the orthogonal polynomials 
%
%Input Arguments:
%----------------
% w is the column vector containing the frequencies in rad/sec
% H is the column vector containing the values of the measured Receptance
%Nmod is the number of modes to be fitted
%Nnum is the degree of the numerator polynomial (optional)
%Nden is the degree of the denominator polynomial (optional)
%
%Output Arguments:
%----------------
%wnat are the identified natural frequencies (rad/sec)
%zeta are the identified damping ratios
%peaks are the the values of the fitted Receptance FRF at the peaks
%Hrecov is the fitted Receptance FRF
%
%--------------------------------------------------------------------------------------------------------------------------
%Note: 
%1) Please make sure that DC component is not there in the FRF measurements
%2) This is intented for Receptance FRF (when Nnum and Nden are not specified explicitely in the function ; for other types of FRF specify Nnum and Nden accordingly)
%--------------------------------------------------------------------------------------------------------------------------
%
%Ex: [wnat,zeta,peaks] = ratpolft(w(200:530),H(200:530,2),3) ;
%Ex: [wnat,zeta,peaks] = ratpolft(w(200:530),H(200:530,2),3,5,6) ;
%
%References:
%-----------
%Main Formulation and Forsythe Algorithm
%1)Richardson, M.H. and Formenti, D.L., 1982, November. Parameter estimation from frequency response measurements using rational fraction polynomials. In Proceedings of the 1st international modal analysis conference (Vol. 1, pp. 167-186). Union College Schenectady, NY.
%Forsythe polynomials and conversion of orthogonal polynomial coefficients to normal polynomial coefficients
%2)Kelly, L.G., 1967. Handbook of numerical methods and applications. Reading, Mass.: Addison-Wesley, 1967
%3)Cristian Guti�rrez Acu�a (2019). Rational Fraction Polynomial Method (https://www.mathworks.com/matlabcentral/fileexchange/3805-rational-fraction-polynomial-method), MATLAB Central File Exchange. Retrieved October 26, 2019.
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

if nargin == 3             %Default for Receptance
    Nnum = 2*Nmod-1 ;
    Nden = 2*Nmod ;
end

[NH,MH] = size(H) ;
Hrecov = zeros(NH,MH) ;
[Nw,Mw] = size(w) ;
wmax = max(w) ;

if MH~=1 || Mw~=1
    
    error('Both FRF and frequencies vectors should be column vectors') ;
    
end

if NH~=Nw
    
    error('Both FRF and frequencies vectors sholud be of same size') ;
    
end

%===========================================================================================================================

[Phi,PhiTr] = forspol(w,H,1,Nnum) ; %Forsythe algorithm for numerator
[Theta,ThetaTr] = forspol(w,H,2,Nden) ; %Forsythe algorithm for denominator

T = sparse(diag(H))*Theta(:,1:end-1) ;
W = H.*Theta(:,end) ;
clear w H
Xp = -2*real(Phi'*T) ;
Gp = 2*real(Phi'*W) ;

%===========================================================================================================================================================

LHS = [eye(Nnum + 1) Xp ; Xp.' eye(Nden)] ;
RHS = [Gp; zeros(Nden,1)] ;
cd = LHS\RHS ;
c = cd(1:Nnum+1);
d = cd(Nnum+2:end) ;
d1 = [d;1] ;
clear LHS RHS cd T W Xp Gp ;

%===========================================================================================================================================================

% MXp = size(Xp,1) ;
% d = -inv(eye(MXp)- Xp.'*Xp)*Xp.'*Gp ;
% c = Gp - Xp*d ;
% d1 = [d;1] ;

%===========================================================================================================================================================

for i = 1:1:NH
    
    Hrecovnum = sum(c.'.*Phi(i,:)) ;
    Hrecovden = sum(d1.'.*Theta(i,:)) ;
    Hrecov(i) =   Hrecovnum/Hrecovden ;
    
end

%============================================================================================================================

a = PhiTr*c ;
b = ThetaTr*d1 ;

%============================================================================================================================
poldum = roots(b(end:-1:1))*wmax ;
pol = zeros(Nmod,1) ;
j = 1 ;

for i = 1:2:2*Nmod
    pol(j) = poldum(i) ;
    j = j + 1 ;
end

clear poldum

wnat = abs(pol) ;
zeta = -real(pol)./(wnat) ;
peaks = frfrecon((wnat/wmax)',a,b) ;

%==========================================================================================================================

end

function [P,Transmat] = forspol(w,H,opt,Npol) 

%===========================================================================================================================================================
%***********************************************************************************************************************************************************
%
%                                      --------------------------------------------------
%                                                          PREAMBLE
%                                      --------------------------------------------------
%
%
%By Adarsh S ,Ph.D. Candidate IIT Kanpur
%
%This function calculates values of the orthogonal polynomials by using the
%Forsythe recursive algorithm 
% It also returns the transformation matrix for converting the coefficients of the orthogonal polynomials to that of the normal coefficients 
%
% w is the row vector containing the frequencies in rad/sec
% H is the row vector containing the values of the measured Receptance
%opt = 1 for numerator polynomials and opt = 2 for denominator polynomials
%Npol is the maximum degree of the polynomials
%P is the matrix containing the evaluated orthogonal plynomials
%Transmat is the transformation matrix for converting the coefficients of the orthogonal polynomials to that of the normal coefficients 

%References:
%Generation and evaluation of Forsythe polynomials
%1)Richardson, M.H. and Formenti, D.L., 1982, November. Parameter estimation from frequency response measurements using rational fraction polynomials. In Proceedings of the 1st international modal analysis conference (Vol. 1, pp. 167-186). Union College Schenectady, NY.
%Finding the Transformation matrix
%2)Kelly, L.G., 1967. Handbook of numerical methods and applications. Reading, Mass.: Addison-Wesley, 1967.
%3)Cristian Guti�rrez Acu�a (2019). Rational Fraction Polynomial Method (https://www.mathworks.com/matlabcentral/fileexchange/3805-rational-fraction-polynomial-method), MATLAB Central File Exchange. Retrieved October 26, 2019.

%Ex:  [P,Transmat] = forspol(w,H,opt,Npol) 

%                                      -----------------------------------------------
%                            *********|| All rights reserved; Adarsh S; October, 2019 || *********
%                                      -----------------------------------------------


%***********************************************************************************************************************************************************
%===========================================================================================================================================================

wmax = max(w) ;
w = w./wmax ;
Nw = length(w) ;
R = zeros(Nw,Npol+2) ;
P = zeros(Nw,Npol+1) ;
Transmat = zeros(Npol+1,Npol+2);

if opt == 1
    
    q = ones(Nw,1) ;
    
elseif opt == 2
    
    q = (abs(H)).^2 ;
    
else
    
    error('The value of opt must be 1 or 2') ;
    
end

%===========================================================================================================================================================

R(:,1:2) = [zeros(Nw,1), ( 1/sqrt(2*sum(q)) ).*ones(Nw,1)] ;
Transmat(1,2) = 1/sqrt(2*sum(q)) ;

for k = 1 : Npol
    
    Vkm1 = 2*( sum( w.*R(:,k+1).*R(:,k).*q ) ) ;
    Sk = w.*R(:,k+1) - Vkm1.*R(:,k) ;
    Dk = sqrt( 2*sum( (Sk.^2).*q )  ) ;
    R(:,k+2) = Sk./Dk ;
    Transmat(:,k+2) = -Vkm1*Transmat(:,k) ;
	Transmat(2:k+1,k+2) = Transmat(2:k+1,k+2) + Transmat(1:k,k+1) ;
    Transmat(:,k+2) = Transmat(:,k+2)/Dk ;
    
end

R = R(:,2:Npol+2) ;
Transmat = Transmat(:,2:Npol+2) ;

%===========================================================================================================================================================

compstor = zeros(1,Npol+1) ;

for k = 0 : Npol
    
    P(:,k+1) = R(:,k+1)*(1i)^k ;
    compstor(1,k+1) = (1i)^k ;
    
end
 
Transmat = (compstor'*compstor).*Transmat ;  
    
%===========================================================================================================================================================    

end

function [H] = frfrecon(w,a,b)

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
%This function reconstructs the FRF values from the coefficients of the numerator and denominator polynomials
%This can be used to increase the resolution of the FRF
%
%a is the coefficients of numerator polynomials (lower to higher order)
%b is the coefficients of denominator polynomials (lower to higher order)
% w is the normalized frequencies (rad/sec) at which the FRF is to be evaluated 
% H are the values of the FRF evaluated at the given frequencies
%
%                                      -----------------------------------------------
%                            *********|| All rights reserved; Adarsh S; October, 2019 || *********
%                                      -----------------------------------------------
%
%
%***********************************************************************************************************************************************************
%===========================================================================================================================================================

Na = length(a) ;
Nb = length(b) ;
Nw = length(w) ;
num = zeros(Nw,1) ;
den = zeros(Nw,1) ;

%===========================================================================================================================================================

for i = 1:1:Nw
    for j = 1:1:Na
        num(i) = num(i) + a(j)*(1i*w(i))^(j-1) ;
    end
end

for i = 1:1:Nw
    for j = 1:1:Nb
        den(i) = den(i) + b(j)*(1i*w(i))^(j-1) ;
    end
end

H = num./den ;

%===========================================================================================================================================================

end

