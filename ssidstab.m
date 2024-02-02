function [modn,freqn] = ssidstab(Y,s,dt,modn)



%=================================================================================================================================
                                          %Computation of Hankel matrix of past and future observations
                                          %------------------------------------------------------------                                          

[m n] = size(Y) ; % m is the numeber of channels of data and n is the total number of data points
H = bHankel(Y,2*s) ; %Hankel matrix (top half is past and bottom half is future observations)
clear Y ;
[~,nc] = size(H) ; %Number of columns of the Hankel matrix

%=================================================================================================================================
                                               % LQ factorization of the Hankel matrix  
                                               %----------------------------------
                                          
L = triu(qr(H'))' ; % Start of LQ factorization of Hankel matrix (Reduces storage space and the number of flops while calculating the projection) 
clear H ;
L21 = L((m*s+1:m*s+ m),(1:m*s)) ;
L31 = L((m*s+m+1:2*m*s),(1:m*s)) ;
clear L ;

%=================================================================================================================================
                                                   %Projection step
                                                   %--------------
                                                   
O = [L21;L31] ; %Projection of rowspace of future observation matrix to the row space of past observation matrix (Main Theorem of SSID), i.e, E(Yf/Yp) %***
clear L21 L31 ;
[U,Sig,V] = svd(O,'econ') ;  %(W1OW2 = O since the algorithm is unweighted pricnipal component )

%=================================================================================================================================

lmodn = length(modn) ;
freqn = cell(lmodn,1) ;

for i = 1:1:lmodn
    Ob = U(:,1:modn(i))*sqrt(Sig(1:modn(i),1:modn(i))) ;
    [mOb,nOb] = size(Ob) ;
    Adis = pinv(Ob(1:mOb-m,:))*Ob(m+1:mOb,:) ;
    [Psi,Mu] = eig(Adis) ;
    [Psi,Mu] = stablepoles(Psi,Mu,dt) ;
    [fud,~] = dispol2contmod(Mu,1/200);
    freqn{i,1} = fud ;
end

clear O U Sig V Ob ;

for i = 1:1:lmodn
    fud  =  freqn{i,1} ;
    n = modn(i)*ones(length(fud),1) ;
    scatter(fud,n,'+') ;
    hold on
end
 
grid on

end
    













