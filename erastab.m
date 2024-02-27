function [modn,freqn,Phin,zetan] = erastab(R,s,dt,modn)

%=========================================================================================================================

[m,n,o] = size(R) ; %Size of 3-D correlation matrix 


%=========================================================================================================================

H = bHankel(R(:,:,1:o-1),s) ; % Block Hankel matrix
Hp = bHankel(R(:,:,2:o),s) ;  % Shifted block Hankel matrix
clear R
[U,Sig,V] = svd(H) ;
clear H
%=================================================================================================================================

lmodn = length(modn) ;
freqn = cell(lmodn,1) ;
Phin = cell(lmodn,1) ;
zetan = cell(lmodn,1) ;

for i = 1:1:lmodn
    Ob = U(:,1:modn(i))*sqrt(Sig(1:modn(i),1:modn(i))) ;  %Extended Observability matrix
    Co =  sqrt(Sig(1:modn(i),1:modn(i)))*(V(:,1:modn(i)))';
    Ad = pinv(Ob)*Hp*pinv(Co) ; %Discrete state matrix
    Cd = Ob(1:m,:) ;            %Discrete observation matrix
    [Psi,Mu] = eig(Ad) ; %Mu contains the discrete poles 
    [Psi,Mu] = stablepoles(Psi,Mu,dt) ; %Only stable poles and corresponding eigen vectors are kept
    Phin{i,1} = Cd*Psi ; %Complex mode shapes corresponding to the measured locations
    [fud,zeta,~] = dispol2contmod(Mu,dt);
    freqn{i,1} = fud ;
 zetan{i,1} = zeta ;
  
end

clear fud zeta Ob Co Ad Cd Psi Mu ;

for i = 1:1:lmodn
    fud  =  freqn{i,1} ;
    n = modn(i)*ones(length(fud),1) ;
    scatter(fud,n,'d') ;
    hold on
end

title('Stabilization Plot for Natural Frequency ','fontsize',14)
xlabel('Frequency Hz','fontsize',14,'FontWeight','bold')
ylabel('Modal Order','fontsize',14,'FontWeight','bold')
set(gcf,'color','w');

grid on
end