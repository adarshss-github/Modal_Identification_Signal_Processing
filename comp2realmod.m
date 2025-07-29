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