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
