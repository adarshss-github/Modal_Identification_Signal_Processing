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


[m n] = size(compmat) ;

realmat = zeros(m,n) ;

for i = 1:1:n
   realmat(:,i) = comp2realmod(compmat(:,i)) ;
end

end
    