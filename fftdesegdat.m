function [Ytime,t] = fftdesegdat(Yseg,nY,seg,fs)

[m n p] = size(Yseg) ;
nfft = 2*(n-1) ;
dt = 1/fs ;

Ytime = zeros(m,nY) ;
Ydum = zeros(1,n) ;
t = (0:1:nY-1)*dt ;


for i = 1:1:m
    for j = 1:1:p
        Ydum = ifft( [Yseg(i,:,j)] [flipud( conj(Yseg(i,2:n-1,j) ))] ) ;
        Ytime = Ytime(i,(seg(j)-2*n*0.5 + 1):(seg(j)+ 2*n*0.5 )) + Ydum ;
    end
end

end
        
    