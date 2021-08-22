function [ W ] = SparsityW(X,lambda)
%get similarity by sparse representation
% % %  Coded by Chang Lu (lucifer1992@email.swu.edu.cn), College of Computer and Information Science,
% % %   Southwest University.
% % %   version 1.0 date:2016-08-07
%%%%%%%%%%
%   X: a protein-term matrix
%   lambda: a scalar regularization parameter that balances the tradeoff between reconstruction error and sparsity of coefficients
%%%%%%%%%%
if(matlabpool('size')<5)
 matlabpool close force;
matlabpool local 5;
end
%  if isempty(gcp('nocreate'))
%   delete(gcp('nocreate'));
%   parpool(20);
%  end
[D,N]=size(X);
W1=zeros(N,N-1);
W=zeros(N,N);
opts=[];
    parfor kk=1:N
        Xk=X(:,kk);
        tempX=X(:,setdiff(1:N,kk));
        A=tempX;
        [xp,funval]= nnLeastR(A,Xk,lambda,opts);
        xp=abs(xp);
         W1(kk,:)=xp(1:N-1);
    end  
    for ii=1:N
         W(ii,setdiff(1:N,ii))=W1(ii,:);
    end 
%         delete(gcp('nocreate')); 
            matlabpool close; 
    end

