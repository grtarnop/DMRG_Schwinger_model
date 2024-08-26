% This function generates a set of states with the spin projection s
% and chain size N
function [ vec ] = subspace_fixed_Sz(N, s)

% compute all possible configurations of 0 and 1
% with  q  ones and N - q zeros and make a set of vectors vec
q = N/2 + s; 
vec = cell(nchoosek(N,q),1);           % initialize a cell for vectors
count = 1;
   for k = 0:2^N-1
       state = (dec2bin(k,N)-'0');  % represent number k as a binary of length N
       if sum(state) == q           % check if state has charge q (i.e. q ones )
            vec{count} = 1;
                for i = 1:N
                    if state(i) == 1 
                        vec{count} = sparse(kron(vec{count},[1;0]));
                    else
                        vec{count} = sparse(kron(vec{count},[0;1]));  
                    end
                end    
       count = count + 1;
       end
   end
% combine all the vectors into a big matrix 2^N x C(N, N/2 + s)
% so each column has length 2^N and there are  C(N, N/2 + s) columns
vec = sparse([vec{:}]);

end