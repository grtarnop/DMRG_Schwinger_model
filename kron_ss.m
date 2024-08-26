% v2. 08/09/2024.
% This function computes a tensor product of two operators A and B:
% H{s1, s2} = A{sl1, sl2} x B{sr1, sr2} for a given spin projections s1 and s2 
% using the Generalized Kronecker-product rule.
% The tensor product picks only combinations of spin projections such that 
% sl1 + sr1 = s1 and sl2 + sr2 = s2 where sl and sr are    
% spin projections of A and B
% sl1, sl2 are in an acsending order and sr1, sr2 are in a descending order 
% The operators A and B are given by cells 
% A{ql1, ql2} and B{qr1, qr2}, where q = N/2 + s + 1  = 1, 2 ,..., 2s+1
% and ql1, ql2 and qr1, qr2 go over all the spin projections of the given blocks
% N = [N1, N2] is the number of spin-1/2 sites of A and B blocks

function [ H ] = kron_ss(A, B, q1, q2, N)

ss1 = spin_combinations(N, q1 - sum(N)/2 - 1); % combinations of spins for s1 = q1 - (N1+N2)/2 - 1
ss2 = spin_combinations(N, q2 - sum(N)/2 - 1); % combinations of spins for s2 = q2 - (N1+N2)/2 - 1
k1 = size(ss1, 1);                        % number of combinations of spins for s1
k2 = size(ss2, 1);                        % number of combinations of spins for s2


Haux = cell(k1, k2);
% loop over all combinations of the spin projections
for j1 = 1:k1
    for j2 = 1:k2
     Haux{j1, j2} = kron(A{ss1(j1, 1) + N(1)/2 + 1, ss2(j2, 1) + N(1)/2 + 1}, ... 
                         B{ss1(j1, 2) + N(2)/2 + 1, ss2(j2, 2) + N(2)/2 + 1});
    end
end
H = cell2mat(Haux);

end