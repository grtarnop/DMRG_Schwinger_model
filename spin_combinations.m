% v2. 08/09/2024.
% This function finds all possible combinations of spin projections for a given string 
% N = [N(1), N(2), ..., N(p)] of dimensions of the Hamiltonian parts and 
% the total spin projection s
% The output is an array sct(k, p), where k is the number of combinations 
% of the spin projections. 
% The first column of spins sct(:, 1) is in acending order

function [ sct ] = spin_combinations(N, s)

%N = [2, 2, 2];
%s = 0;
p = size(N, 2);

% A loop to find all the possible spin projections q(i) = N(i)/2 + s(i) 
% and filter only thouse which satisfy q(1) + q(2) + ... + q(p) = q_tot
% The loop runs over i = 1 ,..., (N(1) + 1)*...*(N(p) + 1)
% and we convert i to the string (q(1), q(2), ..., q(p)) using that 
% i = (N(2) + 1)*...*(N(p) + 1)*q(1) + (N(3) + 1)*...*(N(p) + 1)*q(2) + ... 
%       + (N(p) + 1)*q(p - 1) + q(p) + 1

sct = [];               % initialize an array for possible spin combinations 
q = zeros(1, p);        % initialize a zero array for spin projections (charges)
q_tot = sum(N)/2 + s;   % total spin projection (charge) of the chain 
for i = 1:prod(N + 1) 
    % find the first element q(1)
    q(1) = floor((i - 1)/prod(N(2:end) + 1));
    % find the middle elements q(2), q(3),..., q(p - 1)
    for j = 2:(p - 1)
        sub = 0;
        for l = 1:(j - 1)
            sub = sub + prod(N(l + 1:end) + 1)*q(l);
        end
        q(j) = floor((i - 1 - sub)/prod(N(j + 1:end) + 1));
    end
    % find the last element q(m)
    sub = 0;
    for l = 1:(p - 1)
        sub = sub + prod(N(l + 1:end) + 1)*q(l);
    end
    q(p) = i - 1 - sub;  
    if sum(q) == q_tot 
        sct = horzcat(sct, (q - N/2)');
    end
end
sct = sct';
end
