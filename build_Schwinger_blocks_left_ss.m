% v1. 11/12/2023
% v3. 05/14/2024  (included theta angle)
% v4. 08/15/2024  (added spin projections domains and parallelization)
% This function builds left blocks in spin-spin sectors for the Schwinger model
% Input:
%   n: is the number of spin-1/2 sites in the final left blocks
%   N: is the total number of spin-1/2 sites in the left block at the end of
%       infinite DMRG procedure
% the output is cells by (n + 1) x (n + 1) size (all spin sectors)
function [ HL_new, QL_new, SL_new, SpL_new, IdL_new, TL_new ] = build_Schwinger_blocks_left_ss(HLss, ...
                                  QLss, SLss, SpLss, IdLss, H2ss, Q2ss, S2ss, M2ss, ...
                                  HXY2ss, Sp2Lss, Sp2Rss, I4ss, n, x, mu, y, theta, s, N)

% spin projections domain (represented in the form of integer indices)
if s >= 0
    if (n >= 2) && (n <= N - s)
        sdn = [-n/2: n/2] + n/2 + 1;
    elseif (n > N - s) && (n <= N + s)
        sdn = [n/2 - N + s: n/2] + n/2 + 1;
    elseif (n > N + s) && (n <= 2*N - 2)
        sdn = [n/2 - N + s : N + s - n/2] + n/2 + 1;
    end
elseif s < 0
    if (n >= 2) && (n <= N - abs(s))
        sdn = [-n/2: n/2] + n/2 + 1;
    elseif (n > N - abs(s)) && (n <= N + abs(s))
        sdn = [-n/2: N - abs(s) - n/2] + n/2 + 1;
    elseif (n > N + abs(s)) && (n <= 2*N - 2)
        sdn = [n/2 - N - abs(s): N - abs(s) - n/2] + n/2 + 1;
    end    
end

% Pre-allocate the cells for better performance
HL_new = cell(n + 1, n + 1);
QL_new = cell(n + 1, n + 1);
SL_new = cell(n + 1, n + 1);
SpL_new = cell(n + 1, n + 1);
IdL_new = cell(n + 1, n + 1);

% Pre-allocate the arrays to store intermediate results
HL_temp = cell((n + 1)^2, 1);
QL_temp = cell((n + 1)^2, 1);
SL_temp = cell((n + 1)^2, 1);
SpL_temp = cell((n + 1)^2, 1);
IdL_temp = cell((n + 1)^2, 1);

% Total number of iterations
total_iterations = (n + 1) * (n + 1);

% define auxiliarly operators                              
Sm2Rss = cellfun(@ctranspose,Sp2Rss','UniformOutput',false);
SmLss = cellfun(@ctranspose, SpLss','UniformOutput',false);
QLSqss = cellfun( @(x) x.^2 + (1/2 + theta/pi)*x, QLss, 'UniformOutput', false);

parfor idx = 1:total_iterations
% Calculate the corresponding q1 and q2 values
    [q1, q2] = ind2sub([n + 1, n + 1], idx);
    % check that q1 and q2 in the domain
    if ismember(q1, sdn) && ismember(q2, sdn)  
     % Compute the results and store them in the temporary arrays
        % build new HL   
        HL_temp{idx} = sparse(kron_ss(HLss, I4ss, q1, q2, [n - 2, 2]) ...
         + x*(kron_ss(SpLss, Sm2Rss, q1, q2, [n - 2, 2]) ...
              + kron_ss(SmLss, Sp2Rss, q1, q2, [n - 2, 2])...
              + kron_ss(IdLss, HXY2ss, q1, q2, [n - 2, 2])) ...
         + (mu/2)*kron_ss(IdLss, M2ss, q1, q2, [n - 2, 2]) ...
         + y*(kron_ss(IdLss, H2ss, q1, q2, [n - 2, 2]) ... 
              + kron_ss(QLss, S2ss, q1, q2, [n - 2, 2]) ...
              + 2*kron_ss(QLSqss, I4ss, q1, q2, [n - 2, 2]) ...
              + (theta^2/(2*pi)^2)*kron_ss(IdLss, I4ss, q1, q2, [n - 2, 2])));            
        % build new QL 
        QL_temp{idx} = sparse(kron_ss(QLss, I4ss, q1, q2, [n - 2, 2]) ... 
                            + kron_ss(IdLss, Q2ss, q1, q2, [n - 2, 2]));
                        
                     
        % build new SL
        SL_temp{idx} = sparse(kron_ss(SLss, I4ss, q1, q2, [n - 2, 2])... 
                         + kron_ss(IdLss, S2ss, q1, q2, [n - 2, 2]) ...
                         + 4*kron_ss(QLss, I4ss, q1, q2, [n - 2, 2]));
        % build new SpL and IdL
        SpL_temp{idx} = sparse(kron_ss(IdLss, Sp2Lss, q1, q2, [n - 2, 2]));
        IdL_temp{idx} = sparse(kron_ss(IdLss, I4ss, q1, q2, [n - 2, 2]));
    end
end

% Reconstruct the 2D cell arrays from the linearized temporary arrays
for idx = 1:total_iterations
    [q1, q2] = ind2sub([n + 1, n + 1], idx);
    % check that q1 and q2 in the domain
    if ismember(q1, sdn) && ismember(q2, sdn)  
        HL_new{q1, q2} = HL_temp{idx};
        SL_new{q1, q2} = SL_temp{idx};
        QL_new{q1, q2} = QL_temp{idx};
        SpL_new{q1, q2} = SpL_temp{idx};
        IdL_new{q1, q2} = IdL_temp{idx};
    end
end

% generate cell of ledgers TL
for q1 = 1:(n + 1)
    if ismember(q1, sdn)       % check that q1 in the domain
        sct = spin_combinations([n - 2 ,2], q1 - n/2 - 1);
            for j = 1:size(sct, 1)
                TL_new{q1}(j, 1) = sct(j, 1); 
                TL_new{q1}(j, 2) = sct(j, 2); 
                TL_new{q1}(j, 3) = size(HLss{sct(j, 1) + (n - 2)/2 + 1, sct(j, 1) + (n - 2)/2 + 1}, 1);
                TL_new{q1}(j, 4) = size(I4ss{sct(j, 2) + 2/2 + 1, sct(j, 2) + 2/2 + 1}, 1);
            end
    end
end

end