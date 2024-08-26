% v1. 11/12/2023
% v3. 05/14/2024  (included theta angle)
% v4. 08/15/2024  (added spin projections domains and parallelization)
% This function builds right blocks in spin-spin sectors for the Schwinger model
% n is the number of spin-1/2 sites in the final right blocks
% the output is cells by (n + 1) x (n + 1) size (all spin sectors)
function [ HR_new, QR_new, SR_new, SpR_new, IdR_new, TR_new ] = build_Schwinger_blocks_right_ss(HRss, ...
                                  QRss, SRss, SpRss, IdRss, H2ss, Q2ss, S2ss, M2ss, ...
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
HR_new = cell(n + 1, n + 1);
QR_new = cell(n + 1, n + 1);
SR_new = cell(n + 1, n + 1);
SpR_new = cell(n + 1, n + 1);
IdR_new = cell(n + 1, n + 1);

% Pre-allocate the arrays to store intermediate results
HR_temp = cell((n + 1)^2, 1);
QR_temp = cell((n + 1)^2, 1);
SR_temp = cell((n + 1)^2, 1);
SpR_temp = cell((n + 1)^2, 1);
IdR_temp = cell((n + 1)^2, 1);

% Total number of iterations
total_iterations = (n + 1) * (n + 1);

% define auxiliarly operators 
Sm2Lss = cellfun(@ctranspose,Sp2Lss','UniformOutput',false);     
SmRss = cellfun(@ctranspose, SpRss','UniformOutput',false);
Q2Sqss = cellfun( @(x) x.^2 + (1/2 + theta/pi)*x, Q2ss, 'UniformOutput', false);                           
% run over all spin sectors ( q = s + n/2 + 1 = 1,...,n + 1 )
parfor idx = 1:total_iterations
% Calculate the corresponding q1 and q2 values
    [q1, q2] = ind2sub([n + 1, n + 1], idx);
     if ismember(q1, sdn) && ismember(q2, sdn)  
     % Compute the results and store them in the temporary arrays

        % build right blocks
        % build new HR       
        HR_temp{idx} = sparse(kron_ss(I4ss, HRss, q1, q2, [2, n - 2]) ...
         + x*(kron_ss(Sm2Lss, SpRss, q1, q2, [2, n - 2]) ...
              + kron_ss(Sp2Lss, SmRss, q1, q2, [2, n - 2])...
              + kron_ss(HXY2ss, IdRss, q1, q2, [2, n - 2])) ...
         + (mu/2)*kron_ss(M2ss, IdRss, q1, q2, [2, n - 2]) ...
         + y*(kron_ss(H2ss, IdRss, q1, q2, [2, n - 2]) ... 
              + kron_ss(Q2ss, SRss, q1, q2, [2, n - 2]) ...
              + (n - 2)*kron_ss(Q2Sqss, IdRss, q1, q2, [2, n - 2]) ...
              + (theta^2/(2*pi)^2)*kron_ss(I4ss, IdRss, q1, q2, [2, n - 2]))); 
                     
        % build new QR             
        QR_temp{idx} = sparse(kron_ss(I4ss, QRss, q1, q2, [2, n - 2]) ... 
                            + kron_ss(Q2ss, IdRss, q1, q2, [2, n - 2])); 
                     
        % build new SR
        SR_temp{idx} = sparse(kron_ss(I4ss, SRss, q1, q2, [2, n - 2])... 
                         + kron_ss(S2ss, IdRss, q1, q2, [2, n - 2]) ...
                         + 2*(n - 2)*kron_ss(Q2ss, IdRss, q1, q2, [2, n - 2]));
        % build new SpR and IdR
        SpR_temp{idx} = sparse(kron_ss(Sp2Rss, IdRss, q1, q2, [2, n - 2]));
        IdR_temp{idx} = sparse(kron_ss(I4ss, IdRss, q1, q2, [2, n - 2]));
      end
end

% Reconstruct the 2D cell arrays from the linearized temporary arrays
for idx = 1:total_iterations
    [q1, q2] = ind2sub([n + 1, n + 1], idx);
    % check that q1 and q2 in the domain
    if ismember(q1, sdn) && ismember(q2, sdn)  
        HR_new{q1, q2} = HR_temp{idx};
        SR_new{q1, q2} = SR_temp{idx};
        QR_new{q1, q2} = QR_temp{idx};
        SpR_new{q1, q2} = SpR_temp{idx};
        IdR_new{q1, q2} = IdR_temp{idx};
    end
end

% generate cell of ledgers TR
for q1 = 1:(n + 1)
    if ismember(q1, sdn)   % check that q1 in the domain
        sct = spin_combinations([2, n - 2], q1 - n/2 - 1);
        for j = 1:size(sct, 1)
           TR_new{q1}(j, 1) = sct(j, 1);
           TR_new{q1}(j, 2) = sct(j, 2);
           TR_new{q1}(j, 3) = size(I4ss{sct(j, 1) + 2/2 + 1, sct(j, 1) + 2/2 + 1}, 1);
           TR_new{q1}(j, 4) = size(HRss{sct(j, 2) + (n - 2)/2 + 1, sct(j, 2) + (n - 2)/2 + 1}, 1);
        end
    end
end

end

