function [funML, s, times] = SOLVIT(array, realization, s, max_iter)
%% Initialization

N = array.N; n = array.n; I = array.I; ML = array.ML;
sensors = realization.sensors; d = realization.d;
array_no_ref = sensors(:, 2:N); ref = sensors(:, 1);  % the network without the reference

funML = zeros(max_iter+1,1); times = zeros(max_iter, 1);
funML(1) = ML(s,d,array_no_ref,ref);

iter = 0;


%%  Algorithm

while iter < max_iter
    iter = iter + 1; tic

    % constants for the matrix Q and vector b
    Q1 = s - ref;
    Q2 = s - array_no_ref;
    normQ2 = sqrt(sum(Q2.^2));
    Q3 = Q2./normQ2;

    % the Qi matrices as a row block matrix
    Q = (Q1/norm(Q1)).*kron(reshape(Q3, 1, n*(N-1)), ones(n,1));
    QT = reshape(permute(reshape(Q, n, n, []),[1,3,2]), n*(N-1), [])';  % taking transpose over the blocks

    % update the vector b and the matrix M
    alpha = d(2:N)./norm(s - ref);
    b = sum(array_no_ref, 2) + (N-1)*ref + sum(Q3.*d(2:N), 2) + ...
        sum(alpha.*ref, 2) - sum(Q.*reshape(array_no_ref, 1, n*(N-1)), 2) - ...
        sum(QT.*repmat(ref', [1 N-1]), 2);
    M = kron(2 + alpha, I) - Q - QT; % the Mi matrices as a row block matrix

    % update the next iteration
    s = sum(reshape(M, n, n, N-1), 3)\sum(b, 2);
    times(iter) = toc;
    funML(iter + 1) = ML(s,d,array_no_ref,ref);
end