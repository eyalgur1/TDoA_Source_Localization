function [funv, funML, s, times] = TNAM(array, realization, sAG, rAG, s, max_iter, Lip)
%% Initialization

N = array.N; I = array.I; e1 = array.e1; E1 = array.E1;
ML = array.ML; % Psi = array.Psi; % Psi is the smooth reformulation
sensors = realization.sensors;
array_no_ref = sensors(:, 2:N); ref = sensors(:, 1);  % the network without the reference
p = sum(array_no_ref, 2) + (N-1)*ref;
d = realization.d; d_sum = sum(d);

% set the Lipschitz constant function
if contains(Lip, 'd')
    Lip_const = @(A)4*(N-1); % distributed
elseif contains(Lip, 'c')
    Lip_const = @(A)2*norm(A);  % centralized
end

funv = [];  % for the smooth reformulation
funML = zeros(max_iter + 1,1); times = zeros(max_iter, 1);
funML(1) = ML(s,d,array_no_ref,ref);
tot_iter = 0;  % total iteration counter
out_iter = 0;


%%  Algorithm

while tot_iter < max_iter
    out_iter = out_iter + 1;
    time_a_start = tic;

    % update u
    u = (s - sensors).*E1;  % cahnge the sign of all non-reference sensors
    norm_u = sqrt(sum(u.^2)); nz_norm = (norm_u > 0);
    if sum(nz_norm) < N  % if s = pi
        u(:, nz_norm) = u(:, nz_norm)./norm_u(nz_norm);
        u(:, not(nz_norm)) = e1;
    else  % s is not any sensor
        u = u./norm_u;
    end

    % update A and L
    A = (N-1)*I + sum(u(:, 2:N), 2)*u(:, 1)'; A = 0.5*(A + A');
    L = Lip_const(A);

    % update z
    z1 = u(:, 1)*(sum(array_no_ref.*u(:, 2:N))); z2 = ref'*u(:, 1)*u(:, 2:N); z3 = d(2:N).*u(:, 2:N);
    z = p + sum(z1 + z2 - z3, 2);

    % CVX update (instead of FISTA)
    % cvx_begin
    % tot_iter = tot_iter + 1
    % variable s(n)
    % minimize quad_form(s, A) - (z')*s + d_sum*norm(s-ref)
    % cvx_end
    % funv(end+1) = Psi(s,u);

    % FISTA initilization
    sj = s; yj = s; tj = 1;
    iterAG = sAG + 2^floor(out_iter/rAG) - 1;
    times(tot_iter + 1) = toc(time_a_start);
    for j = 1:iterAG  % FISTA iterations

        if tot_iter >= max_iter
            break
        end
        tot_iter = tot_iter + 1;

        time_b_start = tic;
        sj_prev = sj;
        tj_prev = tj;

        aj = -(2/L)*A*yj + yj +(1/L)*z - ref; L_norma = L*norm(aj);
        sj = max(0, ((L_norma - d_sum)/L_norma))*aj + ref;
        tj = (1 + sqrt(1 + 4*tj^2))/2;
        yj = sj + ((tj_prev - 1)/tj)*(sj - sj_prev);

        times(tot_iter) = toc(time_b_start) + times(tot_iter);
        % funv(tot_iter+1) = Psi(sj,u,p,d,array_no_ref,ref);  % for the smooth reformulation
        funML(tot_iter + 1) = ML(sj, d, array_no_ref, ref);
    end
    s = sj;
end

% %% Cumulative Times Calculation
% cum_times = cumsum(times);
% time_x_axis = 0:0.00002:0.04;
% time_iter=zeros(length(time_x_axis) ,1);
% 
% for l = 1:length(time_x_axis)
%     ll = 1;
%     while cum_times(ll) < time_x_axis(l)
%         ll = ll+1;
%         if ll >= length(cum_times)
%             break
%         end
%     end
%     time_iter(l) = ll;
% end
% funML_cum = funML(time_iter);