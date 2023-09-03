function [funML, s, times] = FP(array, realization, s, max_iter)
%% Initialization

N = array.N; ML = array.ML;
sensors = realization.sensors; d = realization.d;
array_no_ref = sensors(:, 2:N); ref = sensors(:, 1);  % the network without the reference

% vectors and constancs for the algorithm
b = sensors;
b(:, 1) = (N - 1)*b(:, 1);
w = @(s)[sum(sqrt(sum((s - array_no_ref).^2)) - d(2:N)), norm(s - ref) + d(2:N)];
g = @(s)(s - sensors)./sqrt(sum((s - sensors).^2));
coeff = 1/(2*(N-1));

funML = zeros(max_iter + 1, 1); times = zeros(max_iter, 1);
funML(1) = ML(s, d, array_no_ref, ref);
iter = 0;


%%  Algorithm
while iter < max_iter
    iter = iter + 1; tic
    s = coeff*(sum(b + w(s).*g(s), 2));  % update position s
    times(iter) = toc;
    funML(iter + 1) = ML(s, d, array_no_ref, ref);
end

% %%
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