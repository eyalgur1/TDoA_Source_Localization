function realization = create_realization(array, sigma)

sensors = array.positions;
N = array.N;
s_real = array.s_real;


% set parameters for realization r and noise sigma (across all methods)
noises = sigma*randn(1,N);
dist = sqrt(sum((s_real-sensors).^2)) + noises;  % initial noised distances ||s-pi|| + eps_i
[~, ind_ref] = min(dist);  % index of reference sensor


% rearrange the network such that the first index corresponds to the reference sensor
if ind_ref == N
    sensors = [sensors(:,N) sensors(:,1:N-1)];
    dist = [dist(N) dist(1:N-1)];
    noises = [noises(N) noises(1:N-1)];
elseif ind_ref > 1
    sensors = [sensors(:, ind_ref) sensors(:, 1:ind_ref-1) sensors(:, ind_ref+1:N)];
    dist = [dist(ind_ref) dist(1:ind_ref-1) dist(ind_ref+1:N)];
    noises = [noises(ind_ref) noises(1:ind_ref-1) noises(ind_ref+1:N)];
end 


% save to structure 
realization.sigma = sigma;
realization.sensors = sensors;
realization.noises = noises;
realization.noised_distances = dist;
realization.d = dist - dist(1);  % vector of measurements ||s-pi|| - ||s-p1|| + eps_i - eps_1 (the first coordiante is 0)