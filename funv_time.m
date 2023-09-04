function [funv_cumulative, len_iter] = funv_time(intervals, times, funv)

% Inputs:
% intervals: line segment of a:b:c where a > 0
% times: array of run time of each iteration (output of an iterative algorithm)
% funv: array of function value in each iteration (output of an iterative algorithm)

cumulative_times = cumsum(times);  % array of the size max_iter with cumulative running times
len_cum = length(cumulative_times);
len_intervals = length(intervals);
time_iter = zeros(length(cumulative_times) ,1);  % pre-allocation for the iterations that correspond with cumulative times
funv_cumulative = zeros(len_intervals, 1); 

for i = 1:length(intervals) 
    j = find(cumulative_times < intervals(i), 1, 'last');  % for each time segment [0, intervals(i)], find the max iteration j for which cumulative_times(j) < intervals(i)
    if ~isscalar(j)  % can only happen in the first iterations: if j is not a scalar, then the algorithm's accumulate time up to this point is less than intervals(1)
        j = 1;
    end
    time_iter(i) = j;  % set the correpsoding iteration index
    if j >= len_cum  % stop if reached the last iteration (the algorithm run time is shorter than intervals(end))
        break
    end
end

time_iter = [1; time_iter];  % add the first iteration (time 0)
time_iter(time_iter == 0) = [];  % remove all indices for which the time intervals are greater than the run time (removes all last iterations if needed)
if length(time_iter) > len_intervals
    time_iter = time_iter(1:len_intervals);
end
len_iter = length(time_iter);
funv_cumulative(1:len_iter) = funv(time_iter);  % array of function values that correpsond to the time segments