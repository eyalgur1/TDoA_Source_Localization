%% Set Hyper-Parameters

clear; close all; cvx_solver SeDuMi     % close all and set CVX solver for SDP                         
N = 4;                                  % number of sensors                             % set the seed according to past tests
rng(312)
n = 2;                                  % dimension
Lip = 'c';                              % centralized ('c') or distribured ('d') Lipschitz constant for TNAM
num_nets = 200;                         % number of arrays to run
R = 100;                                % number of realizations (Monte Carlo) for each net
Sigma = logspace(-4, -1, 10)';          % values of noise
sF = [1,2,3,4,5]; rF = 1000;            % parameters of FISTA inner iterations
methods_to_run = ...                    % set of all required methods from {'TNAM', 'FP', 'SOLVIT', 'SLS', 'SCWLS', 'SDP'}
    {'SDP'};           % the correct order is: TNAM -> FP -> SOLVIT -> others (the others can also be without the former three)
near_field = 1;                         % if near field, then s_real is in the vicinity of sensors; otherwise it is outside the convex hull
approx_circular = 0;                    % approximately circular array or not
if near_field                           % set the maximum number of iterations according to past tests
    max_iter = 2000;
else
    max_iter = 5000;                        
end
Q_mat = 'Id';                           % covariance matrix fo SLS and SCWLS methods; 'Id' is the identity and 'sig' depends on sigma
eta = 10^(-7);                          % regularizer for SDP
intervals = 0:0.00005:0.04;             % time segments intervals for function value plot over time


% set output structure and names to plot legend (this section is set
% automotically - do not set manullay)
output = struct; output.stat.near_field = near_field; output.stat.approx_circular = approx_circular; output.stat.cov_mat_SLS = Q_mat;
TNAM_exists = sum(contains(methods_to_run, 'TNAM'));
num_of_methods = TNAM_exists*(length(methods_to_run) - 1 + length(sF)) + not(TNAM_exists)*length(methods_to_run);  % the number of total methods (inclusing TNAM with different sF)
names_to_legend = cell(num_of_methods, 1);  % methods names to graph legends
SLS_indicator = (sum(contains(methods_to_run, 'SLS')) == 1);
SDP_indicator = (sum(contains(methods_to_run, 'SDP')) == 1);
SCWLS_indicator = (sum(contains(methods_to_run, 'SCWLS')) == 1);
names_to_legend_time = cell(num_of_methods - SLS_indicator - SDP_indicator - SCWLS_indicator, 1);
tableVars = ["Method", "Avg Total Run Time", "Avg Mean Time Per Iter"];  % column names of running time tables
time_plot_length = length(intervals);


%% Run Algorithms

for k = 1:num_nets  % create array positions

    positions = -0.5 + rand(n, N);  % create random positions in the box [-0.5, 0.5]
    positions = approx_circular*(positions./sqrt(sum(positions.^2))) + (1 - approx_circular)*positions;  % approximately circular (by projection onto sphere) or random in the box
    s_real = near_field*(-0.5 + rand(n, 1)) + (1 - near_field)*(0.5 + rand(n, 1));  % real source position outside the convex hull (up to distance 1) or within the vicinity of sensors
    output.(['net', num2str(k)]).sensors = positions; output.(['net', num2str(k)]).s_real = s_real;
    array = create_array(N, n, positions, s_real);


    for j = 1:length(Sigma)  % for each value of sigma, set starting point for all methods (any set of Monte Carlo trials have the same s0)
        sigma = Sigma(j);
        s0 = near_field*(-0.5 + rand(n, 1)) + (1 - near_field)*(0.5 + rand(n, 1));

        output.(['net', num2str(k)]).(['sigma', num2str(j)]).sigma = sigma;
        output.(['net', num2str(k)]).(['sigma', num2str(j)]).s0 = s0;
        ML_val = zeros(length(sF), max_iter + 1, length(methods_to_run)); ML_cum = zeros(length(sF), time_plot_length, length(methods_to_run));  % pre-allocation for function values (along the rows), where each method is in dimension 3 and TNAM has several rows for different values of sF
        method_len_iter = zeros(length(sF), R, length(methods_to_run));
        s_out = zeros(n, R, length(sF), length(methods_to_run)); times_out = zeros(2, length(sF), length(methods_to_run));  % pre-allocation for the output positions of the source, and mean (first row) and total (second row) execution time

        for r = 1:R  % create a realization of the measurements for each value of sigma, for all methods
            realization = create_realization(array, sigma);
            fprintf(['net=',num2str(k),'/',num2str(num_nets),' || sigma=',num2str(sigma),' (', num2str(j),'/',num2str(length(Sigma)),') || realization=',num2str(r),'/',num2str(R),'\n'])

            % generate a different starting point for any realization according to SLS solution
            % h = 0.5*(realization.d(2:N) - sum(realization.sensors(:, 2:end).^2) - norm(realization.sensors(:, 1)^2))';
            % G = -[realization.sensors(:, 2:end) - realization.sensors(:, 1); d(2:N)]';
            % s0 = G\h; s0 = s0(1:n);

            for m = 1:length(methods_to_run)  % run all methods for each realization
                method = methods_to_run{m};

                switch method

                    case 'TNAM'
                        for ss = 1:length(sF)
                            ssF = sF(ss);
                            [~, ML_fun, s, times] = TNAM(array, realization, ssF, rF, s0, max_iter, Lip);
                            ML_val(ss, :, m) = ML_val(ss, :, m) + (1/R)*ML_fun';
                            [funv_cumulative, len_iter] = funv_time(intervals, times, ML_fun);
                            ML_cum(ss, :, m) = ML_cum(ss, :, m) + (1/R)*funv_cumulative'; method_len_iter(ss, r, m) = len_iter;
                            s_out(1:n, r, ss, m) = s;
                            times_out(1, ss, m) = times_out(1, ss, m) + (1/R)*sum(times); times_out(2, ss, m) = times_out(2, ss, m) + (1/R)*mean(times);
                        end

                    case {'SLS','SDP','SCWLS'}  % non-iterative methods
                        if contains(method, 'SLS')
                            [ML_fun, s, times] = SLS(array, realization, Q_mat);
                        elseif contains(method, 'SDP')
                            [ML_fun, s, times] = SDP(array, realization, eta);
                        elseif contains(method, 'SCWLS')
                            [ML_fun, s, times] = SCWLS(array, realization, Q_mat);
                        end
                        ML_val(1, :, m) = ML_val(1, :, m) + (1/R)*ML_fun;
                        s_out(1:n, r, 1, m) = s;
                        times_out(1, 1, m) = times_out(1, 1, m) + (1/R)*times; times_out(2, 1, m) = 0;

                    otherwise  % iterative methods (except for T-NAM)
                        if contains(method, 'FP')
                            [ML_fun, s, times] = FP(array, realization, s0, max_iter);
                        elseif contains(method, 'SOLVIT')
                            [ML_fun, s, times] = SOLVIT(array, realization, s0, max_iter);
                        end
                        ML_val(1, :, m) = ML_val(1, :, m) + (1/R)*ML_fun';
                        [funv_cumulative, len_iter] = funv_time(intervals, times, ML_fun);
                        ML_cum(1, :, m) = ML_cum(1, :, m) + (1/R)*funv_cumulative'; method_len_iter(1, r, m) = len_iter;
                        s_out(1:n, r, 1, m) = s;
                        times_out(1, 1, m) = times_out(1, 1, m) + (1/R)*sum(times); times_out(2, 1, m) = times_out(2, 1, m) + (1/R)*mean(times);
                end
            end
        end


        % save arrays of average function values, location outputs and running times for each value of sigma
        for m = 1:length(methods_to_run)
            method = methods_to_run{m};
            switch method

                case 'TNAM'
                    for ss = 1:length(sF)
                        ssF = sF(ss);
                        output.(['net', num2str(k)]).(['sigma', num2str(j)]).(method).(['sF',num2str(ssF)]).avg_ML_val =  ML_val(ss, :, m);
                        output.(['net', num2str(k)]).(['sigma', num2str(j)]).(method).(['sF',num2str(ssF)]).avg_ML_cum =  ML_cum(ss, :, m);
                        output.(['net', num2str(k)]).(['sigma', num2str(j)]).(method).(['sF',num2str(ssF)]).s_out =  s_out(:, :, ss, m);
                        output.(['net', num2str(k)]).(['sigma', num2str(j)]).(method).(['sF',num2str(ssF)]).avg_mean_time_iter =  times_out(2, ss, m);
                        output.(['net', num2str(k)]).(['sigma', num2str(j)]).(method).(['sF',num2str(ssF)]).avg_tot_time =  times_out(1, ss, m);
                    end

                otherwise
                    output.(['net', num2str(k)]).(['sigma', num2str(j)]).(method).avg_ML_val =  ML_val(1, :, m);
                    if not(contains(method, 'SLS')) && not(contains(method, 'SDP')) && not(contains(method, 'SCWLS'))
                        output.(['net', num2str(k)]).(['sigma', num2str(j)]).(method).avg_ML_cum =  ML_cum(1, :, m);
                    end
                    output.(['net', num2str(k)]).(['sigma', num2str(j)]).(method).s_out =  s_out(:, :, 1, m);
                    output.(['net', num2str(k)]).(['sigma', num2str(j)]).(method).avg_mean_time_iter =  times_out(2, 1, m);
                    output.(['net', num2str(k)]).(['sigma', num2str(j)]).(method).avg_tot_time =  times_out(1, 1, m);
            end
        end
    end
end


%% Plot avearge ML function values (vs. iterations and run time) and run time tables

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultLineLineWidth',1)
set(groot,'defaultAxesFontSize',14)


value_vs_sigma = zeros(length(sF), length(Sigma), length(methods_to_run));
for j = 1:length(Sigma)
    current_method = 0; current_method_time = 0;
    figure(j); figure(j + 100);  % one plot for function value along iterations and another along accumulated run time
    sigma = Sigma(j);
    ML_to_plot = zeros(length(sF), max_iter + 1, length(methods_to_run)); ML_cum_to_plot = zeros(length(sF), time_plot_length, length(methods_to_run));
    table_columns = zeros(num_of_methods, 2);


    % plot function values and create columns of run times table to each value of sigma
    for m = 1:length(methods_to_run)
        method = methods_to_run{m};
        switch method

            case 'TNAM'
                for ss = 1:length(sF)
                    current_method = current_method + 1; current_method_time = current_method_time + 1;
                    ssF = sF(ss);
                    for k = 1:num_nets
                        ML_to_plot(ss, :, m) = ML_to_plot(ss, :, m) + (1/num_nets)*output.(['net', num2str(k)]).(['sigma', num2str(j)]).(method).(['sF',num2str(ssF)]).avg_ML_val;
                        ML_cum_to_plot(ss, :, m) = ML_cum_to_plot(ss, :, m) + (1/num_nets)*output.(['net', num2str(k)]).(['sigma', num2str(j)]).(method).(['sF',num2str(ssF)]).avg_ML_cum;
                        table_columns(current_method, :) = table_columns(current_method, :) + (1/num_nets)*[output.(['net', num2str(k)]).(['sigma', num2str(j)]).(method).(['sF',num2str(ssF)]).avg_tot_time;...
                            output.(['net', num2str(k)]).(['sigma', num2str(j)]).(method).(['sF',num2str(ssF)]).avg_mean_time_iter]';
                    end
                    value_vs_sigma(ss, j, m) = ML_to_plot(ss, end, m);
                    names_to_legend{current_method} = ['$\mathrm{TNAM}, s=',num2str(ssF),'$']; names_to_legend_time{current_method_time} = names_to_legend{current_method};
                    figure(j); semilogy(ML_to_plot(ss, :, m)); hold on  % function value vs. iteration plots
                    figure(j + 100); semilogy(ML_cum_to_plot(ss, 1:min(method_len_iter(ss, :, m)), m)); hold on  % function value vs. run time plots
                end

            otherwise
                current_method = current_method + 1;
                if not(contains(method, 'SLS')) && not(contains(method, 'SDP')) && not(contains(method, 'SCWLS'))
                    current_method_time = current_method_time + 1;
                    names_to_legend_time{current_method_time} = method;
                end
                for k = 1:num_nets
                    ML_to_plot(1, :, m) = ML_to_plot(1, :, m) + (1/num_nets)*output.(['net', num2str(k)]).(['sigma', num2str(j)]).(method).avg_ML_val;
                    if not(contains(method, 'SLS')) && not(contains(method, 'SDP')) && not(contains(method, 'SCWLS'))
                        ML_cum_to_plot(1, :, m) = ML_cum_to_plot(1, :, m) + (1/num_nets)*output.(['net', num2str(k)]).(['sigma', num2str(j)]).(method).avg_ML_cum;
                    end
                    table_columns(current_method, :) = table_columns(current_method, :) + (1/num_nets)*[output.(['net', num2str(k)]).(['sigma', num2str(j)]).(method).avg_tot_time;...
                        output.(['net', num2str(k)]).(['sigma', num2str(j)]).(method).avg_mean_time_iter]';
                end
                value_vs_sigma(1, j, m) = ML_to_plot(1, end, m);
                names_to_legend{current_method} = method;
                figure(j); semilogy(ML_to_plot(1, :, m)); hold on  % function value vs. iteration plots
                if not(contains(method, 'SLS')) && not(contains(method, 'SDP')) && not(contains(method, 'SCWLS'))  % function value vs. run time plots
                    figure(j + 100); semilogy(ML_cum_to_plot(1, 1:min(method_len_iter(1, :, m)), m)); hold on
                end
        end
    end

    % create title and legend to each function value plot
    figure(j)  % function value vs. iteration plots
    legend(names_to_legend); set(gca, 'YScale', 'log'); grid on; xlim([1 max_iter + 1])
    title(['$\sigma=',num2str(sigma),'$']); ylabel('$\mathrm{Average\ ML\ function\ value}$'); xlabel('$\mathrm{Iterations}$')
    output.stat.(['sigma', num2str(j)]).ML_plot = figure(j);
    hold off

    figure(j + 100)  % function value vs. run time plots
    legend(names_to_legend_time); set(gca, 'YScale', 'log'); grid on; xlim([1 length(intervals)]); xticks(intervals(1):length(intervals)/10:length(intervals)); xticklabels([intervals(1:length(intervals)/10:end)])
    title(['$\sigma=',num2str(sigma),'$'])
    ylabel('$\mathrm{Average\ ML\ function\ value}$'); xlabel(['$\mathrm{Accumulated\ run\ time\ (in\ seconds)\ until\ }', num2str(max_iter), '\mathrm{\ iterations\ are\ reached}$'])
    output.stat.(['sigma', num2str(j)]).ML_cum_plot = figure(j + 100);
    hold off


    % running times table to each value of sigma
    run_time_tab = table(string(names_to_legend), table_columns(:, 1), table_columns(:, 2), 'VariableNames', tableVars);
    run_time_tab = table(run_time_tab, 'VariableNames', {['Run Time, sigma=',num2str(sigma)]});  % table title
    output.stat.(['sigma', num2str(j)]).run_time_table = run_time_tab;
end

%%
% plot final function values
figure(200); hold on
for m = 1:length(methods_to_run)
    method = methods_to_run{m};
    switch method
        case 'TNAM'
            for ss = 1:length(sF)
                plot(log10(Sigma), value_vs_sigma(ss,:,m))
            end

        otherwise
            plot(log10(Sigma), value_vs_sigma(1,:,m))
    end
end
legend(names_to_legend); set(gca, 'YScale', 'log'); grid on;
ylabel('$\mathrm{Average\ ML\ function\ value}$'); xlim([log10(min(Sigma)) log10(max(Sigma))]); xlabel('$\log_{10}(\sigma)$', 'Interpreter', 'latex');
title('$\mathrm{Output\ Function\ Value}$');
output.stat.final_funv_plot = figure(200);


%% Calculate and save CRLB, RMSE and norm of bias

% function handles of Jacobian, covarince matrix, RMSE and norm of bias
J_fun = @(s, sensors)((s - sensors(:, 2:N))./sqrt(sum((s - sensors(:, 2:N)).^2)) - (s - sensors(:, 1))/norm(s - sensors(:, 1)))';  % Jacobiam matrix
cov_fun = @(sigma)(sigma^2)*ones(N-1) + (sigma^2)*eye(N-1);  % covariance matrix
RMSE_fun = @(s_out)sqrt(trace((1/R)*sum(reshape(kron((s_out - mean(s_out ,2)), ones(1, n))...
    .*reshape(s_out - mean(s_out ,2), 1 ,[]), n, n, []), 3)));  % RMSE function
bias_fun = @(s_out)norm((1/R)*sum(s_out - mean(s_out ,2), 2));  % norm of bias function


% calculate and save CRLB, RMSE and bias for every net, for each value of sigma
for j = 1:length(Sigma)
    sigma = Sigma(j);
    cov = cov_fun(sigma);  % covariance matrix (constant over the nets for each value of sigma)

    for k = 1:num_nets
        % CRLB for every net
        J = J_fun(output.(['net', num2str(k)]).s_real, output.(['net', num2str(k)]).sensors);
        FIM =  J'*(cov\J);
        output.(['net', num2str(k)]).(['sigma', num2str(j)]).CRLB = sqrt(trace(inv(FIM)));

        % RMSE and bias for every net and method
        for m = 1:length(methods_to_run)
            method = methods_to_run{m};
            switch method

                case 'TNAM'
                    for ss = 1:length(sF)
                        ssF = sF(ss);
                        s_out = output.(['net', num2str(k)]).(['sigma', num2str(j)]).(method).(['sF',num2str(ssF)]).s_out;
                        output.(['net', num2str(k)]).(['sigma', num2str(j)]).(method).(['sF',num2str(ssF)]).RMSE =  RMSE_fun(s_out);
                        output.(['net', num2str(k)]).(['sigma', num2str(j)]).(method).(['sF',num2str(ssF)]).bias =  bias_fun(s_out);
                    end

                otherwise
                    s_out = output.(['net', num2str(k)]).(['sigma', num2str(j)]).(method).s_out;
                    output.(['net', num2str(k)]).(['sigma', num2str(j)]).(method).RMSE =  RMSE_fun(s_out);
                    output.(['net', num2str(k)]).(['sigma', num2str(j)]).(method).bias =  bias_fun(s_out);
            end
        end
    end
end


% calculate and save average CRLB, RMSE and bias over the nets, for each value of sigma
for j = 1:length(Sigma)
    sigma = Sigma(j);
    avg_CRLB = 0;

    % calculate average CRLB
    for k = 1:num_nets
        avg_CRLB = avg_CRLB + (1/num_nets)*output.(['net', num2str(k)]).(['sigma', num2str(j)]).CRLB;
    end
    output.stat.(['sigma', num2str(j)]).sigma = sigma;
    output.stat.(['sigma', num2str(j)]).avg_CRLB = avg_CRLB;

    % calculate average RMSE and bias
    for m = 1:length(methods_to_run)
        method = methods_to_run{m};
        avg_RMSE = 0; avg_bias = 0;
        switch method

            case 'TNAM'
                for ss = 1:length(sF)
                    ssF = sF(ss);
                    for k = 1:num_nets
                        avg_RMSE = avg_RMSE + (1/num_nets)*output.(['net', num2str(k)]).(['sigma', num2str(j)]).(method).(['sF',num2str(ssF)]).RMSE;
                        avg_bias = avg_bias + (1/num_nets)*output.(['net', num2str(k)]).(['sigma', num2str(j)]).(method).(['sF',num2str(ssF)]).bias;
                    end
                    output.stat.(['sigma', num2str(j)]).(method).(['sF',num2str(ssF)]).avg_RMSE = avg_RMSE;
                    output.stat.(['sigma', num2str(j)]).(method).(['sF',num2str(ssF)]).avg_bias = avg_bias;
                end

            otherwise
                for k = 1:num_nets
                    avg_RMSE = avg_RMSE + (1/num_nets)*output.(['net', num2str(k)]).(['sigma', num2str(j)]).(method).RMSE;
                    avg_bias = avg_bias + (1/num_nets)*output.(['net', num2str(k)]).(['sigma', num2str(j)]).(method).bias;
                end
                s_out = output.(['net', num2str(k)]).(['sigma', num2str(j)]).(method).s_out;
                output.stat.(['sigma', num2str(j)]).(method).avg_RMSE = avg_RMSE;
                output.stat.(['sigma', num2str(j)]).(method).avg_bias = avg_bias;
        end
    end
end


%% Plot avearge norm of bias, RMSE and CRLB

CRLB_line = zeros(length(Sigma), 1);
RMSE_line = zeros(length(Sigma), num_of_methods); bias_line = RMSE_line;

for j = 1:length(Sigma)
    sigma = Sigma(j);
    current_method = 0;
    CRLB_line(j) = output.stat.(['sigma', num2str(j)]).avg_CRLB;

    for m = 1:length(methods_to_run)
        method = methods_to_run{m};


        switch method

            case 'TNAM'
                for ss = 1:length(sF)
                    current_method = current_method + 1;
                    ssF = sF(ss);
                    RMSE_line(j, current_method) = output.stat.(['sigma', num2str(j)]).(method).(['sF',num2str(ssF)]).avg_RMSE;
                    bias_line(j, current_method) = output.stat.(['sigma', num2str(j)]).(method).(['sF',num2str(ssF)]).avg_bias;
                    names_to_legend{current_method} = ['$\mathrm{TNAM}, s=',num2str(ssF),'$'];
                end

            otherwise
                current_method = current_method + 1;
                RMSE_line(j, current_method) = output.stat.(['sigma', num2str(j)]).(method).avg_RMSE;
                bias_line(j, current_method) = output.stat.(['sigma', num2str(j)]).(method).avg_bias;
                names_to_legend{current_method} = method;
        end
    end
end

figure(1000)
plot(bias_line, 'LineWidth', 1)
legend(names_to_legend, 'Location', 'northwest'); grid on; xlabel('$\sigma$', 'Interpreter', 'latex'); xticklabels(Sigma);
title('$\mathrm{Norm\ of\ Bias}$')
output.stat.bias_plot = figure(1000);

figure(1001)
plot(log10(Sigma),RMSE_line); hold on
plot(log10(Sigma),CRLB_line, 'black')
legend([names_to_legend; 'CRLB'], 'Location', 'northwest'); grid on; xlabel('$\log_{10}(\sigma)$', 'Interpreter', 'latex');
title('$\mathrm{RMSE\ (linear\ scale)}$');
output.stat.RMSE_lin_plot = figure(1001);

figure(1002)
plot(log10(Sigma),log10(RMSE_line)); hold on
plot(log10(Sigma),log10(CRLB_line), 'black')
legend([names_to_legend; 'CRLB'], 'Location', 'northwest'); grid on; xlabel('$\log_{10}(\sigma)$', 'Interpreter', 'latex');
title('$\mathrm{RMSE\ (}\log_{10}\mathrm{\ scale)}$')
output.stat.RMSE_log_plot = figure(1002);


%% Save the Output Structure
if ~near_field
    geo = "far";
elseif approx_circular
    geo = "circ";
else 
    geo = "near";
end
save("output_TDOA_"+string(N)+"sen_"+geo+"_"+string(char(datetime('now','Format','yyyy_MM_dd_hh_mm_ss')))+".mat", 'output')