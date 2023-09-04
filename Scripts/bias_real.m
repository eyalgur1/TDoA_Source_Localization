num_nets = 200;
R = 100;
Sigma = logspace(-4, -1, 10)';
sF = [1,2,3,4,5];
methods_to_run = {'TNAM', 'FP', 'SOLVIT'};


for j = 1:length(Sigma)
    sigma = Sigma(j);

    for k = 1:num_nets
        for m = 1:length(methods_to_run)
            method = methods_to_run{m};
            switch method

                case 'TNAM'
                    for ss = 1:length(sF)
                        ssF = sF(ss);
                        %s_out = output.(['net', num2str(k)]).(['sigma', num2str(j)]).(method).(['sF',num2str(ssF)]).s_out;
                        output.stat.(['sigma', num2str(j)]).(method).(['sF',num2str(ssF)]).avg_out =  [0;0];
                    end

                otherwise
                    %s_out = output.(['net', num2str(k)]).(['sigma', num2str(j)]).(method).s_out;
                    output.stat.(['sigma', num2str(j)]).(method).avg_out =  [0;0];
                    %output.(['net', num2str(k)]).(['sigma', num2str(j)]).(method).RMSE =  RMSE_fun(s_out);
                    %output.(['net', num2str(k)]).(['sigma', num2str(j)]).(method).bias =  bias_fun(s_out);
            end
        end
    end
end