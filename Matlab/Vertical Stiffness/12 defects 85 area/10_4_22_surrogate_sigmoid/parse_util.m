  clear all      
%%
load('training_result.mat')
vdata = parse_input(x, 12);
writematrix(vdata, 'result_surrogate_sigmoid.txt')

plot_fun('result_surrogate_sigmoid.txt')
% note: change in structure toplogy around 6000-6500 sim.txt
%%
function vdata = parse_input(v_input,num_defects)            
% start parse input from heavyside function
    tf_stats = v_input(1:25);
    tf_stats = round(tf_stats);
    % end heavyside parsing

    % start defects testing and parsing
    if sum(tf_stats) ~= num_defects
        vdata = 0; % indicating violation of constraint
    else
        vdata = zeros(1,100);
        vdata(1:25) = tf_stats;
        count_defects = 1;
        for i = 1:25
            if tf_stats(i) == 1 %defect in place    
                vdata(2*i+24 : 2*i+25) = [v_input(25+count_defects) v_input(25+num_defects+count_defects)];
                vdata(75+i) = v_input(25+2*num_defects+count_defects);
                count_defects = count_defects + 1;
            else % no defect in place
                vdata(2*i+24 : 2*i+25) = [0,0];
                vdata(75+i) = 0;
            end
        end
    end
end