  clear all      
%%
load('continuous_surrogate_opt_result.mat')
vdata = parse_input(x, 12);
writematrix(vdata, 'result_surrogate_cont.txt')

plot_fun('result_surrogate_cont.txt')
% note: change in structure toplogy around 6000-6500 sim.txt
%% parse input
parse_input([1,2,3,4,6,6,7,8,9,10,11,12,13],13)
function vdata = parse_input(v_input)
    %start parsing first 12
    index_array = 1:25; 
    % prev line: array we are deleting from, the remaining is the 0, deleted is 1
    void_indices = zeros(1,num_defects);
    v25 = zeros(1,25);
    v_def = v_input(1:num_defects); % frist #def elements from our input
    for i = 1: num_defects
        void_indices(i) = index_array(v_def(i));
        index_array(v_def(i)) = [];
    end
    v25(void_indices) = 1;
    % end parsing first 12
    vdata = zeros(1,100);
    count_defects = 1;
    for i = 1:25
        if v25(i) >= 0.5 %defect in place    
            vdata(2*i+24 : 2*i+25) = [v_input(25+count_defects) v_input(25+obj.num_defects+count_defects)];
            vdata(75+i) = v_input(25+2*obj.num_defects+count_defects);
            count_defects = count_defects + 1;
        else % no defect in place
            vdata(2*i+24 : 2*i+25) = [0,0];
            vdata(75+i) = 0;
        end
    end
end