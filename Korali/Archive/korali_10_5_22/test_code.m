clear all;
clc;
diary logfile1
global nex ney LX LY defects_num
nex = 300;
ney = 300;
LX = 2;
LY = 2;
defects_num = 12;


optimize_struct = VM_structure(nex,ney,LX,LY,defects_num);

fun = @(v) -1*optimize_struct.compute(v);

%% surrogate
%rng default % For reproducibility
% nonlinear constraint and shift, no first 25, 
options = optimoptions('surrogateopt','PlotFcn','surrogateoptplot');
options.MaxFunctionEvaluations = 30000;
%% 25 + 12*3 = 61 inputs
lb = zeros(defects_num*3+25,1);
lb(26:2*defects_num+25) = 0.02;
lb(2*defects_num+26) = -pi/2;

ub = ones(defects_num*3+25,1);
ub(26:2*defects_num+25) = 0.18;
ub(2*defects_num+26) = pi/2;
%% linear constr
A = [];
b = [];

% Aeq = [ones(1,25), zeros(1,defects_num*3)];
% beq = 12;
Aeq = [];
beq = [];
% Linear constr attempt

%% nonlinear constr

%pack nonlinear constraints
objconstr = packfcn(fun,@constr);
% intcon=1:25;
intcon = [];
[x,fval,exitflag,output] = surrogateopt(objconstr,lb,ub,intcon,A,b,Aeq,beq,options)

diary off
%%
find_solid_area(a(1:100))
%%
function solid_area_percentage = find_solid_area(v1)
defects_num = 12;
void_area = pi*v1(26:25+defects_num)*v1(26+defects_num:25+2*defects_num)';
total_area = 2*2;
solid_area_percentage = (total_area - void_area)/total_area;
end

function num_defects_true = find_num_defects_from_input(v1)
    temp = round(v1(1:25));
    num_defects_true = ones(1,25)*temp';
end

% function [c,ceq] = constr(v1)
% c = [find_solid_area(v1) - 0.85
%     -1*find_solid_area(v1) + 0.69
%     find_num_defects_from_input(v1) - 12.1
%     -find_num_defects_from_input(v1) + 11.9];
% ceq = [];
% disp(['solid area percentage: ',num2str(find_solid_area(v1))])
% end

function [c,ceq] = constr(v1)
c = [find_solid_area(v1) - 0.85
    -1*find_solid_area(v1) + 0.69
    ];
ceq = [];
disp(['solid area percentage: ',num2str(find_solid_area(v1))])
end
