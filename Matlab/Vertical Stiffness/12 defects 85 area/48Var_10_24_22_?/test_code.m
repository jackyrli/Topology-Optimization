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
% 48 inputs
options = optimoptions('surrogateopt','PlotFcn','surrogateoptplot');
options.MaxFunctionEvaluations = 20000;
%% 41 inputs
lb = zeros(defects_num*4,1);
lb(defects_num + 1 : defects_num*3) = 0.02;
lb(defects_num*3 + 1: defects_num*4) = 0;

ub = zeros(defects_num*4,1);
ub(defects_num + 1 : defects_num*3) = 0.18;
ub(defects_num*3 + 1: defects_num*4) = pi/4;

for i = 1:defects_num
    lb(i) = 1;
    ub(i) = 26 - i;
end
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
intcon = 1:defects_num;
[x,fval,exitflag,output] = surrogateopt(objconstr,lb,ub,intcon,A,b,Aeq,beq,options)

diary off
%%
function solid_area_percentage = find_solid_area(v1)
defects_num = 12;
void_area = pi*v1(defects_num+1 : defects_num*2)*v1(defects_num*2 + 1 : defects_num*3)';
total_area = 2*2;
solid_area_percentage = (total_area - void_area)/total_area;
end

function [c,ceq] = constr(v1)
c = find_solid_area(v1) - 0.85;
ceq = [];
disp(['solid area percentage: ',num2str(find_solid_area(v1))])
end
