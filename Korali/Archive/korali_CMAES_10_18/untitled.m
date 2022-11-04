rng default % For reproducibility
objconstr = @(x)(4*x(:,1).^2 - 2.1*x(:,1).^4 + x(:,1).^6/3 ...
    + x(:,1).*x(:,2) - 4*x(:,2).^2 + 4*x(:,2).^4);
lb = [-2.1,-2.1];
ub = -lb;
options = optimoptions('surrogateopt','MaxFunctionEvaluations',20);
[x,fval,exitflag,output,trials] = surrogateopt(objconstr,lb,ub,options);


%%
angle = pi*(rand(30,25)-0.5);
v = ones(30,50);% size 30 * 50
v(:,1:2:50) = 0.17;
v(:,2:2:50) = 0.1247;
% first 10 shape
load('first.mat');
for i = 1:2:50
    if first(floor((i+1)/2))==0
        v(1:10,i:i+1) = 0;
    end
end
% second

%third
