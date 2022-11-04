clear all; clc; close
vdata = load('sim_results/sim6.txt');
[c,ceq] = constr(vdata(1:100))

%% plot cost function throughout iterations
result = zeros(3199,1);
for i = 0:3198
    temp = readmatrix(['sim_results/sim',num2str(i),'.txt']);
    result(i + 1) = temp(101);
end
figure
scatter(0:3198, result)
xlabel('iterations')
ylabel('cost function value')
title('Korali CMAES 10/18/22')

%%
function solid_area_percentage = find_solid_area(v1)
%input: v1 is 50*1 array, represents the major axis and minor axis of
%ellipse
% this function is for l1 l2 alternating
%output: solid_area percentage of the 2d structure
void_area = 0;
for i = 1:2:50
    void_area = void_area + pi*v1(i)*v1(i+1);
end
total_area = 2*2;
solid_area_percentage = (total_area - void_area)/total_area;
end

function location = generate_location(v1)
v2 = v1(1:50);
v1 = zeros(1,25);
for i = 1:2:50
    if v2(i)==0 || v2(i+1)==0
       v1((i+1)/2)=0;
    else
       v1((i+1)/2)=1;
    end
end
location = v1;
end

function [c,ceq] = constr(v1)
c = [find_solid_area(v1) - 0.85
    -1*find_solid_area(v1) + 0.69
    ];
ceq = [generate_location(v1)*ones(25,1)-12];
disp(find_solid_area(v1))
end
