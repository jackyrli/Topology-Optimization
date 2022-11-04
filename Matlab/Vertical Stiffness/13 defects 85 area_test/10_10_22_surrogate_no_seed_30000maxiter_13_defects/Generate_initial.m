clear all;
x = zeros(12409,64);
y = zeros(12409,1);
for i = 0:12408
    mat = readmatrix(['sim',num2str(i),'.txt']);
    y(i+1) = mat(101);
    first25 = mat(1:25);
    second = mat(26:end);
    l1 = second(1:2:50);
    l2 = second(2:2:50);
    theta = second(51:75);
    
    index = find(first25);
    l1 = l1(index);
    l2 = l2(index);
    l3 = theta(index);
    
    x(i+1,1:25) = first25;
    x(i+1,26:end) = [l1' l2' l3'];
end
