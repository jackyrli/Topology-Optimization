clc
index_array = 1:25;
void_indices = zeros(1,12);
v1 = zeros(1,25);
vdata = [1,2,3,4,6,6,7,8,9,10,11,12]
for i = 1: 12
    void_indices(i) = index_array(vdata(i));
    index_array(vdata(i)) = [];
end
v1(void_indices) = 1;
v1
void_indices
index_array