vdata = readmatrix('original.txt'); %input vector
rotated_index = [21:-5:1,22:-5:2,23:-5:3,24:-5:4,25:-5:5];
vdata90 = zeros(1,100); %output vector, need to check colum or row
for i = 1:25
    cur_index = rotated_index(i);
    vdata90(i) = vdata(cur_index); % T/F
    vdata90(i*2+24) = vdata(cur_index*2+24); % l1 TODO rotate
    vdata90(i*2+25) = vdata(cur_index*2+25); % l2 TODO rotate
    vdata90(i+75) = -pi/2+vdata(cur_index+75); % theta TODO rotate
end
writematrix(vdata90','rotated.txt')