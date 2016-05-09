function list = list_of_all(m)
% usage: list = list_of_all(m)
% simple function to produce of list of all pairs of m elements

nk = m*(m-1)/2;
list = zeros(nk,2);
k = 1;
for i=1:m
    for j=(i+1):m
        list(k,1) = i;
        list(k,2) = j;
        k = k+1;
    end
end
