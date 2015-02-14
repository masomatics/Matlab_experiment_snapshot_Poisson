NN = 100000
hoge = rand(1,NN);
%rate = [2, 5, 4, 10]; 
rate = rand(1,4).^2;
cnt = [0,0,0,0];
    cumrate(1) = rate(1);
    cumrate(2) = cumrate(1) +rate(2);
    cumrate(3) = cumrate(2) +rate(3);
for k = 1:NN
    rxn = 1 + sum(cumrate < (hoge(k)*sum(rate))); 
    switch rxn
        case 1
            cnt(1) = cnt(1) + 1;
        case 2
            cnt(2) = cnt(2) + 1;
        case 3
            cnt(3) = cnt(3) + 1;
        case 4
            cnt(4) = cnt(4) + 1;
    end
end
        
cnt/ NN
rate/ sum(rate)