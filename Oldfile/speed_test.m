N = 5000; 
rep= 4500;
rxn_mat =    [1     0       -1      0       0; 
                  0     1       0       -2      2;
                  0     0       0       1       -1]; 
              
piyo= kron(rxn_mat, eye(N));

hoge = rand([N*5,1]);
tic,
hh = waitbar(0, '‚Ç‚ç‚¦‚à‚ñ')
for k = 1:rep
    waitbar(k/rep) 
    piyo* hoge;
end 
toc