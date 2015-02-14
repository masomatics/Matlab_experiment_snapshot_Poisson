%Example in matlab

rng('default');
mu1 = [1,2];
sigma1 = [3, .2;.2, 2];
mu2 = [-1,-2];
sigma2 = [2 0 ; 0 1];
datXY = [mvnrnd(mu1, sigma1, 200); mvnrnd(mu2, sigma2, 100)]


scatter(datXY(:,1), datXY(:,2), 10, 'ko');


options = statset('Display', 'final');
gm = gmdistribution.fit(datXY,2,'OPtions', options)

hold on
ezcontour(@(x,y)pdf(gm,[x y]), [-8,6], [-8, 6])
hold off;
gm.mu
gm.Sigma
gm.PComponents

%Snapshot, x and y

datamat = data_generation_N(1,1,1,1,1,1,100,true,true)
snapshot_y = snapshot_data(datamat, 45, true)
snapshot_x = snapshot_data(datamat, 45, false)

scatter(snapshot_x, snapshot_y)
hold on;
datXY = [snapshot_x', snapshot_y']
gm = gmdistribution.fit(datXY,2,'OPtions', options)
ezcontour(@(x,y)pdf(gm,[x y]), [-25,15], [-15, 30])
hold off;


gm.mu
gm.Sigma
gm.PComponents


%Snapshot, x 


gm = gmdistribution.fit(snapshot_x',2,'OPtions', options)
scatter(snapshot_x, pdf(gm,snapshot_x'))
gm.mu
gm.Sigma
gm.PComponents



for(k = 1:50)
    snapshot_x = snapshot_data(datamat, k, true)    ;
    
    
    MaxComp =2
    AIC = zeros(1,MaxComp);
    obj = cell(1,MaxComp);
    for p = 1:2
        obj{p} = gmdistribution.fit(snapshot_x',p);
        AIC(p)= obj{p}.AIC;
    end
    dbstop if warning
    [minAIC,numComponents] = min(AIC);
    display(['#components: ', num2str(numComponents)])

    gm = obj{numComponents};
    datsurf(1,:,k) = snapshot_x;
    datsurf(2,:,k) =pdf(gm,snapshot_x');
end


%Snapshot, y 


gm = gmdistribution.fit(snapshot_y',2,'OPtions', options)
scatter(snapshot_y, pdf(gm,snapshot_y'))
gm.mu
gm.Sigma
gm.PComponents


for(k = 1:50)
    snapshot_y = snapshot_data(datamat, k, true)    ;
    
    
    MaxComp =2
    AIC = zeros(1,MaxComp);
    obj = cell(1,MaxComp);
    for p = 1:2
        obj{p} = gmdistribution.fit(snapshot_y',p);
        AIC(p)= obj{p}.AIC;
    end
  %  dbstop if warning
    [minAIC,numComponents] = min(AIC);
    display(['#components: ', num2str(numComponents)])

    gm = obj{numComponents};
    datsurf(1,:,k) = snapshot_y;
    datsurf(2,:,k) =pdf(gm,snapshot_y');
end

