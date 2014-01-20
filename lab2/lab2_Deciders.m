%% Set up

 hold off
 clc
 clear all
 close all

alfa=5*10^(5);
c = 300000;

%r = c*t/2;

v= 10^-4;
r=0:15000;
eta = alfa./(2*r.^2); % solved for on paper
Xr = 1/(sqrt(2*pi));

%% Part 2.1

Pfa = 1-normcdf(eta,0,sqrt(v));

plot(r,Pfa,'k-');
xlabel('Distance from Radar (m)')
ylabel('Prbability')
title('Probability vs Distance from Radar')

Pd=1-normcdf(eta,alfa./r.^2,sqrt(v));

hold on
plot(r,Pd,'b-')
legend('Pm','Pd')

pause

%% Part 2.2

NPcnst = 10^(-3);
NPcnstnew = 1-NPcnst

NPThreshhold = norminv(NPcnstnew,0,sqrt(v))

Pfa2 = 1-normcdf(NPThreshhold,0,sqrt(v));
Pm2=normcdf(NPThreshhold,alfa./r.^2,sqrt(v));
Pd2 = 1-Pm2;

figure(2)
plot(r,Pd2,'k-')
xlabel('Distance from Radar (m)')
ylabel('Prbability of Detection')
title('P_D vs Distance from Radar with Neymam-Pearson')

location = find(Pd2>=0.9);
maxDistance = location(end) % this is the maximum distance

% distance=norminv(0.9,alfa./r.^2,sqrt(v));

pause

%% Part 2.3

figure(3)
PTest = linspace(0,1,1000); % take 1000 linearly spaced values for Pfa

etaNew = norminv(PTest,0,sqrt(v)); % now get values for eta that result
                                   % in the previously created set for Pfa
colors = {'k','r','b'}; % generates three different line colors to be used
count = 1;
for k = [2000 5000 10000]
    PfaNew = 1-normcdf(etaNew,0,sqrt(v));
    PdNew = 1-normcdf(etaNew,alfa./k.^2,sqrt(v));
    hold on
    h = plot(PfaNew,PdNew,colors{count});
    count = count + 1;
end

title('ROC Plot for Different Radius Values')
xlabel('P_{FA}');
ylabel('P_{D}');
legend('r = 2000', 'r = 5000', 'r = 10000',0);

% the following code expands the axis so the data is not on the edge
f=gca; n=0.1;
w=axis;
xr = (w(2)-w(1));
yr = (w(4)-w(3));
axis(f, [w + n*[-xr xr -yr yr]]);
hold off
%% Part 3.2

rmax = 15000;
count = 1;
figure
color = {'k','r','b','g','y','m'}
for l=[1 3 10 30 100 300]
    etaNew = norminv(PTest,0,sqrt(v/l));
    PfaNew = 1-normcdf(etaNew,0,sqrt(v/l));
    PdNew = 1-normcdf(etaNew,alfa./rmax.^2,sqrt(v/l));
    hold on
    h = plot(PfaNew,PdNew,color{count});
    count = count + 1;
end
title('ROC Plot for different l values')
xlabel('P_{FA}');
ylabel('P_{D}');
legend('l= 1', 'l = 3', 'l = 10', 'l = 30', 'l = 100', 'l = 300', 0);



% the following code expands the axis so the data is not on the edge
f=gca; n=0.1;
w=axis;
xr = (w(2)-w(1));
yr = (w(4)-w(3));
axis(f, [w + n*[-xr xr -yr yr]]);

%% Part 3.3

NPcnst = 10^(-3);
NPcnstnew = 1-NPcnst;
r = 10000;
l=[1:300]
NPThreshhold = norminv(NPcnstnew,0,sqrt(v./l));

Pfa2 = 1-normcdf(NPThreshhold,0,sqrt(v./l));
Pm2=normcdf(NPThreshhold,alfa./r.^2,sqrt(v./l));
Pd2 = 1-Pm2;

figure()
plot(l,Pd2)
xlabel('Number of measurements l')
ylabel('Prbability of Detection')
title('P_D vs l from Radar with Neymam-Pearson')

location = find(Pd2>=0.9);
lo = location(1) % this is the maximum distance

%% Part 3.4


NPcnst = 10^(-3);
NPcnstnew = 1-NPcnst;
r = 1:15000;
figure()
color = {'k','r'};
count = 1;
for l=[1, lo]
NPThreshhold = norminv(NPcnstnew,0,sqrt(v./l));

Pfa2 = 1-normcdf(NPThreshhold,0,sqrt(v./l));
Pm2=normcdf(NPThreshhold,alfa./r.^2,sqrt(v./l));
Pd2 = 1-Pm2;
hold on
plot(r,Pd2,color{count})
count = count +1;
end


xlabel('Distance r')
ylabel('Prbability of Detection')
legend('l= 1','l=l0');
title('P_D vs r from Radar with Neymam-Pearson')

%% Part 3.5

NPcnst = 10^(-3);
NPcnstnew = 1-NPcnst;
r = 10000;
factor = 0:0.001:100; 
count = 1;
Pd2new = -1000;

NPThreshhold = norminv(NPcnstnew,0,sqrt(v./lo));
Pfa2 = 1-normcdf(NPThreshhold,0,sqrt(v./lo));
Pm2=normcdf(NPThreshhold,alfa./r.^2,sqrt(v./lo));
Pd2 = 1-Pm2;
count =1;

while (abs(Pd2-Pd2new)>0.00001 && count<length(factor))
    NPThreshhold = norminv(NPcnstnew,0,sqrt(v));
    Pfa2new = 1-normcdf(NPThreshhold,0,sqrt(v));
    Pm2new=normcdf(NPThreshhold,factor(count).*alfa./r.^2,sqrt(v));
    Pd2new = 1-Pm2new;
    count = count +1;
end

Pd2new
Pd2
factor(count-1)

%% Part 4.1
load data_P2.mat
figure()
imagesc(X)
xlabel('Distance')

%% Part 4.2

T = zeros(100,1000);
for k=1:100
    
    T(k,:)= mean(X(1:k,:));
    
end

figure()
imagesc(T)
title('Value of T')
ylabel('Number of Observations')
xlabel('Distance from sensor (m)')
%% Part 4.3

v= 10^-4;
NPVect = zeros(100,1);

for k=1:100
    NPcnst = 10^(-3);
    NPcnstnew = 1-NPcnst;
    NPVect(k) = norminv(NPcnstnew,0,sqrt(v)/k);
end

% NPVect = NPVect+.00511; % this will increase the value of eta

NPMatrix = NPVect * ones(1,1000);

Decisions = (T>NPMatrix);

figure()
imagesc(Decisions)
title('Output of Decider (Red == 1)')
ylabel('Number of Observations')
xlabel('Distance from sensor (m)')
