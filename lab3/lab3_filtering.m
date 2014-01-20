close all
clear all
clc
load data.mat

%% Practice 3

%% Exercise 1

var_s = 1;
var_n = 0.2;
N = 100;

s = [0.54 1.83 -2.26 0.86 0.32]';
randn('seed',0);rand('seed',0);

%% 1.1

un=randn(1,N);
c=[un(1);zeros(length(s)-1,1)];
noise=sqrt(var_n)*randn(1,N);
U = toeplitz(c, un);

x1=U'*s+noise';

xfilter=(filter(s,1,un)+noise)';

pause 
%% 1.2


P = inv(U*U'+(var_n/var_s)*eye(5,5));

mean_s =  P*U*x1 % mean for the guassian distrubutio of 
                 % the posterior p(s|x)
                 
cov_s = var_n*P  % covarianza for the guassian distrubutio of 
                 % the posterior p(s|x)
                
%% 1.3

ustar=[0;U(1:end-1,end)];
mu_star = ustar'*P*U*x1 % mean of the new output
v_star =var_n+var_n*ustar'*P*ustar % variance of the new output

%% 2.1

var_s = 0.2;
var_n = 4e-5;

xnew =voicein;
unew = voiceout;
Unew = zeros(7,length(unew));

count = 1;
for k=([3122 5953 9999 14999 18999 29295 39385]+1)
    Unew(count,1:k-1) = zeros(1,k-1);
    Unew(count,k:end) = unew(1:length(unew)-(k-1));
    count = count + 1;
end


Pnew = inv(Unew*Unew'+(var_n/var_s)*eye(7,7));
mean_s =  Pnew*Unew*xnew 
cov_s = var_n*Pnew

%% 2.2

filtered = Unew'*mean_s;
voiceremote = voicein-filtered;

fs= 22.05e3;

%% 2.3
sLMS = zeros(7,1);
msemin = 100000000;
mumin=0;


% this code finds mu value that minimizes the error (noise)
% mu = 6.3

% for mu=5:.1:7
%     sLMS = zeros(7,1);
%     voiceremoteLMS = zeros(100000,1);
%     
%     for k=1:100000
%         voiceremoteLMS(k) = voicein(k)-Unew(:,k)'*sLMS;       
%         sLMS = sLMS + mu*(voicein(k)-Unew(:,k)'*sLMS)*Unew(:,k);
%     end
%     
%     mseNew = mean((voiceremoteLMS).^2);
%     
%     
%     if mseNew<msemin
%         msemin=mseNew;
%         mumin=mu;
%     end
% end

mu = 6.3;
for k=1:100000
    voiceremoteLMS(k) = voicein(k)-Unew(:,k)'*sLMS;
    sLMS = sLMS + mu*(voicein(k)-Unew(:,k)'*sLMS)*Unew(:,k);
end

mumin
sLMS

% an objective way to see that sLMS is better than the batch
% verson is just to compute the MSE mean(voiceremote.^2)
% here we sea that voiceremoteLMS is smaller than voiceremose



