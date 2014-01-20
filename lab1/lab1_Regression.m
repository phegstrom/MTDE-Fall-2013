% lab1.m
% [Parker Hegstrom]
% [October 3, 2013]



clear
load DatosP1
format long

nTrain=length(xTrain);
nVal=380;
nTest=190;
mx = mean(xTrain); stdx = std(xTrain);
xTrain = (xTrain - ones(nTrain,1)*mx)./(ones(nTrain, 1)*stdx);
xVal = (xVal - ones(nVal,1)*mx)./(ones(nVal,1)*stdx);
xTest = (xTest - ones(nTest,1)*mx)./(ones(nTest,1)*stdx);


%% Create Figure

figure(1)

for k=1:9
    subplot(3,3,k)
    plot(xTrain(:,k),sTrain,'.')
    ylabel('s')
    xlabel(strcat('x_{',int2str(k),'}'))
end


%% Estimation

on = ones(380,1);
xTraine = [on,xTrain]; % add row of 1's to xTrain
Xe = (xTraine' * xTraine);
We = inv(Xe) * xTraine' * sTrain % column vec holding coefs
Wet = We';
xTrainet = xTraine';
sHAT = Wet*xTrainet;
MSE = mean((sTrain-sHAT').^2) % value for MSE using all vars


Wo=We(1,1)*ones(380,1);
MSEO = mean((sTrain-Wo).^2) % value for MSE using base estimator


%% Calculate for 2.b (without the Wo part)
idx = zeros(2^9-1,9);
MSEAll = zeros(511,1);
WeShort = We(2:end);


Best = [10000;10000;10000;10000;10000;10000;10000;10000;10000];
BestLogical = zeros(9,9);



for k = 1:2^9-1
    idx = bitget(k,1:9);
    numOfVars = sum(idx);
    idxLogical = logical(bitget(k,1:9));
    XTrainNew = [ones(nTrain,1) xTrain(:,idxLogical)];
    WeNew = [We(1) ; WeShort(idxLogical)];
    sHatNew = WeNew' * XTrainNew';
    MSE = mean((sTrain-sHatNew').^2);
    if MSE < Best(numOfVars)
        Best(numOfVars) = MSE;
        BestLogical(numOfVars,:) = idx;
    end
    
end



figure(2)

t=1:9;
for k=1:9
    BestLogical(k,:) = BestLogical(k,:).*t;
end

plot(t, BestLogical,'.')
axis([1 9 1 9])
ylabel('Selected Variables')
xlabel('Number of Variables Used')

% note can also use spy(matrix with ones and 0s) to make the figure


%% Part 3

x1Train = xTrain(:,1);
x1TrainAll = [x1Train x1Train.^2 x1Train.^3 x1Train.^4 ...
    x1Train.^5 x1Train.^6 x1Train.^7 x1Train.^8];
x1WeAll = zeros(8,8);
x1MSEAll = zeros(8,1); %%%% now fill this out


for k=1:8
    x1TrainNew = [on x1TrainAll(:,1:k)];
    x1WeNew = inv(x1TrainNew' * x1TrainNew) * x1TrainNew' * sTrain;
    for j=1:length(x1WeNew)
        x1WeAll(j,k) = x1WeNew(j);
    end
    sHat = x1WeNew' * x1TrainNew';
    x1MSE = mean((sTrain - sHat').^2);
    x1MSEAll(k) = x1MSE;
end



% x1TrainOnes = on;
% x1WeOnes = inv(x1TrainOnes' * x1TrainOnes) * x1TrainOnes' * sTrain;
% sHatOnes =

t = 0:8;

figure(3)
plot(t,[MSEO;x1MSEAll],'o')             % here concatenate the MSEO val
title('MSE vs. Degree of Polynomial')
xlabel('Degree of Polynomial')
ylabel('MSE')


%% Part 3 b

% use x1Train(randperm(n,k))

ALLMSEVal = zeros(8,1000);
ALLMSETrain = zeros(8,1000);
allValToPlot = zeros(8,8);
allTrainToPlot = zeros(8,8);

% choosing i number of samples to make the model

count = 1;
for i=20:20:160

    
    for n=1:1000 % because random, we need to average over 1000 trials
        randSample = randperm(380,i); % getting i random numbers
        x1Train = xTrain(randSample,1);
        sTrainRand = sTrain(randSample,1);
        x1TrainAll = [x1Train x1Train.^2 x1Train.^3 x1Train.^4 ...
            x1Train.^5 x1Train.^6 x1Train.^7 x1Train.^8];
        x1WeAll = zeros(8,8);
        x1MSETrainAll = zeros(8,1);
        x1MSEValAll = zeros(8,1); %%%% now fill this out
        x1Val = xVal(:,1);
        x1ValAll = [x1Val x1Val.^2 x1Val.^3 x1Val.^4 x1Val.^5 x1Val.^6 ...
            x1Val.^7 x1Val.^8];
        
        
        % determine coefs for different order of function
        for k=1:8
            x1TrainNew = [on(1:i) x1TrainAll(:,1:k)]; % concatenating column of 1's
            x1ValNew = [on x1ValAll(:,1:k)];
            x1WeNew = inv(x1TrainNew' * x1TrainNew) * x1TrainNew' * sTrainRand;
            
            x1WeAll(1:length(x1WeNew),k) = x1WeNew;
            
            sHat = x1WeNew' * x1TrainNew';
            
            sHatVal = x1WeNew' * x1ValNew';
            x1MSETrain = mean((sTrain(randSample) - sHat').^2);
            x1MSEVal = mean((sVal - sHatVal').^2);
            x1MSETrainAll(k) = x1MSETrain;
            x1MSEValAll(k) = x1MSEVal;
        end
        
        ALLMSETrain(:,n) = x1MSETrainAll;
        ALLMSEVal(:,n) = x1MSEValAll;
        
    end
    
    averageMSEVal = zeros(8,1);
    averageMSETrain = zeros(8,1);
    
    for k=1:8
        averageMSEVal(k) = mean(ALLMSEVal(k,:));
        averageMSETrain(k) = mean(ALLMSETrain(k,:));
    end
    
    % these vectors contain the avg MSE calculated for each value of i
    % (the number of random samples chosen)
    allValToPlot(:,count) =  averageMSEVal;
    allTrainToPlot(:,count) = averageMSETrain;
    count = count +1;
end


% Create the last figure 

figure(4)
K = 1:1:8;
h = plot(K',allTrainToPlot, '-', K', allValToPlot, '--')
% set marker styles
set(h(9),'Marker','+');
set(h(10),'Marker','o');
set(h(11),'Marker','*');
set(h(12),'Marker','x');
set(h(13),'Marker','s');
set(h(14),'Marker','d');
set(h(15),'Marker','p');
set(h(16),'Marker','h');
% create legend
legend(h(9:16),'K = 20','K = 40','K = 60','K = 80',...
                'K = 100','K = 120','K = 140','K = 160')
title('MSE vs. Degree of the Model')
xlabel('Degree of Model')
ylabel('MSE')
axis([1 8 0 50])


