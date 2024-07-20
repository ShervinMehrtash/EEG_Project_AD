% This MATLAB .m file contains all codes written for Signals and Systems
% final project --------------------------------------------------------
% Student Name  =  Shervin Mehrtash  -----------------------------------
% Student Number  =  400102052  ----------------------------------------



% Preprocessing for Subject 1 ---->
% ----------------------------------------------------------------------



%% Epoching the Data -->

Data = EEG.data; % Loading the Data from EEG Workspace
Start = 14*200; % Starting sample - from 14 sec after start of the event
End = 1214*200; % Ending Sample - Neglecting the data after 120 trials
Clean_Data = Data(:, Start:End);

%% Loop to fill the 3D epoch Matrix -->

Epoch = zeros(19,600,120);

for i = 1:19
    for j = 1:120
        Epoch( i , : , j) = Clean_Data(i, (j-1)*2000 + 601 : (j-1)*2000 + 1200);
    end
end

%% Deleting the unnecessary channels and keeping the subsamples -->

Cleaned_Epoch = EEG.data;
Subsampled_Epoch = Cleaned_Epoch([1, 5, 10, 15], : , : ); % The Final Epoch

% Converting the Subsampled Epoch array to struct -->
Epoch_Struct = struct('Epoch_Subsampled', Subsampled_Epoch);






% Preprocessing for Subject 2 ---->
% ----------------------------------------------------------------------



%% Epoching the Data -->

Data = EEG.data; % Loading the Data from EEG Workspace
Start = 14*200; % Starting sample - from 14 sec after start of the event
End = 1214*200; % Ending Sample - Neglecting the data after 120 trials
Clean_Data = Data(:, Start:End);

%% Loop to fill the 3D epoch Matrix -->

Epoch = zeros(19,600,120);

for i = 1:19
    for j = 1:120
        Epoch( i , : , j) = Clean_Data(i, (j-1)*2000 + 601 : (j-1)*2000 + 1200);
    end
end

%% Deleting the unnecessary channels and keeping the subsamples -->

Cleaned_Epoch = EEG.data;
Subsampled_Epoch = Cleaned_Epoch([1, 5, 10, 15], : , : ); % The Final Epoch

% Converting the Subsampled Epoch array to struct -->
Epoch_Struct = struct('Epoch_Subsampled', Subsampled_Epoch);






% Process Sections (From section 4 on) ---->
%-----------------------------------------------------------------------



%% Initial Commands
clear; clc;

%% Loading the required Datasets
load('Codes\AD.mat');
load('Codes\normal.mat');
load('Codes\MCI.mat');

%% Eliciting the Epoch for each patient
Normal_Data = zeros(15,2); % First column for Frequent and Second for Rare
AD_Data = zeros(13,2); % First column for Frequent and Second for Rare
MCI_Data = zeros(7,2); % First column for Frequent and Second for Rare

% Loop for narmal dataset
for n = 1:15
    Epoch = normal(n).epoch;
    Odor = normal(n).odor;
    Normal_Data(n, : ) = PLV_Calculator(Epoch, Odor);
end

% Loop for AD dataset
for n = 1:13
    Epoch = AD(n).epoch;
    Odor = AD(n).odor;
    AD_Data(n, : ) = PLV_Calculator(Epoch, Odor);
end

% Loop for MCI dataset
for n = 1:7
    Epoch = MCI(n).epoch;
    Odor = MCI(n).odor;
    MCI_Data(n, : ) = PLV_Calculator(Epoch, Odor);
end

%% Finding the Boxplots for each group and each Odor
figure
subplot(2,2,1)
boxplot(Normal_Data(:,1));
title('PLV Data Summary for Normal patients exposed to Frequent Odor','Interpreter','latex');
subplot(2,2,2)
boxplot(Normal_Data(:,2));
title('PLV Data Summary for Normal patients exposed to Rare Odor','Interpreter','latex');
subplot(2,2,3)
boxplot(AD_Data(:,1));
title('PLV Data Summary for AD patients exposed to Frequent Odor','Interpreter','latex');
subplot(2,2,4)
boxplot(AD_Data(:,2));
title('PLV Data Summary for AD patients exposed to Rare Odor','Interpreter','latex');

%% Plotting the Distribution of PLVs for 4 different groups
Fit_Dist1 = fitdist(Normal_Data(:,1),"Normal");
Fit_Dist2 = fitdist(Normal_Data(:,2),"Normal");
Fit_Dist3 = fitdist(AD_Data(:,1),"Normal");
Fit_Dist4 = fitdist(AD_Data(:,2),"Normal");
X = -0.5:0.01:2;
PDF1 = pdf(Fit_Dist1,X); PDF2 = pdf(Fit_Dist2,X); PDF3 = pdf(Fit_Dist3,X); PDF4 = pdf(Fit_Dist4,X);
figure
subplot(2,2,1)
plot(X,PDF1);
title('Normal Distribution fit to the PLV (Normal patient - Frequent Odor)','Interpreter','latex');
subplot(2,2,2)
plot(X,PDF2);
title('Normal Distribution fit to the PLV (Normal patient - Rare Odor)','Interpreter','latex');
subplot(2,2,3)
plot(X,PDF3);
title('Normal Distribution fit to the PLV (AD patient - Frequent Odor)','Interpreter','latex');
subplot(2,2,4)
plot(X,PDF4);
title('Normal Distribution fit to the PLV (AD patient - Rare Odor)','Interpreter','latex');

%% Finding the P-Values --->
[h1,p1] = ttest2(Normal_Data(:,1),AD_Data(:,1));
[h2,p2] = ttest2(Normal_Data(:,2),AD_Data(:,2));

%% Finding Phase Difference --->
% Choosing a Random Guy (Number 12 in normal and 6 in AD patients)
Epoch_Random_Normal = normal(2).epoch;
Odor_Random_Normal = normal(2).odor;
Epoch_Random_AD = AD(1).epoch;
Odor_Random_AD = AD(1).odor;

% Finding the number of frequent odors for each group
FreqNum_Normal = size(Odor_Random_Normal,1) - sum(Odor_Random_Normal);
FreqNum_AD = size(Odor_Random_AD,1) - sum(Odor_Random_AD);

% Defining the arrays to store the value of Phase Difference for each group
PhaseDiff_Normal = zeros(FreqNum_Normal,1);
PhaseDiff_AD = zeros(FreqNum_AD,1);

Counter1 = 0; Counter2 = 0;

% Loop to find and store the phase difference of normal random guy
for i = 1 : size(Odor_Random_Normal,1)
    if Odor_Random_Normal(i,1) == 0
        Counter1 = Counter1 + 1;
        x1 = Epoch_Random_Normal(2,:,i);
        x2 = Epoch_Random_Normal(3,:,i);
        Result = PDFinder(x1,x2);
        PhaseDiff_Normal(Counter1,1) = Result;
    end
end

% Loop to find and store the phase difference of AD random guy
for i = 1 : size(Odor_Random_AD,1)
    if Odor_Random_AD(i,1) == 0
        Counter2 = Counter2 + 1;
        x1 = Epoch_Random_AD(2,:,i);
        x2 = Epoch_Random_AD(3,:,i);
        Result = PDFinder(x1,x2);
        PhaseDiff_AD(Counter2,1) = Result;
    end
end

% Plotting the polar histogram
figure
subplot(1,2,1)
polarhistogram(PhaseDiff_Normal);
title('Polar Histogram for Phase Difference of a Randomly selected Normal patient','Interpreter','latex');
subplot(1,2,2)
polarhistogram(PhaseDiff_AD);
title('Polar Histogram for Phase Difference of a Randomly selected A.D. patient','Interpreter','latex');

%% Finding the Mean value of Phase Difference for each group --->
% Defining the arrays to store the mean value of Phase Difference
Mean_PhaseDiff_Normal = zeros(15,1);
Mean_PhaseDiff_AD = zeros(13,1);
for n = 1 : 15
    Mean_PhaseDiff_Normal(n,1) = PDPP(normal(n).epoch, normal(n).odor);
end

for n = 1 : 13
    Mean_PhaseDiff_AD(n,1) = PDPP(AD(n).epoch, AD(n).odor);
end

Mean_Normal = sum(Mean_PhaseDiff_Normal) / size(Mean_PhaseDiff_Normal,1);
Mean_AD = sum(Mean_PhaseDiff_AD) / size(Mean_PhaseDiff_AD,1);

% Plotting the polar histogram
figure
subplot(1,2,1)
polarhistogram(Mean_PhaseDiff_Normal);
title('Polar Histogram for Mean of Phase Difference in Normal dataset','Interpreter','latex');
subplot(1,2,2)
polarhistogram(Mean_PhaseDiff_AD);
title('Polar Histogram for Mean of Phase Difference in AD dataset','Interpreter','latex');

%% Finding the PLV between pair of channels

% For this section, I used the code that finds the PLV for a given pair
% of channels, And Manually calculate the PLV and store them in the below
% arrays ---->

% Code Used to find the proper values of PLV -->

% S = sum(Normal_Data(:,1)); % Or (:,2) for Rare Odor exposure
% S = S/15; % Or S/13 For AD group

% Normal Group Data -->
Normal_PLV = [
    0.5897 0.5920;
    0.5550 0.5616;
    0.6175 0.6219;
    0.8064 0.8071;
    0.8169 0.8189;
    0.7969 0.7973
    ];

% AD Group Data -->
AD_PLV = [
    0.6029 0.6026;
    0.5467 0.5488;
    0.6253 0.6264;
    0.6458 0.6362;
    0.7978 0.8008;
    0.7254 0.7235
    ];

% Cleared Data for the Heatmap -->
Final_PLV = [
    0.5897 0.6029;
    0.5550 0.5467;
    0.6175 0.6253;
    0.8064 0.6458;
    0.8169 0.7978;
    0.7969 0.7254;
    0.5920 0.6026;
    0.5616 0.5488;
    0.6219 0.6264;
    0.8071 0.6362;
    0.8189 0.8008;
    0.7973 0.7235
    ];

% Plotting the heatmap of the elicited data
XAxis = {'Normal Group','AD Group'};
YAxis = {'Fp1/Fz - Frequent','Fp1/Cz - Frequent','Fp1/Pz - Frequent','Fz/Cz - Frequent','Fz/Pz - Frequent','Cz/Pz - Frequent','Fp1/Fz - Rare','Fp1/Cz - Rare','Fp1/Pz - Rare','Fz/Cz - Rare','Fz/Pz - Rare','Cz/Pz - Rare'};
figure
heatmap(XAxis, YAxis, Final_PLV);


% Finding the P-Values for Channel
% This section uses the function of other previous parts, so the code is
% commented to avoid any violence in other sections of the code -->

% [h3,p3] = ttest2(Normal_Data(:,1),AD_Data(:,1));
% [h4,p4] = ttest2(Normal_Data(:,2),AD_Data(:,2));

h3 = 0; h4 = 0; % Elicited Values for h of channels Cz and Pz
p3 = 0.3647; p4 = 0.3733; % Elicited Values for p-value of channels Cz and Pz

%% ---------- BONUS SECTION 5.1 ----------
% As mentioned earlier, the PLV values of Two different pair of channels
% are found manually, to avoid any violence in the code :)

% Code Used to find the proper values of PLV for MCI -->
% S = sum(MCI_Data(:,1));
% S = S/7;

% Cleared Data for the Heatmap of pair of channels (2,3 and 3,4) -->
Final_PLV_3State = [
    0.8064 0.7087 0.6458;
    0.8071 0.6995 0.6362;
    0.7969 0.7124 0.7254;
    0.7973 0.7105 0.7235
    ];

% Plotting the heatmap of the elicited data
XAxis = {'Normal Group','MCI Group','AD Group'};
YAxis = {'Fz/Cz - Frequent','Fz/Cz - Rare','Cz/Pz - Frequent','Cz/Pz - Rare'};
figure
heatmap(XAxis, YAxis, Final_PLV_3State);

% Elicited values for p-value
% Similarly, we used the code of the previous sections to avoid any
% violation
pVal1 = 0.3221; % Normal/MCI - Frequent - Fz/Cz
pVal2 = 0.5154; % MCI/AD - Frequent - Fz/Cz
pVal3 = 0.3015; % Normal/MCI - Rare - Fz/Cz
pVal4 = 0.5347; % MCI/AD - Rare - Fz/Cz
pVal5 = 0.3614; % Normal/MCI - Frequent - Cz/Pz
pVal6 = 0.8908; % MCI/AD - Frequent - Cz/Pz
pVal7 = 0.3656; % Normal/MCI - Rare - Cz/Pz
pVal8 = 0.8935; % MCI/AD - Rare - Cz/Pz


%% ---------- BONUS SECTION 5.2 ----------
% Finding the MVL for signals from channels Fz and Cz -->
% Defining vectors to store the MVLs
MVL_Normal = zeros(15,1);
MVL_MCI = zeros(7,1);
MVL_AD = zeros(13,1);

% Loop for normal dataset
for n = 1:15
    Epoch = normal(n).epoch;
    Odor = normal(n).odor;
    MVL_Normal(n,1) = MVL_Calculator(Epoch, Odor);
end

% Loop for MCI dataset
for n = 1:7
    Epoch = MCI(n).epoch;
    Odor = MCI(n).odor;
    MVL_MCI(n, 1) = MVL_Calculator(Epoch, Odor);
end

% Loop for AD dataset
for n = 1:13
    Epoch = AD(n).epoch;
    Odor = AD(n).odor;
    MVL_AD(n, 1) = MVL_Calculator(Epoch, Odor);
end
Sum1 = sum(MVL_Normal/size(MVL_Normal,1));
Sum2 = sum(MVL_MCI/size(MVL_MCI,1));
Sum3 = sum(MVL_AD/size(MVL_AD,1));
%--------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------



%% Function used to find the Mean value for Phase Difference --->

function PD = PDPP(Epoch_Random_Normal, Odor_Random_Normal)
% Finding the number of frequent odors 
FreqNum_Normal = size(Odor_Random_Normal,1) - sum(Odor_Random_Normal);

% Defining the arrays to store the value of Phase Difference
PhaseDiff_Normal = zeros(FreqNum_Normal,1);

Counter1 = 0;

% Loop to find and store the phase difference of normal random guy
for i = 1 : size(Odor_Random_Normal,1)
    if Odor_Random_Normal(i,1) == 0
        Counter1 = Counter1 + 1;
        x1 = Epoch_Random_Normal(2,:,i);
        x2 = Epoch_Random_Normal(3,:,i);
        Result = PDFinder(x1,x2);
        PhaseDiff_Normal(Counter1,1) = Result;
    end
end

PD = sum(PhaseDiff_Normal) / size(PhaseDiff_Normal,1);

end

%% Function used to find the Phase Difference --->
function PhaseDiff = PDFinder(x1,x2)
    dot_product = dot(x1, x2);
    norm_product = (norm(x1) * norm(x2));
    PhaseDiff = acos(dot_product / norm_product);
end


%% Function used to find the PLV value --->
function PLV = PLV_Finder(data, channel1, channel2)

    % Extract the time series data for the two channels
    signal1 = data(channel1, :);
    signal2 = data(channel2, :);

    % Compute the Fourier transform of the signals
    fftSignal1 = fft(signal1);
    fftSignal2 = fft(signal2);

    % Define the frequency axis for FFT
    freqAxis = (0:length(signal1)-1)*(200/length(signal1));

    % Find indices corresponding to the desired frequency range
    freqIndices = find(freqAxis >= 35 & freqAxis <= 40);

    % Extract the phase angles for each frequency bin within the range
    phaseSignal1 = angle(fftSignal1(freqIndices));
    phaseSignal2 = angle(fftSignal2(freqIndices));

    % Compute PLV as mean of phase differences across frequencies
    PLV = abs(mean(exp(1i * (phaseSignal1 - phaseSignal2))));
end




%% Function to process and manage the storing of PLVs --->
function PLVPP = PLV_Calculator(Epoch, Odor)
    Num_of_Frequent = size(Odor,1) - sum(Odor);
    
    % Loop to find the PLV for each trial -->

    % Defining an array to store the PLV values for each patient
    PLV_Array = zeros(1,size(Epoch,3));

    for i = 1 : size(Epoch,3)
        PLV_Array(1,i) = PLV_Finder(Epoch(:,:,i),2,3);
    end

    % Defining vars to store the sum of different exposures
    Sum_Frequent = 0;
    Sum_Rare = 0;

    % Dividing Frequent and Rare exposures of each trial
    for i = 1:size(Odor,1)
        if Odor(i,1) == 0
            Sum_Frequent = Sum_Frequent + PLV_Array(1,i);
        end
        if Odor(i,1) == 1
            Sum_Rare = Sum_Rare + PLV_Array(1,i);
        end
    end

    % Finding the Average PLV for Frequent and Rare exposures
    PLV_Frequent = Sum_Frequent / Num_of_Frequent;
    PLV_Rare = Sum_Rare / (size(Odor,1) - Num_of_Frequent);
    PLVPP = [PLV_Frequent, PLV_Rare];
end




%% Function to Find the MVL of two signals --->
function MVL = MVLFinder(data, channel1, channel2)
% Extract the time series data for the two channels
signal1 = data(channel1, :);
signal2 = data(channel2, :);

% Compute the analytic signals of the two signals
analytic_signal1 = hilbert(signal1);
analytic_signal2 = hilbert(signal2);

% Compute the pairwise multiplication of the analytic signals of the two signals
pairwise_multiplication = analytic_signal1 .* conj(analytic_signal2);

% Compute the mean vector resulting from the pairwise multiplication of the analytic signals of the two signals
mean_vector = mean(pairwise_multiplication);

% Compute the length of the mean vector resulting from the pairwise multiplication of the analytic signals of the two signals
MVL = abs(mean_vector);

end

%% Function to process and manage the storing of MVLs --->
function MVLPP = MVL_Calculator(Epoch, Odor)    
    % Loop to find the PLV for each trial -->
    % Defining an array to store the PLV values for each patient
    MVL_Array = zeros(1,size(Epoch,3));

    for i = 1 : size(Epoch,3)
        MVL_Array(1,i) = MVLFinder(Epoch(:,:,i),2,3);
    end

    % Defining vars to store the sum
    Sum = 0;

    % Finding the sum
    for i = 1:size(Odor,1)
        Sum = Sum + MVL_Array(1,i);
    end

    % Finding the Average MVL for the person
    MVLPP = Sum / size(Odor,1);
end