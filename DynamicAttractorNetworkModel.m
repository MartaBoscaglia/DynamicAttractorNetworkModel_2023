%% Code related to "A dynamic attractor network model of memory formation, reinforcement and forgetting"

% This code allows to implement the model described in the paper "A dynamic
% attractor network model of memory formation, reinforcement and
% forgetting" by Marta Boscaglia, Chiara Gastaldi, Wulfram Gerstner and
% Rodrigo Quian Quiroga.

%This specific script applies to the case of a single-memory network. For
%additional information about this script and regarding the modifications to be applied in order to implement the other
%cases (e.g. network with 2 patterns) read the file README.txt, stored in the same repository of this
%script.

%This script mainly consists of three sections:
%1. Definition of parameters and directories (to be checked before running the script);
%2. Evolution of the network dynamics;
%3. Saving of the variables for further analysis.
%(Note that, during the evolution of the network dynamics, the script automatically
%saves the weight matrix at certain timesteps previous to the stimulation phase, during the
%stimulation phase (at each stim onset) and after the stimulation phase)

%%

clear all
close all
clc

%% Definition of parameters and directories

dirmain=''; %folder in which you want to save the results of the simulation

SavingResults=1; %1 if you want to save, in dirmain, the matlab variables at the end of the simulation (weight matrices are automatically saved during the simulation, independently of the value of SavingResults)

cd(dirmain)

% Defining parameters

HowManyStimNeurons=10; %number of stimulated neurons

N=100; %size of the network

dt=0.1; %timestep

r_max=1; %maximal firing rate
r0=0; %baseline firing rate
tau=1; %timescale determining the neuronal activation
h0=0.15; %base firing threshold in the absence of firing
b=100; %slope parameter of the sigmoidal transfer function

Wmax=3/HowManyStimNeurons; %maximal synaptic weight
Wmin=-0.5/HowManyStimNeurons; %minimal synaptic weight

Noise_Factor_RATE=0.2/HowManyStimNeurons; %constant defining the noise
C_noise=sqrt(dt/tau)*Noise_Factor_RATE; %coupling constant of the Gaussian noise. It contains a factor sqrt(dt/tau) as explained in the description of eq. 6 in the manuscript

SaturationWeights=1; %set to 1 for having hard bounds [Wmin, Wmax]
NoNegativeWeights=0; %1 if you want Wmin=0
autoconnections=0; % if 1 then a neuron can be connected to itself (i.e. diagonal of the weight matrix is not all zeros); if 0 then the diagonal of the weight matrix is all zeros

tmax=480000; %total time of the simulation (time is in a.u.) %note that, in this script settings, we consider 50000 a.u. of pre-stimulation phase + 420000 a.u. of stimulation phase + 10000 a.u. of after-stimulation phase
meanWin=15/dt; %length of the window adopted for calculating the learning rule's running average (This is in a.u.*10)

Factor_SR=2; %divisive normalization constant
Factor_SW=1; %synaptic normalization constant
SR_thr=-10000; %this value is very low because the heaviside function, used later in the code, does not have any effect on the divisive normalization
SW_thr=Wmax/6; %synaptic normalization depends on the summed weight of the strong connections

tau_teta=7*tau; %neural adaptation time constant
D_teta=1; %constant determining the adaptation strength

LR=1; %learning rate
beta=0.0025; %forgetting rate
tau_w=50; %timescale determining the learning

% Definition of the external stimulation and of the neurons to be
% stimulated

I=zeros(N,1);

vector_neurons=[];
vector_neurons=1:N;
vector_neurons=vector_neurons(randperm(length(vector_neurons)));

neurons_assembly=vector_neurons(1:HowManyStimNeurons); %vector defining the neurons belonging to the stimulated assembly
neurons_NOTassembly=vector_neurons((HowManyStimNeurons+1):end); %vector defining the neurons that will not be directly stimulated

I1=I;
I1(neurons_assembly)=1; %vector of the external stimulation for stimulating "neurons_assembly"

%Definition of the number, time and duration of the stimulations

HowManyStims=7000; %number of stimulations to "neurons_assembly"

tonset_first=50000; %time (in a.u.) of the first stimulation onset to "neurons_assembly"
toffset_first=50005; %time of the first stimulation offset to "neurons_assembly"

tsep=60; %time distance (in a.u.) between two consecutive stimulation to "neurons_assembly"

tonsets=[];
toffsets=[];

count=1;
for i=1:HowManyStims
    eval(['tonsets(' num2str(i) ')=tonset_first+(' num2str(count) '-1)*' num2str(tsep) ';']);
    eval(['toffsets(' num2str(i) ')=toffset_first+(' num2str(count) '-1)*' num2str(tsep) ';']);
    count=count+1;
end


% definition of variables for saving the weight matrices after the end of the stimulation phase
% (this is to check the dynamics of the network after the stimulation has
% finished)

HowManyAfterStim=100;
IntervalTest_AfterStim=100; %in a.u. %note that, with these values, we are considering 10000 a.u. of time after stimulation (100x100 a.u.)


t_AfterStim=[];

for i=1:HowManyAfterStim

    if i==1
        eval(['t_AfterStim(' num2str(i) ')=toffsets(end)+' num2str(IntervalTest_AfterStim) ';']);
    else
        eval(['t_AfterStim(' num2str(i) ')=t_AfterStim(' num2str(i) '-1)+' num2str(IntervalTest_AfterStim) ';']);
    end
    count=count+1;

end

%t_AfterStim is a vector of times defining when, after the stimulation phase, we will automatically save the weight matrix

% Definition of the variables for saving the results of the simulation

mean_weight_outsideTOinside=[];
mean_weights_assembly=[];
mean_weights_NOTassembly=[];
std_weights_assembly=[];
std_weights_NOTassembly=[];

r_ALL=zeros(N,(tmax/dt+1));
Input_saved_WithI_WithSc=zeros(N,(tmax/dt+1));

SF_time=zeros(N,(tmax/dt+1));
SF_time_weights=zeros(N,(tmax/dt+1));
SF_time_rates=zeros(N,(tmax/dt+1));

teta_adapt_ALL=zeros(N,(tmax/dt+1));


%% Evolution of the network dynamics

teta_adapt_before=ones(N,1)*h0;

W_before=zeros(N,N);

Tot_Synapsis_before=sum(W_before,2);

r_ALL(:,1)=ones(N,1)*r0;

Input_saved_WithoutI_before=r_ALL(:,1).*0;
Input_saved_WithI_before=r_ALL(:,1).*0;

count_WM=1;

WM_saved=[];

teta_adapt_ALL(:,1)=teta_adapt_before;

D_teta_temp=[];
D_teta_temp=D_teta;

t_temp=dt;
count_temp=2;

count_temp_first=count_temp;
t_temp_first=t_temp;
countForSavingWeightMatrix=0;

countMatSaved=0;


while t_temp<tmax

    countForSavingWeightMatrix=countForSavingWeightMatrix+1;

    teta_adapt=zeros(N,1);

    %count_temp %remove comment if you want to display the iteration number in the command window

    flag=0;
    checkflag=0;
    h=zeros(N,1);
    fi=zeros(N,1);


    if ~isempty(find(tonsets<=t_temp&toffsets>t_temp)) %checking if we are in the time in which we give stimulation to "neurons_assembly"
        if count_temp>2

            SF=[];
            SF_rates=[];
            SF_weights=[];

            r_forscaling=[];
            r_withoutnn=[];

            r_forscaling=r_ALL(:,count_temp-1);

            r_forscaling_heaviside=[];
            r_forscaling_heaviside=r_forscaling.*heaviside(r_forscaling-SR_thr);
            r_forscaling=[];
            r_forscaling=r_forscaling_heaviside;

            r_withoutnn=sum(r_forscaling)*ones(N,1);
            r_withoutnn=r_withoutnn-(r_ALL(:,count_temp-1).*heaviside(r_ALL(:,count_temp-1)-SR_thr));

            W_before_heaviside=[];
            W_before_heaviside=W_before.*heaviside(W_before-SW_thr);

            sum_IncomingW=[];
            sum_IncomingW=sum(W_before_heaviside,2);

            SW_denominator=[];
            SW_denominator=1+(sum_IncomingW).*Factor_SW;

            SR_denominator=[];
            SR_denominator=1+((r_withoutnn)./HowManyStimNeurons).*Factor_SR;

            SF_weights=1./sqrt(SW_denominator);
            SF_rates=1./(SR_denominator);
            SF=SF_weights.*SF_rates;
        else
            SF=ones(N,1);
            SF_rates=ones(N,1);
            SF_weights=ones(N,1);
        end

        h=(dotprod(W_before,r_ALL(:,count_temp-1))+I1).*SF;

        Input_saved_WithI_WithSc(:,count_temp)=(dotprod(W_before,r_ALL(:,count_temp-1))+I1).*SF;

        SF_time(:,count_temp)=SF;
        SF_time_rates(:,count_temp)=SF_rates;
        SF_time_weights(:,count_temp)=SF_weights;

        flag=1;
        checkflag=checkflag+1;

        if checkflag>1
            error('there is an error!')
        end
    end

    noise=[];
    % noise=(sqrt(dt/tau))*Noise_Factor_RATE*randn(N,1);
    noise=C_noise*randn(N,1);

    r_i_sub=[];


    if flag==0 %if flag is 0 it means that, at the current time, there is no external stimulation (so, the vector 'I' will be used)
        if count_temp>2
            SF=[];
            SF_rates=[];
            SF_weights=[];

            r_forscaling=[];
            r_withoutnn=[];

            r_forscaling=r_ALL(:,count_temp-1);

            r_forscaling_heaviside=[];
            r_forscaling_heaviside=r_forscaling.*heaviside(r_forscaling-SR_thr);
            r_forscaling=[];
            r_forscaling=r_forscaling_heaviside;

            r_withoutnn=sum(r_forscaling)*ones(N,1);
            r_withoutnn=r_withoutnn-(r_ALL(:,count_temp-1).*heaviside(r_ALL(:,count_temp-1)-SR_thr));

            W_before_heaviside=[];
            W_before_heaviside=W_before.*heaviside(W_before-SW_thr); %CHANGED WHEN WMAX=0.3

            sum_IncomingW=[];
            sum_IncomingW=sum(W_before_heaviside,2);

            SW_denominator=[];
            SW_denominator=1+(sum_IncomingW).*Factor_SW;

            SR_denominator=[];
            SR_denominator=1+((r_withoutnn)./HowManyStimNeurons).*Factor_SR;

            SF_weights=1./sqrt(SW_denominator);

            SF_rates=1./(SR_denominator);
            SF=SF_weights.*SF_rates;
        else
            SF=1;
            SF_rates=1;
            SF_weights=1;
        end

        h=(dotprod(W_before,r_ALL(:,count_temp-1))+I).*SF;

        Input_saved_WithI_WithSc(:,count_temp)=(dotprod(W_before,r_ALL(:,count_temp-1))+I).*SF;

        SF_time(:,count_temp)=SF;
        SF_time_rates(:,count_temp)=SF_rates;
        SF_time_weights(:,count_temp)=SF_weights;

    end

    fi=r_max./(1+exp(-b*(h-teta_adapt_before)));

    r_ALL(:,count_temp)=r_ALL(:,count_temp-1)+(-r_ALL(:,count_temp-1)+r0+fi)*(dt/tau)+noise;

    teta_adapt=teta_adapt_before+(-teta_adapt_before+h0+D_teta_temp*(r_ALL(:,count_temp-1)-r0))*(dt/tau_teta);
    teta_adapt_ALL(:,count_temp)=teta_adapt;

    meantemp_i=[];

    if (count_temp-1)>meanWin
        meantemp_i=mean(r_ALL(:,(count_temp-(meanWin)):(count_temp-1)),2); %Corrected 1209
    else
        meantemp_i=mean(r_ALL(:,(count_temp-(count_temp-1)):(count_temp-1)),2); %Corrected 1209
    end

    r_i_sub=r_ALL(:,count_temp)-meantemp_i;

    teta_adapt_before=teta_adapt;

    r_j_sub=[];
    r_j_sub=r_i_sub';

    % Learning rule

    W=W_before+LR*(r_i_sub*r_j_sub)*(dt/tau_w)-beta*W_before*(dt/tau_w);

    flag_max=0;
    flag_min=0;
    flag_both=0;

    if SaturationWeights==1||NoNegativeWeights==1

        indexes_overMax=[];
        indexes_overMax=find(W>Wmax);

        indexes_underMin=[];
        indexes_underMin=find(W<Wmin);

        indexes_negativeValues=[];
        indexes_negativeValues=find(W<0);

        W_reshaped=[];
        W_reshaped=reshape(W,[size(W,1)*size(W,2),1]);

        if SaturationWeights==1
            W_reshaped(indexes_overMax)=Wmax;
            W_reshaped(indexes_underMin)=Wmin;
        end

        if NoNegativeWeights==1
            W_reshaped(indexes_negativeValues)=0;
        end

        W=[];
        W=reshape(W_reshaped,[N,N]);

    end

    if autoconnections==0

        W_temp=[];
        W_temp=W;
        W=[];

        W = W_temp - diag(diag(W_temp));

    end

    Tot_Synapsis_before=sum(W_before,2);

    W_before=W;

    % mean of the weight matrix

    mean_weights(count_temp)=mean(W,'all');
    std_weights(count_temp)=std(W,0,'all');

    % mean weight within the assembly

    summ=0;
    countersteps=0;
    collecting=[];
    rowtemp=[];
    columntemp=[];

    for k=1:length(neurons_assembly)
        rowtemp=neurons_assembly(k);
        for l=1:length(neurons_assembly)
            countersteps=countersteps+1;
            columntemp=neurons_assembly(l);
            summ=summ+W(rowtemp,columntemp);
            collecting(countersteps)=W(rowtemp,columntemp);
        end
    end

    mean_weights_assembly(count_temp)=summ/(length(neurons_assembly)*length(neurons_assembly));
    std_weights_assembly(count_temp)=std(collecting);

    % mean weight outside the assembly

    summ=0;
    countersteps=0;
    collecting=[];
    rowtemp=[];
    columntemp=[];

    for k=1:length(neurons_NOTassembly)
        rowtemp=neurons_NOTassembly(k);
        for l=1:length(neurons_NOTassembly)
            countersteps=countersteps+1;
            columntemp=neurons_NOTassembly(l);
            summ=summ+W(rowtemp,columntemp);
            collecting(countersteps)=W(rowtemp,columntemp);
        end
    end

    mean_weights_NOTassembly(count_temp)=summ/(length(neurons_NOTassembly)*length(neurons_NOTassembly));
    std_weights_NOTassembly(count_temp)=std(collecting);

    % Mean connection FROM assembly neurons TO each neuron outside the
    % assembly (note that connections are symmetric)

    for n_neurons_outside=1:length(neurons_NOTassembly)
        summ=0;
        for i=1:length(neurons_assembly)
            summ=summ+W(neurons_assembly(i),neurons_NOTassembly(n_neurons_outside));
        end
        mean_weight_outsideTOinside(n_neurons_outside,count_temp)=summ/length(neurons_assembly);
    end

    % saving weight matrix

    flag_stim_FirstOnset_1=0; %weight matrices will be saved at each stim onset
    if ~isempty(find((t_temp-tonsets)<dt&(t_temp-tonsets)>=0))
        flag_stim_FirstOnset_1=1;
    end

    flag_stim_AfterStim_1=0; %weight matrices will be saved at the defined t_AfterStim
    if ~isempty(find((t_temp-t_AfterStim)<dt&(t_temp-t_AfterStim)>0))
        flag_stim_AfterStim_1=1;
    end

    flag_Before_Stimulation=0;
    if count_temp>45000/dt&&count_temp<tonset_first/dt %the weight matrix is saved starting from 45000 a.u. (this is to save it for some time before the stimulation phase starts, at 50000 a.u.)
        if countForSavingWeightMatrix>=2000
            flag_Before_Stimulation=1;
            countForSavingWeightMatrix=0;
        end
    end

    

    if flag_Before_Stimulation==1||flag_stim_FirstOnset_1==1||flag_stim_AfterStim_1==1


            clear WeightsMatrix_TP
            clear WeightsMatrix_TP_temp

            % this part of the script is meant to arrange the weight matrix
            % such that we save it
            % with the connections of the stimulated neurons in the
            % top-left of the matrix

            WeightsMatrix_TP=[];
            WeightsMatrix_TP=W(:,:);

            Ordered_WM=nan(N,N);

            for i = 1 : length(neurons_assembly)
                for j =1 : length(neurons_assembly)
                    Ordered_WM(i,j)=WeightsMatrix_TP(neurons_assembly(i),neurons_assembly(j));
                end
            end

            offset=[];
            offset=(length(neurons_assembly));
            counti=0;
            countj=0;

            for i = (1+offset) : (offset+length(neurons_NOTassembly))
                counti=counti+1;
                for j =(1+offset) : (offset+length(neurons_NOTassembly))
                    countj=countj+1;
                    Ordered_WM(i,j)=WeightsMatrix_TP(neurons_NOTassembly(counti),neurons_NOTassembly(countj));
                end
                countj=0;
            end

            offset=[];
            offset=(length(neurons_assembly));
            counti=0;
            countj=0;

            for i = 1 : length(neurons_assembly)
                counti=counti+1;
                for j =(1+offset) : (offset+length(neurons_NOTassembly))
                    countj=countj+1;
                    Ordered_WM(i,j)=WeightsMatrix_TP(neurons_assembly(counti),neurons_NOTassembly(countj));
                end
                countj=0;
            end

            offset=[];
            offset=(length(neurons_assembly));
            counti=0;
            countj=0;

            for i = (1+offset) : (offset+length(neurons_NOTassembly))
                counti=counti+1;
                for j =1 : length(neurons_assembly)
                    countj=countj+1;
                    Ordered_WM(i,j)=WeightsMatrix_TP(neurons_NOTassembly(counti),neurons_assembly(countj));
                end
                countj=0;
            end

            for i=1:N
                for j=1:N
                    if j==i
                        Ordered_WM(i,j)=0;
                    end
                end
            end

            countMatSaved=countMatSaved+1;

            cd(dirmain)

            Ordered_WM_forHist=[];
            Ordered_WM_forHist=reshape(Ordered_WM, [N*N 1]);

            eval(['save Mat' num2str(countMatSaved) '_counttemp' num2str(count_temp) ' Ordered_WM'])

            WM_saved(:,:,count_WM)=Ordered_WM;
            count_WM=count_WM+1;


        countForSavingWeightMatrix=0;

    end

    % update current time

    t_temp=t_temp+dt;
    count_temp=count_temp+1;

end


%% Saving the variables to be possibly checked for further analysis

% Display the firing rate map (with stimulated neurons at the top)

firingRates_neuronsAssembly=[];
firingRates_neuronsAssembly=r_ALL(neurons_assembly,:);

firingRates_neuronsNOTAssembly=[];
firingRates_neuronsNOTAssembly=r_ALL(neurons_NOTassembly,:);

firingRates_ordered=[];
firingRates_ordered=[firingRates_neuronsAssembly;firingRates_neuronsNOTAssembly];


figure
imagesc(firingRates_ordered)
colormap(jet)
colorbar

if SavingResults

    cd(dirmain)

    %Variables are splitted into 2 parts in order to be saved because of memory
    %reasons

    r_ALL_part1=[];
    r_ALL_part1=r_ALL(:,1:2500000);

    save AllRates_1 r_ALL_part1
    clear r_ALL_part1

    r_ALL_part2=[];
    r_ALL_part2=r_ALL(:,2500001:size(r_ALL,2));

    save AllRates_2 r_ALL_part2
    clear r_ALL_part2


    teta_adapt_ALL_part1=[];
    teta_adapt_ALL_part1=teta_adapt_ALL(:,1:2500000);
    save AllTeta_1 teta_adapt_ALL_part1
    clear teta_adapt_ALL_part1

    teta_adapt_ALL_part2=[];
    teta_adapt_ALL_part2=teta_adapt_ALL(:,2500001:size(teta_adapt_ALL,2));

    save AllTeta_2 teta_adapt_ALL_part2
    clear teta_adapt_ALL_part2


    SF_time_rates_part1=[];
    SF_time_rates_part1=SF_time_rates(:,1:2500000);
    SF_time_rates_part2=[];
    SF_time_rates_part2=SF_time_rates(:,2500001:size(SF_time_rates,2));

    save SF_time_rates_1 SF_time_rates_part1
    save SF_time_rates_2 SF_time_rates_part2

    clear SF_time_rates_part1
    clear SF_time_rates_part2


    SF_time_weights_part1=[];
    SF_time_weights_part1=SF_time_weights(:,1:2500000);

    save SF_time_weights_1 SF_time_weights_part1
    clear SF_time_weights_part1

    SF_time_weights_part2=[];
    SF_time_weights_part2=SF_time_weights(:,2500001:size(SF_time_weights,2));

    save SF_time_weights_2 SF_time_weights_part2
    clear SF_time_weights_part2


    mean_weight_outsideTOinside_1_part1=[];
    mean_weight_outsideTOinside_1_part1=mean_weight_outsideTOinside(:,1:2500000);

    save mean_weight_outsideTOinside_1_part1 mean_weight_outsideTOinside_1_part1
    clear mean_weight_outsideTOinside_1_part1

    mean_weight_outsideTOinside_1_part2=[];
    mean_weight_outsideTOinside_1_part2=mean_weight_outsideTOinside(:,2500001:size(mean_weight_outsideTOinside,2));

    save mean_weight_outsideTOinside_1_part2 mean_weight_outsideTOinside_1_part2
    clear mean_weight_outsideTOinside_1_part2


    param=[];

    param.N=N;
    param.tmax=tmax;

    param.beta=beta;
    param.LR=LR;

    param.dt=dt;
    param.tau=tau;

    param.HMS=HowManyStims;
    param.toffset_first=toffset_first;
    param.tonset_first=tonset_first;
    param.toffsets=toffsets;
    param.tonsets=tonsets;
    param.tsep=tsep;

    param.autoconnections=autoconnections;
    param.neurons_assembly=neurons_assembly;
    param.neurons_NOTassembly=neurons_NOTassembly;

    param.NoiseFactorRATE=Noise_Factor_RATE;
    param.NoNegativeWeights=NoNegativeWeights;
    param.SaturationWeights=SaturationWeights;
    param.Wmax=Wmax;
    param.Wmin=Wmin;
    param.r_max=r_max;

    param.h0=h0;
    param.b=b;

    param.count_temp_first=count_temp_first;
    param.t_temp_first=t_temp_first;

    param.D=D_teta;
    param.tau_teta=tau_teta;
    param.teta_adapt=teta_adapt;

    param.meanWin=meanWin;

    param.mean_weights_assembly=mean_weights_assembly;
    param.std_weights_assembly=std_weights_assembly;

    param.mean_weights_NOTassembly=mean_weights_NOTassembly;
    param.std_weights_NOTassembly=std_weights_NOTassembly;

    param.mean_weights=mean_weights;
    param.std_weights=std_weights;


    param.Factor_SR=Factor_SR;
    param.Factor_SW=Factor_SW;
    param.SR_thr=SR_thr;
    param.SW_thr=SW_thr;



    param.HowManyAfterStim=HowManyAfterStim;
    param.IntervalTest_AfterStim=IntervalTest_AfterStim; %in a.u.

    save param param


    Inputparam_WithI_WithSc_part1=[];
    Inputparam_WithI_WithSc_part1=Input_saved_WithI_WithSc(:,1:2500000);

    save Inputparam_WithI_WithSc_1 Inputparam_WithI_WithSc_part1
    clear Inputparam_WithI_WithSc_part1

    Inputparam_WithI_WithSc_part2=[];
    Inputparam_WithI_WithSc_part2=Input_saved_WithI_WithSc(:,2500001:size(Input_saved_WithI_WithSc,2));

    save Inputparam_WithI_WithSc_2 Inputparam_WithI_WithSc_part2
    clear Inputparam_WithI_WithSc_part2



end