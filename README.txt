The script "DynamicAttractorNetworkModel.m" contains the code used to implement the rate attractor model described in the paper: "A dynamic attractor network model of memory formation, reinforcement and forgetting" by Marta Boscaglia, Chiara Gastaldi, Wulfram Gerstner and Rodrigo Quian Quiroga.

The script is set to run the case of "single-memory network" with the formation and reinforcement of an assembly of 10 neurons in a network of N=100. The assembly is stimulated with frequency=1/60 a.u. for 7000 times, starting from 50000 a.u. to 420000 a.u. (then, in the final part of the simulation, there is no stimulation for 10000 a.u., so that the dynamics of the network after the stimulation phase, "forgetting", can be checked).


The parameters in the script are defined as in Table 1 of the manuscript:

HowManyStimNeurons=10; %number of stimulated neurons
N=100; %size of the network
r_max=1; %maximal firing rate
r0=0; %baseline firing rate
tau=1; %timescale determining the neuronal activation
h0=0.15; %base firing threshold in the absence of firing
b=100; %slope parameter of the sigmoidal transfer function
Wmax=3/HowManyStimNeurons; %maximal synaptic weight
Wmin=-0.5/HowManyStimNeurons; %minimal synaptic weight
Noise_Factor_RATE=0.2/HowManyStimNeurons; %constant defining the noise
C_noise=sqrt(dt/tau)*Noise_Factor_RATE; %coupling constant of the Gaussian noise. It contains a factor sqrt(dt/tau) as explained in the description of eq. 6 in the manuscript
meanWin=15/dt; %length of the window adopted for calculating the learning rule's running average (This is in a.u.*10)
Factor_SR=2; %divisive normalization constant
Factor_SW=1; %synaptic normalization constant
tau_teta=7*tau; %neural adaptation time constant
D_teta=1; %constant determining the adaptation strength
LR=1; %learning rate
beta=0.0025; %forgetting rate
tau_w=50; %timescale determining the learning

----------------------------------------------------

The script can be set to run the case of "single-memory network" with frequencies different than frequency=1/60 a.u. (as in Fig 6) by changing the definition of the number of the stimulations (i.e. HowManyStims) and the time distance between two consecutive stimulations (i.e. tsep) (see line 89 to line 96).
Specifically,
- if frequency=1/30 a.u. --> 'HowManyStims=14000;', 'tsep=30;'
- if frequency=1/40 a.u. --> 'HowManyStims=10500;', 'tsep=40;'
- if frequency=1/120 a.u. --> 'HowManyStims=3500;', 'tsep=120;'

The code can be adapted to simulate a network with 2 stimulated assemblies (as in Fig 7 and Fig 8) by defining a second vector indicating the neurons to be stimulated (see line 79 to line 85; define a vector 'neurons_assembly_2'; note that, consequently, there will also be a second vector of external stimulation, I2) and, then, defining the characteristics of the stimulation for the second assembly (e.g. in case of two assemblies stimulated with the same frequency, as in Fig 7, 'HowManyStims_2=HowManyStims;' 'tsep_2=tsep;' 'tonset_first_2=50030;' 'toffset_first_2=50035;').
Then, in the "Evolution of the network dynamics", you need to consider the second stimulation (with I2) to be given to 'neurons_assembly_2' (line 196 to line 252 can be taken as an example of how to give stimulation to an assembly, according to the parameters of the stimulation for that specific assembly).
If you still want to display, at multiple times, the weight matrix with the stimulated neurons at the top of the matrix, you will need to consider that there is a second assembly (see line 466 to line 548). The same consideration should be done in every case the first assembly is considered (this means adapting the plot of firing rates, saving of param.mat etc..)

The same reasoning applies for other possible adaptations of the code (e.g. a network with 3 patterns). Specifically, if you want to have two phases of stimulation to the same assemblies, as in Fig.9, you need to define the characteristics of the stimulation (i.e. HowManyStims, tsep, tonsets, toffsets etc..) for each assembly in each stimulation phase, and use them all in the "Evolution of the network dynamics" section. 











 