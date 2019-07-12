function AEALPRS(Global)
% <algorithm> <A-G>
% AEA-LP-RS
% H --- 91 --- Reference vector initialization parameter

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    [H,omega] = Global.ParameterSet(91,100);

    %% Generate the weight vectors
    W = UniformPoint(H,Global.M);   
    
    %% Generate random population
    Population = Global.Initialization();
    
    %% Optimization
    while Global.NotTermination(Population)        
%        MatingPool = randi(Global.N,1,Global.N);
       MatingPool = TournamentSelection(2,Global.N,sum(max(0,Population.cons),2));
       Offspring  = Global.Variation(Population(MatingPool));
       [Population,omega] = EnvironmentalSelection([Population,Offspring],Global,W,omega);   
    end
end