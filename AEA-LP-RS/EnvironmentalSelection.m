function [Population,omega] = EnvironmentalSelection(Population,Global,W,omega)
% The environmental selection of AEA-LP-RS

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Normalize the objective values of the current Population
    fmin  = min(Population.objs,[],1);
    fmax  = max(Population.objs,[],1);
    PopObj = (Population.objs-repmat(fmin,size(Population,1),1))./repmat(fmax-fmin,size(Population,1),1);
    
    %% Subpopulation partition
    cos = 1-pdist2(PopObj,W,'cosine');   
    
    %% Two Stage Ranking and Selection     
    [Population,omega] = TwoStageSelection(Population,cos,W,Global,omega);
   
end


