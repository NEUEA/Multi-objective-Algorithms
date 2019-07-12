function [Population,omega] = TwoStageSelection(Population,cos,W,Global,omega)
N = Global.N;
[~,Region] = max(cos,[],2);
 %% Calculate the distance to the origin of each solution
   PopObj = Population.objs.^2;
   fit = sum(PopObj,2);

%% Two stages sorting
num_SubPop = size(W,1);%Number of individuals per subpopulation 
Sorting    = cell(1,num_SubPop);%Index of each subpopulation from good to bad 
Sub_Pop    = cell(1,num_SubPop);%Index of each subpopulation in the population 
Nr         = zeros(1,num_SubPop);%Number of sub-populations participating in the second-stage sequencing 
back_Pop   = cell(1,1);%Save backup subpopulation
sum_Sort   = 0;
for ii = 1:num_SubPop
     Sub_Pop{1,ii} = find(Region == ii)';
end
while sum_Sort<N
    for i = 1:num_SubPop
        if ~isempty(Sub_Pop{1,i})%Subpopulation is not empty
            %first stage
            first_Sorting = [];
            [FrontNo,MaxFNo,PP] = NDSort2(Population(Sub_Pop{1,i}).objs,Population(Sub_Pop{1,i}).cons,length(Sub_Pop{1,i}));
            for j=1:MaxFNo
                Front_j = find(FrontNo == j);
                if length(Front_j)==1
                    first_Sorting = [first_Sorting,Front_j];
                else% The number of non-dominant elements in this layer is greater than one
                    unique_PP = unique(PP(Front_j));
                    for k = length(unique_PP):-1:1
                        index = find(PP(Front_j)==unique_PP(k));
                        if length(index) == 1%PP value is only one
                            first_Sorting = [first_Sorting,Front_j(index)];
                        else
                            normP =  fit(Sub_Pop{1,i}(Front_j(index)));
                            [~,sort_dis] = sort(normP);
                            first_Sorting = [first_Sorting,Front_j(index(sort_dis))];
                        end
                    end
                end
            end
            %Second stage
            Nr(i) = min(length(Sub_Pop{1,i}),omega);%Global.evaluated/Global.evaluation*floor(length(Sub_Pop{1,i})/2)
            Angle = cos(Sub_Pop{1,i}(first_Sorting(1:Nr(i))),i);
            [~,div_sort] = sort(Angle,'descend');
            Sorting{1,i} = [Sorting{1,i},Sub_Pop{1,i}(first_Sorting(div_sort))];
            sum_Sort = sum_Sort + Nr(i);
            Sub_Pop{1,i}(first_Sorting(div_sort)) = [];
        else %This subpopulation has no individuals
            Nr(i) = 0;
            [~,temp] = sort(cos(:,i),'descend');%Ranking all individuals and vector(i) angles from small to large
            cons_Pop = temp(1:5);%Five of the smaller angles were selected
            normP =  sqrt(fit(cons_Pop));
            [~,sort_Conver] = sort(normP);
            back_Pop{1,1} = [back_Pop{1,1},cons_Pop(sort_Conver(1))];
            %
            %
        end
    end
end



%% Selection
Choose = zeros(1,length(Population));
level = 1;
while sum(Choose) <= N
    Front_k = [];
    for i = 1:num_SubPop
        if level <= length(Sorting{1,i})
             Front_k = [Front_k,Sorting{1,i}(level)];
        end      
    end
    %Niching
    if sum(Choose)+length(Front_k)>N
        omega = level;
        Choose = Niching(Population,Choose,Front_k,N,back_Pop,fit);
        break;
    end
    if sum(Choose)+length(Front_k)==N
         omega = level-1;
         Choose(Front_k) = 1;
         break;
    end
        Choose(Front_k) = 1;
        level = level +1;  
end

Population = Population(Choose==1);
end

