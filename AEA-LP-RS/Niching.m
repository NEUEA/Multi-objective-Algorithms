function Choose = Niching(Population,Choose,Front_k,N,back_Pop,fit)
  %% Association operation in the algorithm
    PopObj = Population.objs;
       
  %% Angle between each two solutions
    angle = acos(1-pdist2(PopObj,PopObj,'cosine'));
    
  %% Niching
     while sum(Choose) < N
        if ~isempty(Front_k)
            % Step 1: Diversity-based selection
            [~,rho] = max(min(angle(Front_k,Choose==1),[],2));
            Choose(Front_k(rho)) = 1;
%             x_k = Front_k(rho);
            Front_k(rho) = [];
             % Step 2: Convergence-based Replacement
             if ~isempty(Front_k) && rand()<0.5
                 Select = find(Choose==1);
                 [~,mu]      = min(min(angle(Front_k,Select),[],2));
                 [theta,r]   = min(angle(Front_k(mu),Select));
                 if theta < (pi/2/(N+1)) && fit(Select(r)) > fit(Front_k(mu))
                     Choose(Select(r))   = 0;
                     Choose(Front_k(mu)) = 1;
                     Front_k(mu) = [];
                 end
             end  
             % step 3: Information Renewing
% %              for j = 1:length(Front_k)
% %                  angle_new = acos(1-pdist2(Population(x_k).objs,Population(Frong_k(j)).objs,'cosine'));
% %                  angle = min(angle_new,angle(x_k,Front_k(j)));
% %              end
        else % Front_k end of use
            if ~isempty(back_Pop{1,1})
                if Choose(back_Pop{1,1}(1)) == 0
                Choose(back_Pop{1,1}(1)) = 1;
                end
                back_Pop{1,1}(1) = [];
            else %back_Pop end of use
                Remain = find(Choose == 0);
                [~,cho] = min(fit(Remain));
                Choose(Remain(cho)) = 1;
            end
        end
     end
    
end

