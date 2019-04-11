function point_sim = comb_gen(Nseed,Nsim)
point_sim = [];
%% iterate till point_sim size reaches Nsim
while size(point_sim,1) < Nsim
    point_sim_tmp = randi(Nseed,Nsim,3);
    cont_list1    = 1:Nseed;
    index         = [];
    for cont1 = cont_list1
        tmp1               = perms(point_sim_tmp(cont1,:));
        cont_list11        = cont_list1;
        cont_list11(cont1) = [];
        for cont11 = cont_list11
            tmp2 = unique([tmp1; point_sim_tmp(cont11,:)],'rows');
            if size(tmp2,1) < 7
                index = [index; cont1];
            end
        end
    end
    index                  = unique(index,'rows'); 
    point_sim_tmp(index,:) = []; 
    point_sim              = [point_sim; point_sim_tmp];
end
point_sim((Nsim+1):end,:)  = [];
end
