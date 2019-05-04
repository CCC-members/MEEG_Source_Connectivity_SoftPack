function indms = best_ranked(J,areas,n_rank)
indms = cell(1,length(areas));
for cont = 1:length(areas)
    Jtmp        = J(areas{cont});
    rank_cont   = tiedrank(Jtmp);
    ind         = find(rank_cont > (length(Jtmp)-n_rank));
    indms{cont} = areas{cont}(ind);
end
end