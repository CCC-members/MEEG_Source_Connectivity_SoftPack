function assign_labels(labels_conn,ngen,dx,dy,slide)
for k_label = 1:length(labels_conn)
    text(-dx,ngen/slide + (k_label - 1)*ngen,labels_conn{k_label},'HorizontalAlignment','center')
    text(ngen/slide + (k_label - 1)*ngen,length(labels_conn)*ngen + dy,labels_conn{k_label},'HorizontalAlignment','center')
end
end