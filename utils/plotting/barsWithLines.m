function barsWithLines(these_labels,these_data,these_cols)

numBars = length(these_labels);

hold on
h=bar(1:numBars,diag(nanmean(these_data,1)),'stacked');
for i=1:numBars
    h(i).FaceColor=these_cols{i};
end
for j=1:size(these_data,1)
    plot([1:numBars],these_data(j,:),'o-k')
end

xlim([0,numBars+1])
xticks(1:numBars)
xticklabels(these_labels)

end