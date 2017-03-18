var tsvcolfile2 = './testcol3.tsv';
d3.tsv(tsvcolfile2, function(tsvcolumns) {
  var tsvfile = './GM12878_Rao_subset.plot_table.out';
  interactive_multiline_plot(tsvfile, tsvcolumns, 3, 7, -13, -8, 'log10 distance', 'log10prob', 'd3div_s5_GM12878_Rao_subset___');
});
