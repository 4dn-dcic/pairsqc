var tsvcolfile2 = './testcol3.tsv';
d3.tsv(tsvcolfile2, function(tsvcolumns) {
  var tsvfile = './K562_Rao.plot_table.out';
  interactive_multiline_plot(tsvfile, tsvcolumns, 3, 7, -9, -5, 'log10 distance', 'log10prob', 'd3div_s5_K562_Rao___');
});
