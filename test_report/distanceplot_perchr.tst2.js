var tsvcolfile2 = './testcol2.tsv';
d3.tsv(tsvcolfile2, function(tsvcolumns) {
  var tsvfile = './tst2.plot_table.out';
  interactive_multiline_plot(tsvfile, tsvcolumns, 3, 7, -13, -9, 'log10 distance', 'log10prob', 'd3div_s5_tst2___');
});
