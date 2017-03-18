var tsvfile = './test_report/tst.plot_table.out';
var tsvcolfile2 = './testcol2.tsv';
d3.tsv(tsvcolfile2, function(tsvcolumns) {
  interactive_multiline_plot(tsvfile, tsvcolumns, 3, 7, -13, -9, 'log10 distance', 'log10prob', 'd3div_hahahaha1');
});
d3.tsv(tsvcolfile2, function(tsvcolumns) {
  interactive_multiline_plot(tsvfile, tsvcolumns, 3, 7, -13, -9, 'log10 distance', 'log10prob', 'd3div_hahahaha2');
});
d3.tsv(tsvcolfile2, function(tsvcolumns) {
  interactive_multiline_plot(tsvfile, tsvcolumns, 3, 7, -13, -9, 'log10 distance', 'log10prob', 'd3div_hahahaha3');
});
