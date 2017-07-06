d3.tsv('./tsvcol.d3div_contact_frequency_vs_genomic_separation_per_chr_sample___.tsv', function(tsvcolumns) {
interactive_multiline_plot('./sample.plot_table.out', tsvcolumns, 3.5,7,-12.307,-8.371, 'log10 genomic separation', 'log10 contact frequency', 'd3div_contact_frequency_vs_genomic_separation_per_chr_sample___');
});
