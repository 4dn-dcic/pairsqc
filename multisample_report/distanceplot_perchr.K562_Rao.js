d3.tsv('./tsvcol.d3div_contact_frequency_vs_genomic_separation_per_chr_K562_Rao___.tsv', function(tsvcolumns) {
interactive_multiline_plot('./K562_Rao.plot_table.out', tsvcolumns, 3.5,7,-8.969,-5.265, 'log10 genomic separation', 'log10 contact frequency', 'd3div_contact_frequency_vs_genomic_separation_per_chr_K562_Rao___');
});
