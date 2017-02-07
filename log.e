Traceback (most recent call last):
  File "pairsqc", line 172, in <module>
    distance_histogram ( args.pairs, args.chrsize, max_logdistance = 8 , min_logdistance = 1 )
  File "pairsqc", line 109, in distance_histogram
    count[orientation][bin_number] += 1
KeyError: '04'
awk: !=4
awk: ^ syntax error
/n/data1/hms/transion/park/juicer-sample-dir/SRR1658650/aligned/merged_nodups.-MT.txt does not exist or does not contain any reads.
  File "pairsqc", line 77
    count = { a: [0] * (max_bin_number + 1) for a in orientation_list }
                                              ^
SyntaxError: invalid syntax
Traceback (most recent call last):
  File "pairsqc", line 172, in <module>
    distance_histogram ( args.pairs, args.chrsize, max_logdistance = 8 , min_logdistance = 1 )
  File "pairsqc", line 109, in distance_histogram
    count[orientation][bin_number] += 1
KeyError: '04'
