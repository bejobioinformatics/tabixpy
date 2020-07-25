# Tests

## genVcfGzPy

### Pypy

`pypy3 -c 'import tabixpy;    tabixpy.genVcfGzPy("tests/annotated_tomato_150.vcf.bgz")'`
real    23m22,428s
user    22m57,712s
sys     00m23,493s

### Python

`python3 -c 'import tabixpy;    tabixpy.genVcfGzPy("tests/annotated_tomato_150.vcf.bgz")'`
real    19m54,345s
user    19m03,969s
sys     00m49,606s

## loadVcfGzPy

### Pypy

`pypy3 -c 'import tabixpy; _= tabixpy.loadVcfGzPy("tests/annotated_tomato_150.vcf.bgz")'`
real    00m1,843s
user    00m1,609s
sys     00m0,145s

### Python

`python3 -c 'import tabixpy; _= tabixpy.loadVcfGzPy("tests/annotated_tomato_150.vcf.bgz")'`
real    00m2,969s
user    00m2,813s
sys     00m0,096s


```sh
$ pypy3 -c 'import tabixpy;    tabixpy.genVcfGzPy("tests/annotated_tomato_150.vcf.bgz")'
INFO    : reading tests/annotated_tomato_150.vcf.bgz
INFO    : getAllPositions :: new chrom :: SL2.50ch00
INFO    : getAllPositions :: reading block        1,000        1,000
INFO    :   lastReal    7,183,757 first_pos    1,424,049 last_pos    1,425,379 block_len 6,930 block_size 65,280 num_cols    93 num_rows   108
INFO    : getAllPositions :: reading block       12,000       12,000
INFO    :   lastReal   95,089,867 first_pos   19,651,978 last_pos   19,652,441 block_len 6,425 block_size 65,280 num_cols    93 num_rows    71
INFO    : getAllPositions :: read    block       12,950       12,950 SL2.50ch00
INFO    : getAllPositions :: new chrom :: SL2.50ch01
INFO    : getAllPositions :: reading block       13,000           50
INFO    :   lastReal  103,387,646 first_pos       45,478 last_pos       46,117 block_len 6,483 block_size 65,280 num_cols    93 num_rows    65
INFO    : getAllPositions :: reading block      113,000      100,050
INFO    :   lastReal  858,118,353 first_pos   98,119,893 last_pos   98,120,374 block_len 6,647 block_size 65,280 num_cols    93 num_rows    56
INFO    : getAllPositions :: read    block      113,410      100,460 SL2.50ch01
INFO    : getAllPositions :: new chrom :: SL2.50ch02
INFO    : getAllPositions :: reading block      114,000          590
INFO    :   lastReal  865,561,494 first_pos      719,446 last_pos      720,521 block_len 6,462 block_size 65,280 num_cols    93 num_rows   109
INFO    : getAllPositions :: reading block      168,000       54,590
INFO    :   lastReal 1,269,244,987 first_pos   54,535,544 last_pos   54,536,225 block_len 7,567 block_size 65,280 num_cols    93 num_rows    65
INFO    : getAllPositions :: read    block      168,895       55,485 SL2.50ch02
INFO    : getAllPositions :: new chrom :: SL2.50ch03
INFO    : getAllPositions :: reading block      169,000          105
INFO    :   lastReal 1,276,532,263 first_pos       90,190 last_pos       91,118 block_len 7,518 block_size 65,280 num_cols    93 num_rows    62
INFO    : getAllPositions :: reading block      746,000       52,346
INFO    :   lastReal 5,660,159,559 first_pos   55,588,654 last_pos   55,589,816 block_len 6,981 block_size 65,280 num_cols    93 num_rows   103
INFO    : getAllPositions :: read    block      746,961       53,307 SL2.50ch11
INFO    : getAllPositions :: new chrom :: SL2.50ch12
INFO    : getAllPositions :: reading block      747,000           39
INFO    :   lastReal 5,667,499,242 first_pos       28,899 last_pos       30,035 block_len 7,681 block_size 65,280 num_cols    93 num_rows    62
INFO    : getAllPositions :: reading block      811,000       64,039
INFO    :   lastReal 6,153,516,785 first_pos   66,519,014 last_pos   66,519,382 block_len 7,078 block_size 65,280 num_cols    93 num_rows    57
INFO    : getAllPositions :: read    block      811,925       64,964 SL2.50ch12
INFO    :  saving tests/annotated_tomato_150.vcf.bgz.tbk
INFO    :  compressing
INFO    :  writing realPositions    -            811,411
INFO    :    fmt H min               0 max          13,027 len             12,949 cdsum       655,029,062,828 fmts <qc12949H
INFO    :    fmt L min           1,617 max     103,024,580 len            100,414 cdsum    48,944,202,966,369 fmts <qc100414L
INFO    :    fmt L min           4,108 max     861,164,995 len             55,470 cdsum    59,393,460,097,143 fmts <qc55470L
INFO    :    fmt L min           2,327 max   1,275,757,047 len             68,677 cdsum   105,198,202,501,533 fmts <qc68677L
INFO    :    fmt L min           3,449 max   1,785,635,506 len             71,011 cdsum   146,348,990,545,680 fmts <qc71011L
INFO    :    fmt L min           2,453 max   2,332,343,627 len             67,193 cdsum   174,044,363,385,931 fmts <qc67193L
INFO    :    fmt L min           2,141 max   2,849,186,456 len             50,764 cdsum   154,342,153,054,861 fmts <qc50764L
INFO    :    fmt L min           1,982 max   3,228,443,696 len             67,818 cdsum   236,620,999,143,857 fmts <qc67818L
INFO    :    fmt L min           1,577 max   3,746,913,672 len             64,038 cdsum   255,593,153,410,401 fmts <qc64038L
INFO    :    fmt L min           2,759 max   4,232,755,670 len             69,441 cdsum   312,502,266,251,081 fmts <qc69441L
INFO    :    fmt Q min           1,632 max   4,766,712,643 len             65,385 cdsum   327,871,852,359,839 fmts <qc65385Q
INFO    :    fmt Q min           3,702 max   5,261,129,562 len             53,303 cdsum   291,277,925,166,372 fmts <qc53303Q
INFO    :    fmt Q min           3,382 max   5,667,228,573 len             64,948 cdsum   384,105,489,815,498 fmts <qc64948Q
INFO    :  writing firstPositions   -            811,411
INFO    :    fmt H min             216 max          34,944 len             12,949 cdsum       127,248,454,214 fmts <qc12949H
INFO    :    fmt L min              25 max       2,281,164 len            100,414 cdsum     5,292,293,795,478 fmts <qc100414L
INFO    :    fmt L min             169 max       3,052,830 len             55,470 cdsum     1,782,315,375,611 fmts <qc55470L
INFO    :    fmt L min             101 max       2,582,078 len             68,677 cdsum     2,613,359,183,140 fmts <qc68677L
INFO    :    fmt L min              48 max       1,526,187 len             71,011 cdsum     2,395,545,129,807 fmts <qc71011L
INFO    :    fmt L min              17 max         856,885 len             67,193 cdsum     2,218,276,400,005 fmts <qc67193L
INFO    :    fmt L min              31 max       2,345,921 len             50,764 cdsum     1,408,571,368,855 fmts <qc50764L
INFO    :    fmt L min              10 max       2,294,715 len             67,818 cdsum     2,360,263,307,701 fmts <qc67818L
INFO    :    fmt L min              12 max         974,596 len             64,038 cdsum     2,227,577,918,367 fmts <qc64038L
INFO    :    fmt L min              71 max       1,472,425 len             69,441 cdsum     2,615,720,583,909 fmts <qc69441L
INFO    :    fmt L min               2 max         398,540 len             65,385 cdsum     2,200,671,373,829 fmts <qc65385L
INFO    :    fmt L min             154 max       1,312,213 len             53,303 cdsum     1,466,218,128,825 fmts <qc53303L
INFO    :    fmt L min             106 max         557,070 len             64,948 cdsum     2,235,800,132,812 fmts <qc64948L
INFO    :  writing lastPositions    -            811,411
INFO    :    fmt H min             210 max          34,935 len             12,949 cdsum       127,269,951,038 fmts <qc12949H
INFO    :    fmt L min              11 max       2,281,169 len            100,414 cdsum     5,292,391,271,726 fmts <qc100414L
INFO    :    fmt L min             170 max       3,052,815 len             55,470 cdsum     1,782,370,110,774 fmts <qc55470L
INFO    :    fmt L min             123 max       2,583,877 len             68,677 cdsum     2,613,426,599,885 fmts <qc68677L
INFO    :    fmt L min              43 max       1,526,172 len             71,011 cdsum     2,395,610,867,908 fmts <qc71011L
INFO    :    fmt L min              16 max         856,893 len             67,193 cdsum     2,218,341,531,065 fmts <qc67193L
INFO    :    fmt L min              27 max       2,345,938 len             50,764 cdsum     1,408,620,563,242 fmts <qc50764L
INFO    :    fmt L min               5 max       2,294,723 len             67,818 cdsum     2,360,330,583,947 fmts <qc67818L
INFO    :    fmt L min              25 max         974,609 len             64,038 cdsum     2,227,643,103,031 fmts <qc64038L
INFO    :    fmt L min              74 max       1,472,432 len             69,441 cdsum     2,615,792,306,282 fmts <qc69441L
INFO    :    fmt L min               3 max         398,532 len             65,385 cdsum     2,200,736,145,629 fmts <qc65385L
INFO    :    fmt L min             159 max       1,312,214 len             53,303 cdsum     1,466,273,809,167 fmts <qc53303L
INFO    :    fmt L min             105 max         557,055 len             64,948 cdsum     2,235,866,545,562 fmts <qc64948L
INFO    :  writing numberRows       -            811,411
INFO    :    fmt b min            -107 max              78 len             12,949 cdsum             1,205,613 fmts <qc12949b
INFO    :    fmt b min             -81 max              90 len            100,414 cdsum             8,870,687 fmts <qc100414b
INFO    :    fmt b min             -68 max             112 len             55,470 cdsum             4,647,715 fmts <qc55470b
INFO    :    fmt b min            -104 max             101 len             68,677 cdsum             6,010,106 fmts <qc68677b
INFO    :    fmt b min             -62 max              67 len             71,011 cdsum             6,289,795 fmts <qc71011b
INFO    :    fmt b min             -96 max              86 len             67,193 cdsum             6,066,674 fmts <qc67193b
INFO    :    fmt b min             -59 max              82 len             50,764 cdsum             4,322,527 fmts <qc50764b
INFO    :    fmt b min             -72 max              94 len             67,818 cdsum             6,112,010 fmts <qc67818b
INFO    :    fmt b min             -77 max              85 len             64,038 cdsum             5,779,794 fmts <qc64038b
INFO    :    fmt b min             -86 max              92 len             69,441 cdsum             6,298,745 fmts <qc69441b
INFO    :    fmt b min             -85 max              97 len             65,385 cdsum             5,978,975 fmts <qc65385b
INFO    :    fmt b min             -72 max             109 len             53,303 cdsum             4,811,114 fmts <qc53303b
INFO    :    fmt b min             -60 max              73 len             64,948 cdsum             5,955,424 fmts <qc64948b
INFO    : 3ab54eb34ba4f1848c9904cd336b761a3b2b1ac947195adc10e669ab717bc325
```


```sh
$ pypy3 -c 'import tabixpy; _= tabixpy.loadVcfGzPy("tests/annotated_tomato_150.vcf.bgz")'
INFO    :  loading tests/annotated_tomato_150.vcf.bgz.tbk
INFO    :  decompressing
INFO    :  reading realPositions
INFO    :         12,949 values
INFO    :        100,414 values
INFO    :         55,470 values
INFO    :         68,677 values
INFO    :         71,011 values
INFO    :         67,193 values
INFO    :         50,764 values
INFO    :         67,818 values
INFO    :         64,038 values
INFO    :         69,441 values
INFO    :         65,385 values
INFO    :         53,303 values
INFO    :         64,948 values
INFO    :  reading firstPositions
INFO    :         12,949 values
INFO    :        100,414 values
INFO    :         55,470 values
INFO    :         68,677 values
INFO    :         71,011 values
INFO    :         67,193 values
INFO    :         50,764 values
INFO    :         67,818 values
INFO    :         64,038 values
INFO    :         69,441 values
INFO    :         65,385 values
INFO    :         53,303 values
INFO    :         64,948 values
INFO    :  reading lastPositions
INFO    :         12,949 values
INFO    :        100,414 values
INFO    :         55,470 values
INFO    :         68,677 values
INFO    :         71,011 values
INFO    :         67,193 values
INFO    :         50,764 values
INFO    :         67,818 values
INFO    :         64,038 values
INFO    :         69,441 values
INFO    :         65,385 values
INFO    :         53,303 values
INFO    :         64,948 values
INFO    :  reading numberRows
INFO    :         12,949 values
INFO    :        100,414 values
INFO    :         55,470 values
INFO    :         68,677 values
INFO    :         71,011 values
INFO    :         67,193 values
INFO    :         50,764 values
INFO    :         67,818 values
INFO    :         64,038 values
INFO    :         69,441 values
INFO    :         65,385 values
INFO    :         53,303 values
INFO    :         64,948 values
INFO    : digestHex  3ab54eb34ba4f1848c9904cd336b761a3b2b1ac947195adc10e669ab717bc325
```