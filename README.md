# PathwayDenester
In current pathway enrichment methods each pathway is tested independently, without taking into consideration that these pathways are structured in hierarchical terms where the largest terms (composed of thousands of genes) contain smaller terms that partially intersect each other and in succession contain even smaller terms. This Nested pattern continues until the smallest terms that may have less than a dozen genes. Because of that some highly significant pathways will inevitably share all, or some, of their genes with other pathways that would be also considered enriched just based on the genes in its intersection.
This hitchhiker affects all terms that share a subset of their genes, and inflates their significance, often giving them better statistical significance than other relevant terms, therefore misleading the interpretation of the results.
This effect also produces results that are very repetitive with many pathways that can have more than 90% of their genes in common.

We created PathwayDenester an algorithm that takes the other pathways into consideration and assigns pathways a p-value based on the likelihood that it was enriched by its own and not by the genes in the intersection with more significant pathways.
This is achieved by testing each pathway (A) with each of the other pathways with higher significance (B), it takes into account the size of pathway A, its number of selected genes, the number of genes it has in common with pathway B, and the number of selected genes that fall in this intersection.

## Paper to be published.

## How to run:
Created and tested using version Python 3.11.5.

Reqquires python to run, can be run using command line arguments:
```
usage: PathwayDenester.py [-h] [--output_address [output_address]] [--to_test_threshold TO_TEST_THRESHOLD] [--pval_treshold PVAL_TRESHOLD]
                                      [--tranlator_gene_names [TRANLATOR_GENE_NAMES]]
                                      [paths_address] [gmt_address]

Look for pathways that are found only because of another more significant one:

positional arguments:
  paths_address         TSV result from pathways analysis like gprofiler, just columns "term_id" - Pathway ID as in the GMT file; "term_name" - Pathway Name, aesthetic only.;
                        "intersection_size" (optional), "term_size" (optional), and "p_value" or "adjusted_p_value" - From pathway enrichment, used for sorting pathways.; "intersection":
                        List of genes that are differentially expressed in pathway, separated by ",", must be quoted in case of CSV.
  gmt_address           GMT file with all pathways, will be used to find which genes are in each enriched pathway

options:
  -h, --help            show this help message and exit
  --output_address [output_address]
                        TSV file name where results will be saved. Default is pathway_list name + '_filtered.tsv'
  --to_test_threshold TO_TEST_THRESHOLD
                        Pathway 'B' will be tested against all pathways 'A' that are more significant AND that the ratio of B-DEGs also found in A to B-DEGs is at least
                        [to_test_threshold].
  --pval_treshold PVAL_TRESHOLD
                        P-value thresold to exclude a pathway. Since each pathway is treated independently, multiple testing corrections shouldn't be applied.
  --tranlator_gene_names [TRANLATOR_GENE_NAMES]
                        If you add a file that can translate the IDs used in gmt files to another name I can translate the genes list of each pathway. Argument are 3 comma separated
                        strings; translator file address, colname of IDs, colname of translated names
```

Example using results from gProfiler:
To run PathwayDenester you need a pathway enrichment analysis result and the GMT file that created it.
In <https://biit.cs.ut.ee/gprofiler/gost> you can click in "random example", then "Run query", click on the "Detail results" tab, and download the csv file.
The GMT file can be downloaded from the "Data sources" dropdown menu (combined ENSG.gmt).

Then run on terminal:

`python folder_where_pathwaydenester_is_saved/PathwayDenester.py folder_where_you_saved_csv_file/results.csv not_necessarily_the_same_folder/gprofiler_full_hsapiens.ENSG.gmt`
