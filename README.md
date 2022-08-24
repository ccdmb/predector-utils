# predector-utils
Utilities for running the predector pipeline.


This is where the models and more extensive utility scripts are developed.

All command line tools are accessible under the main command `predutils`.


## `predutils load_db`

Load a line-delimited JSON into an SQLite database.

Basic usage:

```bash
predutils load_db results.db results.ldjson
```

Options:

```
--replace-name
  Replace record names with the string 'd'.

--drop-null-dbversion
  Don't enter records requiring a database but have no database version specified.

--mem <float>
  Use this much RAM for the SQLite cache. Set this to at most half of the total ram available.

--include ANALYSIS [ANALYSIS ...]
  Only include these analyses, specify multiple with spaces.

--exclude ANALYSIS [ANALYSIS ...]
  Exclude these analyses, specify multiple with spaces. Overrides analyses specified in --include.
```

If you are using this to set up a pre-computed database, specify the `--replace-name`, `--drop-null-dbversion` flags which will make sure any duplicate entries are excluded.
It is also recommended that you exclude results that are faster to recompute than to enter and fetch from the database, including `pepstats`, `kex2_cutsite`, and `rxlr_like_motif`.

This will be relatively fast for small to medium datasets, but can take several hours for many millions of entries.
Setting the `--mem` option is also a good idea to speed up inserting larger datasets.


## `predutils r2js`

Convert the output of one of the analyses into a common [line delimited JSON](http://ndjson.org/) format.
The json records retain all information from the original output files, but are much easier to parse because each line is just JSON.

Basic usage:

```bash
predutils r2js \
  -o outfile.ldjson \
  --software-version 1.0 \
  --database-version 1.0 \
  --pipeline-version 0.0.1 \
  pfamscan \
  pfamscan_results.txt \
  in.fasta
```

Analyses available to parse in place of `pfamscan` are:
`signalp3_nn`, `signalp3_hmm`, `signalp4`, `signalp5`, `deepsig`, `phobius`, `tmhmm`,
`deeploc`, `targetp_plant` (v2), `targetp_non_plant` (v2), `effectorp1`, `effectorp2`,
`apoplastp`, `localizer`, `pfamscan`, `dbcan` (HMMER3 domtab output), `phibase` \*, `pepstats`,
`effectordb`, `kex2_cutsite`, and `rxlr_like_motif`.

\* assumes search with MMseqs with tab delimited output format columns: query, target, qstart, qend, qlen, tstart, tend, tlen, evalue, gapopen, pident, alnlen, raw, bits, cigar, mismatch, qcov, tcov.


## `predutils encode`

Preprocess some fasta files.

1. Strips trailing `*` amino acids from sequences, removes `-` characters, replaces internal `*`s and other redundant/non-standard amino acids with `X`, and converts sequences to uppercase.
2. removes duplicate sequence using a checksum, saving the mapping table to recover the duplicates at the end of the analysis.
3. Replace the names of the deduplicated sequences with a short simple one.


Basic usage:

```bash
predutils encode \
  output.fasta \
  output_mapping.tsv \
  input_fastas/*
```

Note that it can take multiple input fasta files, and the filename is saved alongside the sequences in the output mapping table to recover that information.


By default, the temporary names will be `SR[A-Z0-9]5` e.g. `SR003AB`.
You can change the prefix (default `SR`) with the `--prefix` flag, and the number of id characters (default 5) with the `--length` parameter.


## `predutils split_fasta`

Splits a fasta files into several files each with a maximum of n sequences.

Basic usage:

```bash
predutils split_fasta --template 'chunk{index}.fasta' --size 100 in.fasta
```

The `--template` parameter accepts python `.format` style string formatting, and
is provided the variables `fname` (the input filename) and `index` (the chunk number starting at 1).
To pad the numbers with zeros for visual ordering in directories, use the something like `--template '{fname}.{index:0>4}.fasta'`.
Directories in the template will be created for you if they don't exist.


## `predutils precomputed`

Takes a database and some sequences and uses the sequence checksums to decide what has already been computed.
Outputs the precomputed results if `-o` is set.
Writes fasta for the remaining sequences to be computed.

The analyses and software versions to check for in the database are specified as a tab separated file to `analyses`.

```
usage: predutils precomputed [-h] [-o OUTFILE] [-t TEMPLATE] [--mem MEM] db analyses infasta

positional arguments:
  db                    Where the sqlite database is
  analyses              A 3 column tsv file, no header. 'analysis<tab>software_version<tab>database_version'. database_version should be empty string if None.
  infasta               The fasta file to parse as input. Cannot be stdin.

optional arguments:
  -h, --help            show this help message and exit
  -o OUTFILE, --outfile OUTFILE
                        Where to write the precomputed ldjson results to.
  -t TEMPLATE, --template TEMPLATE
                        A template for the output filenames. Can use python `.format` style variable analysis. Directories will be created.
  --mem MEM             The amount of RAM in gibibytes to let SQLite use for cache.
```


## `predutils decode`

The other end of `predutils encode`.
Takes the common line delimited format from analyses and separates them back
out into the original filenames.

```bash
predutils decode [-h] [-t TEMPLATE] [--mem MEM] db map

positional arguments:
  db                    Where the sqlite database is
  map                   Where to save the id mapping file.

optional arguments:
  -h, --help            show this help message and exit
  -t TEMPLATE, --template TEMPLATE
                        What to name the output files.
  --mem MEM             The amount of RAM in gibibytes to let SQLite use for cache.

predutils decode \
  --template 'decoded/{filename}.ldjson' \
  results.db \
  output_mapping.tsv
```

We use the template flag to indicate what the filename output should be, using python format
style replacement. Available values to `--template` are `filename` and `filename_noext`.
The latter is just `filename` without the last extension.


## `predutils tables`

Take the common line delimited output from `predutils r2js` and recover a tabular version of the raw data.
Output filenames are controlled by the `--template` parameter, which uses python format style replacement.
Currently, `analysis` is the only value available to the template parameter.
Directories in the template will be created automatically.

```
predutils tables [-h] [-t TEMPLATE] [--mem MEM] db

positional arguments:
  db                    Where to store the database

optional arguments:
  -h, --help            show this help message and exit
  -t TEMPLATE, --template TEMPLATE
                        A template for the output filenames. Can use python `.format` style variable analysis. Directories will be created.
  --mem MEM             The amount of RAM in gibibytes to let SQLite use for cache.


predutils tables \
  --template "my_sample-{analysis}.tsv" \
  results.db
```


## `predutils gff`

Take the common line-delimited json output from `predutils r2js` and get a GFF3 formatted
set of results for analyses with a positional component (e.g. signal peptides, transmembrane domains, alignment results).

```
predutils gff \
  --outfile "my_sample.gff3" \
  results.ldjson
```

By default, mmseqs and HMMER search results will be filtered by the built in significance thresholds.
To include all matches in the output (and possibly filter by your own criterion) supply the flag `--keep-all`.


## `predutils rank`

Take a database of results entered by `load_db` and get a summary table
that includes all of the information commonly used for effector prediction, as well as
a scoring column to prioritise candidates.

```
predutils rank \
  --outfile my_samples-ranked.tsv \
  results.db
```


To change that Pfam or dbCAN domains that you consider to be predictive of effectors,
supply a text file with each pfam or dbcan entry on a new line (do not include pfam version number or `.hmm` in the ids) to the parameters `--dbcan` or `--pfam`.

You can also change the weights for how the score columns are calculated.
See `predutils rank --help` for a full list of parameters.


## `predutils regex`

```
predutils regex [-h] [-o OUTFILE] [-l] [-k {kex2_cutsite,rxlr_like_motif,custom}] [-r REGEX] INFILE

positional arguments:
  INFILE                The input fasta file.

optional arguments:
  -h, --help            show this help message and exit
  -o OUTFILE, --outfile OUTFILE
                        Where to write the output to. Default: stdout
  -l, --ldjson          Write the output as ldjson rather than tsv.
  -k {kex2_cutsite,rxlr_like_motif,custom}, --kind {kex2_cutsite,rxlr_like_motif,custom}
                        Which regular expressions to search for. If custom you must specify a regular expression to --regex. Default: custom.
  -r REGEX, --regex REGEX
                        The regular expression to search for. Ignored if --kind is not custom.
```


## `predutils dump_db`

```
predutils dump_db [-h] [-o OUTFILE] [--mem MEM] db

positional arguments:
  db                    The database to dump results from.

optional arguments:
  -h, --help            show this help message and exit
  -o OUTFILE, --outfile OUTFILE
                        Where to write the output to. Default: stdout
  --mem MEM             The amount of RAM in gibibytes to let SQLite use for cache.
```


## `predutils map_to_genome`

This script projects the protein GFF file from the results into genome coordinates based on the GFF used to extract the proteins.
This is intended to support visualisation and selection of candidates in genome browsers like JBrowse or Apollo.

```
usage: predutils map_to_genome [-h] [-o OUTFILE] [--split {source,type}] [--no-filter-kex2] [--id ID_FIELD] genes annotations

positional arguments:
  genes                 Gene GFF to use.
  annotations           Annotation GFF to use (i.e. the GFF3 output of predector)

options:
  -h, --help            show this help message and exit
  -o OUTFILE, --outfile OUTFILE
                        Where to write the output to. If using the --split parameter this
                        becomes the prefix. Default: stdout
  --split {source,type}
                        Output distinct GFFs for each analysis or type of feature.
  --no-filter-kex2      Output kex2 cutsites even if there is no signal peptide
  --id ID_FIELD         What GFF attribute field corresponds to your protein feature seqids?
                        Default uses the Parent field. Because some fields (like Parent)
                        can have multiple values, we'll raise an error if there is more
                        than 1 unique value. Any CDSs missing the specified
                        field (e.g. ID) will be skipped.
```

Pay attention to the `--id` parameter. This determines how the protein ID from the predector output will be matched to the genome GFF3.
By default it looks at the `Parent` attribute, as this is the default name from many GFF protein extraction tools.
But any attribute in the 9th column of the GFF can be used, e.g. to use the CDS feature `ID`, use `--id ID`. 
This script will only consider records with type `CDS`, and matches names exactly (no partial matches).

`--split` provides a convenience function to split the GFFs into multiple GFFs based on type or analysis.
This is to facilitate loading the different features/analyses as separate tracks in the genome browser.


## `predutils score_to_genome`

This script plots numeric scores from the predector ranked table at the CDS coordinates from the GFF used to extract the proteins.
This is intended to support visualisation and selection of candidates from genome browsers like JBrowse or Apollo.

```
usage: predutils scores_to_genome [-h] [-o OUTFILE] [--target TARGET [TARGET ...]] \
  [--id ID_FIELD] [--reducer {min,max,mean,median}] genes annotations

positional arguments:
  genes                 Gene GFF to use.
  annotations           ranking table to use (i.e. the -ranked.tsv output of predector)

options:
  -h, --help            show this help message and exit
  -o OUTFILE, --outfile OUTFILE
                        Where to write the output to. Default: stdout
  --target TARGET [TARGET ...]
                        Only output these columns into the bedgraph file.
                        By default writes a multi bedgraph with all columns.
  --id ID_FIELD         What GFF attribute field corresponds to your protein feature seqids?
                        Default uses the Parent field. Because some fields (like Parent)
                        can have multiple values, we'll raise an error if there is more
                        than 1 unique value. Any CDSs missing the specified
  --reducer {min,max,mean,median}
                        How should we combine scores if two features overlap?
```


As with the `map_to_genome` command, you should ensure that the protein names (in the first column of the `-ranked.tsv` table), match with the attribute in the 9th column of the GFF3 file.

Because multiple CDS features can overlap (e.g. if there are alternative splicing patterns), we must combine those scores somehow so that a single value is associated with the genomic region.
This can be specified with the `--reducer` parameter which takes the maximum value by default (e.g. the highest effectorP score). But if instead you wanted the average use `--reducer mean`.

By default this script outputs a multi-bedgraph file with all numeric columns.
But if you only wish to output a subset of those columns, you can provide them as a space separated list to the `--target` parameter.
Note that an error will be raised if you try to access a non-numeric column.
If the `--target` parameter comes directly before the two positional arguments, you need to indicate where the space separated list stops with a `--`.

e.g.


```
predutils scores_go_genome --id ID --target effector_score apoplastp effectorp2 -- mygenes.gff3 mygenes-ranked.tsv
```


To split the output into single column bedgraphs you can use the following to extract one bedgraph per value column.


```
#!/usr/bin/env bash

set -euo pipefail

# OPTIONAL
# BGTOBW_EXE="bedGraphToBigWig" # This comes from the UCSC toolkit, or here: http://hgdownload.soe.ucsc.edu/admin/exe/
INFILE=$1
FAI=$2
PREFIX=$3

COLUMNS=( $(head -n 1 "${INFILE}" | cut -f4-) )

for i in "${!COLUMNS[@]}"
do
    idx=$(( ${i} + 4 ))
    col="${COLUMNS[${i}]}"
    cut -f1,2,3,${idx} "${INFILE}" | tail -n+2 | awk '$4 != "NA"' > "${PREFIX}${col}.bedgraph"
    # OPTIONAL
    # ${BGTOBW_EXE} "${PREFIX}${col}.bedgraph" "${FAI}" "${PREFIX}${col}.bw"
done
```

This is useful for converting the tracks into bigwig files (which only support a single value), suitable for use with Jbrowse or Apollo.
