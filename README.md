# SeqRuler

<span align='center'>
  <img src="https://github.com/CDCgov/SeqRuler/assets/36117344/e966ff6b-27f4-4521-9e54-eb7d84f3f5e0">
</span>

<p align='center'>A simple GUI interface for calculating TN93 distance or SNP (Hamming) distance between all sequences in an input fasta file.</p>


<div align='center'>
  <img width=80% height=80% alt="Screenshot 2023-05-19 at 2 41 18 PM" src="https://github.com/CDCgov/SeqRuler/assets/36117344/f7afeaad-4b7f-43e8-9d32-992dbf1baebc">
</div>

## Download

- Download:
[SeqRuler.jar](https://github.com/CDCgov/SeqRuler/releases/download/v4.3/SeqRuler.jar)

## Or compile from source
- Compilation:
```bash
mvn clean install
```

## Running

```bash
java -jar SeqRuler.jar
```

## Help

```bash
Usage: SeqRuler [-egGhnprsSV] [-a=<ambiguityHandling>] [-c=<cores>]
                [-d=<distanceMethod>] [-f=<max_ambiguity_fraction>] [-i=FILE]
                [-o=FILE] [-t=<edgeThresholdString>]
  -a, --ambiguity, --ambiguities=<ambiguityHandling>
                             How to handle ambiguous nucleotides. One of [resolve,
                               average, gapmm, skip]
  -c, --cores=<cores>        Number of cores to use for parallel processing.
                               Default: 1
  -d, --distance-method=<distanceMethod>
                             distance metric to use. One of [TN93, SNP]. Default:
                               TN93
  -e, --enumerate_sequences  Enumerate sequences for output file, and produce
                               additional map file giving integer to sequence name
                               mapping. Default: false
  -f, --fraction=<max_ambiguity_fraction>
                             Maximum allowable fraction of ambiguities allowed for
                               'resolve' mode. If exceeded, use 'average' mode.
  -g, --ignore-terminal-gaps Ignore terminal gaps at beginning and end of sequences
                               when calculating distances. [SNP only] Default: true
  -G, --ignore-all-gaps      Ignore all gaps when calculating distances. [SNP only]
                               Default: false
  -h, --help                 Show this help message and exit.
  -i, --inFile=FILE          input file with sequences
  -n, --ignore-ambiguities   Ignore ambiguities when calculating distances. [SNP
                               only] Default: true
  -o, --outFile=FILE         output file with distances
  -p, --pairs                read pairs of sequences from stdin, calculate distance
                               for each pair. format "name1, seq1, name2, seq2\n"
  -r, --run-server           run jetty server
  -s, --stdin                read fasta from stdin. Alternative to reading from a
                               file (-i)
  -S, --stdout               write distances to stdout. Alternative to writing to a
                               file (-o)
  -t, --edge-threshold=<edgeThresholdString>
                             edges above the threshold are not reported in output.
                               {Default: 1.0 (TN93), inf (SNP)}
  -V, --version              Print version information and exit.
```
