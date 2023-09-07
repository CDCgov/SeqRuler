# SeqRuler
<div align='center'>
  <img width=80% height=80% alt="Screenshot 2023-05-19 at 2 41 18 PM" src="https://github.com/CDCgov/SeqRuler/assets/36117344/f7afeaad-4b7f-43e8-9d32-992dbf1baebc">
</div>

GUI interface for calculating TN93 distance or SNP (Hamming) distance between all sequences in the input fasta file.

## Download or Compilation

- Download:
[SeqRuler.jar](https://github.com/CDCgov/SeqRuler/releases/download/v3.5/SeqRuler.jar)


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
java -jar SeqRuler.jar -h
                                    
Usage: SeqRuler [-gGhprsSV] [-a=<ambiguityHandling>] [-c=<cores>]
                [-d=<distanceMethod>] [-f=<max_ambiguity_fraction>] [-i=FILE]
                [-o=FILE] [-t=<edgeThresholdString>]
  -a, --ambiguity, --ambiguities=<ambiguityHandling>
                          How to handle ambiguous nucleotides. One of [resolve,
                            average, gapmm, skip]
  -c, --cores=<cores>     Number of cores to use for parallel processing.
  -d, --distance-method=<distanceMethod>
                          distance metric to use. One of [TN93, SNP]. Default: TN93
  -f, --fraction=<max_ambiguity_fraction>
                          Maximum allowable fraction of ambiguities allowed for
                            'resolve' mode. If exceeded, use 'average' mode.
  -g, --ignore-terminal-gaps
                          Ignore terminal gaps at beginning and end of sequences
                            when calculating distances. [SNP only] Default: true
  -G, --ignore-all-gaps   Ignore all gaps when calculating distances. [SNP only]
                            Default: false
  -h, --help              Show this help message and exit.
  -i, --inFile=FILE       input file with sequences
  -o, --outFile=FILE      output file with distances
  -p, --pairs             read pairs of sequences from stdin, calculate distance for
                            each pair. format "name1, seq1, name2, seq2\n"
  -r, --run-server        run jetty server
  -s, --stdin             read fasta from stdin. Alternative to reading from a file
                            (-i)
  -S, --stdout            write distances to stdout. Alternative to writing to a
                            file (-o)
  -t, --edge-threshold=<edgeThresholdString>
                          edges above the threshold are not reported in output
  -V, --version           Print version information and exit.
```
