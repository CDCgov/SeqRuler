package SNP;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.FileWriter;
import java.io.BufferedReader;
import java.util.*;
import java.util.concurrent.atomic.AtomicReference;
//Parallelization
import java.util.concurrent.*;

public class SNP extends Observable {

    private static class Triplet<A, B, C> {
        final A first;
        final B second;
        final C third;
    
        public Triplet(A first, B second, C third) {
            this.first = first;
            this.second = second;
            this.third = third;
        }
    }

    private File inputFile;
    private File outputFile;
    private float edgeThreshold;
    private int cores = 1;
    private boolean ignoreAmbiguities = true;
    private boolean ignoreTerminalGaps = true;
    private boolean ignoreAllGaps = false;
    private boolean use_stdin = false;
    private boolean use_stdout = false;
    private boolean input_as_pairs = false;

    public void setInputAsPairs(boolean input_as_pairs) {
        this.input_as_pairs = input_as_pairs;
    }

    public void setUseStdin(boolean use_stdin) {
        this.use_stdin = use_stdin;
    }
    public void setUseStdout(boolean use_stdout) {
        this.use_stdout = use_stdout;
    }
    public void setInputFile(File inputFile) {
        this.inputFile = inputFile;
    }
    public void setOutputFile(File outputFile) {
        this.outputFile = outputFile;
    }
    public void setEdgeThreshold(float edgeThreshold) {
        this.edgeThreshold = edgeThreshold;
    }
    public void setCores(int cores) {
        this.cores = cores;
    }
    public void setIgnoreAmbiguities(boolean ignoreAmbiguities) {
        this.ignoreAmbiguities = ignoreAmbiguities;
    }
    public void setIgnoreTerminalGaps(boolean ignoreTerminalGaps) {
        this.ignoreTerminalGaps = ignoreTerminalGaps;
    }
    public void setIgnoreAllGaps(boolean ignoreAllGaps) {
        this.ignoreAllGaps = ignoreAllGaps;
    }


    public void snpFasta() {
        if (!use_stdout) System.out.println("Running SNPFASTA");
        PrintWriter f = null;
        try {
            if (!use_stdout) System.out.println("Reading input file...");
            ArrayList<Seq> seqs;
            if (use_stdin) {
                seqs = read_fasta_stdin();
            } else {
                seqs = read_fasta(inputFile);
            }
            if (!use_stdout) System.out.println("Calculating distances...");
            snp(seqs);
        }
        catch(FileNotFoundException e) {
            e.printStackTrace();
        }
        finally {
            if(f != null) f.close();
        }
    }

    public void snp(ArrayList<Seq> seqs) {
        if (this.cores == 1) {
            snp_sequential(seqs);
        } else {
            snp_parallel(seqs);
        }
    }

    public void snp_sequential(ArrayList<Seq> seqs) {
        //sequential version
        long pairs_count = ((long) seqs.size() * (seqs.size() - 1)) / 2;
        long current_pair = 0;
        long startTime = System.nanoTime(), estimatedTime;

        PrintWriter f = null;
        if (use_stdout)
            System.out.println("Source,Target,Distance");
        else {
            try {
                f = new PrintWriter(new BufferedWriter(new FileWriter(outputFile)));
            } catch (IOException e) {
                e.printStackTrace();
            }
            f.println("Source,Target,Distance");
        }
        for (int i = 1; i < seqs.size(); ++i) {
            if (!use_stdout) System.out.print("Processing " + i + " of " + seqs.size() + " sequences...\r");
            for (int j = 0; j < i; ++ j) {
                int d = snp(seqs.get(i), seqs.get(j));
                if (d <= this.edgeThreshold) {
                    if (use_stdout)
                        System.out.println(String.format("%s,%s,%d", seqs.get(i).getName(), seqs.get(j).getName(), d));
                    else
                        f.println(String.format("%s,%s,%d", seqs.get(i).getName(), seqs.get(j).getName(), d));
                }
                ++current_pair;
                update_percent_complete(pairs_count, current_pair, startTime);
            }
        }
        if (!use_stdout) {
            f.flush();
            f.close();
        }
        setChanged();
        notifyObservers(100);
        return;
    }

    public void snp_parallel(ArrayList<Seq> seqs){
        if (this.cores >= seqs.size()) {
            this.cores = seqs.size()-1;
        }

        if (!use_stdout) System.out.println("Creating thread pool with " + this.cores + " threads...");
        ExecutorService executor = Executors.newFixedThreadPool(this.cores);
        List<Future<Triplet<Integer, Integer, Integer>>> futures = new ArrayList<>();
        AtomicReference<PrintWriter> f = new AtomicReference<>();

        if (use_stdout)
            System.out.println("Source,Target,Distance");
        else {
            try {
                f.set(new PrintWriter(new BufferedWriter(new FileWriter(outputFile))));
            } catch (IOException e) {
                e.printStackTrace();
            }
            f.get().println("Source,Target,Distance");
        }

        long pairs_count = (seqs.size() * (seqs.size() - 1)) / 2;
        long current_pair = 0;
        long startTime = System.nanoTime(), estimatedTime;

        for (int i = 1; i < seqs.size(); ++i) {
            if (!use_stdout) System.out.print("Processing " + i + " of " + seqs.size() + " sequences...\r");
            for (int j = 0; j < i; ++ j) {
                final int row = i;
                final int col = j;
                final Seq seq1 = seqs.get(row);
                final Seq seq2 = seqs.get(col);
                final int L = seq1.getSeq().length();

                futures.add(executor.submit( () -> {
                    int d = snp(seq1, seq2);
                    if (d <= this.edgeThreshold) {
                        if (use_stdout)
                            System.out.println(String.format("%s,%s,%d", seq1.getName(), seq2.getName(), d));
                        else
                            f.get().println(String.format("%s,%s,%d", seq1.getName(), seq2.getName(), d));
                    }
                    return new Triplet<>(row, col, d);
                }));


                if (futures.size() > 1000 || pairs_count < 1000) {
                    for(Future<Triplet<Integer, Integer, Integer>> future: futures) {
                        try {
                            future.get();
                        } catch (InterruptedException | ExecutionException e) {
                            e.printStackTrace();
                        }
                        ++current_pair;
                        update_percent_complete(pairs_count, current_pair, startTime);
                    }
                    futures.clear();
                }
            }
        }

        for (Future<Triplet<Integer, Integer, Integer>> future: futures) {
            try {
                future.get();
            } catch (InterruptedException | ExecutionException e) {
                e.printStackTrace();
            }
            ++current_pair;
            update_percent_complete(pairs_count, current_pair, startTime);
        }

        executor.shutdown();
        if (!use_stdout) {
            f.get().flush();
            f.get().close();
        }
        setChanged();
        notifyObservers(100);
        return;
    }

    public int snp(Seq s1, Seq s2){
        int d = 0;
        int[] seq1 = s1.getSeq_enc();
        int[] seq2 = s2.getSeq_enc();
        int start_pos = 0;
        int end_pos = Math.min(seq1.length, seq2.length);
        if (this.ignoreTerminalGaps) {
            int s1s = findBeginningGapsEnd(seq1);
            int s1e = findTerminalGapsStart(seq1);
            int s2s = findBeginningGapsEnd(seq2);
            int s2e = findTerminalGapsStart(seq2);
            start_pos = Math.max(s1s, s2s);
            end_pos = Math.min(s1e, s2e);
        }
        for (int i = start_pos; i < end_pos; ++i) {
            if (seq1[i] != seq2[i]) {
                if (this.ignoreAllGaps && (seq1[i] == 17 || seq2[i] == 17)) continue;
                if (this.ignoreAmbiguities && (seq1[i] > 4 || seq2[i] > 4)) continue;
                d++;
            }
        }
        return d;
    }

    private void update_percent_complete(long total_pairs_to_compute, long current_pair, long startTime) {
        long estimatedTime;
        int percCompleted;

        if (total_pairs_to_compute < 100) {
            estimatedTime = System.nanoTime() - startTime;
            percCompleted = (int) (current_pair / total_pairs_to_compute);
        } else if (current_pair % (total_pairs_to_compute / 100) == 0) {
            estimatedTime = System.nanoTime() - startTime;
            percCompleted = (int) (current_pair*100/total_pairs_to_compute);
        }
        else 
            return;

        if (!use_stdout) System.out.print(String.format("%d%% completed in ", percCompleted));
        if (!use_stdout) System.out.print(TimeUnit.SECONDS.convert(estimatedTime, TimeUnit.NANOSECONDS));
        if (!use_stdout) System.out.println(" sec                                ");
        setChanged();
        notifyObservers(percCompleted);
    }

    private static int findBeginningGapsEnd(String seq) {
        int i = 0;
        while (i < seq.length() && seq.charAt(i) == '-') {
            ++i;
        }
        return i;
    }

    private static int findBeginningGapsEnd(int[] seq) {
        int i = 0;
        while (i < seq.length && seq[i] == 17) {
            ++i;
        }
        return i;
    }

    private static int findTerminalGapsStart(String seq) {
        int i = seq.length() - 1;
        while (i >= 0 && seq.charAt(i) == '-') {
            --i;
        }
        return i;
    }

    private static int findTerminalGapsStart(int[] seq) {
        int i = seq.length - 1;
        while (i >= 0 && seq[i] == 17) {
            --i;
        }
        return i;
    }

    private static ArrayList<Seq> read_fasta(File inputFile) throws FileNotFoundException {
        Scanner sc = new Scanner(inputFile);
        ArrayList<Seq> a = read_seqs(sc);
        sc.close();
        return a;
    }
    public static ArrayList<Seq> read_seqs(Scanner sc) {
        ArrayList<Seq> seqs = new ArrayList<Seq>();
        String name="", seq="";
        while(sc.hasNextLine()) {
            String line = sc.nextLine().trim();
            if(line.length() == 0) continue;
            if(line.charAt(0)=='>') {
                if (name.length()!=0) seqs.add(new Seq(name, seq));
                name = line.substring(1);
                seq="";
            }
            else seq=seq.concat(line);
        }
        if(name.length()!=0) seqs.add(new Seq(name, seq));
        return seqs;
    }
    private static ArrayList<Seq> read_fasta_stdin() {
        Scanner sc = new Scanner(System.in);
        ArrayList<Seq> a = read_seqs(sc);
        sc.close();
        return a;
    }
}
