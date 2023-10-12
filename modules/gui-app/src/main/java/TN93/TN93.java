package TN93;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.*;
import java.util.concurrent.atomic.AtomicReference;
//Parallelization
import java.util.concurrent.*;




import static java.lang.Math.log;

public class TN93 extends Observable {
    // Used for parallelization
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

    public static class Pair<A, B> {
        final A first;
        final B second;
    
        public Pair(A first, B second) {
            this.first = first;
            this.second = second;
        }
    }


    private float edgeThreshold = 1;
    private File inputFile;
    private File outputFile;
    private String ambiguityHandling;
    private int cores = 1;
    private double max_ambiguity_fraction = -1; // -1 means no limit (default)
    private boolean use_stdin = false;
    private boolean use_stdout = false;
    private boolean input_as_pairs = false;
    private boolean enumerate_sequences = false;
    private HashMap<Integer, String> map_index_to_name = new HashMap<Integer, String>();

    public void setEdgeThreshold(float edgeThreshold) {
        this.edgeThreshold = edgeThreshold;
    }

    public void setUseStdin(boolean use_stdin) {
        this.use_stdin = use_stdin;
    }

    public void setUseStdout(boolean use_stdout) {
        this.use_stdout = use_stdout;
    }

    public void setEnumerateSequences(boolean enumerate_sequences) {
        this.enumerate_sequences = enumerate_sequences;
    }

    public void setInputFile(File inputFile) {
        this.inputFile = inputFile;
    }

    public void setOutputFile(File outputFile) {
        this.outputFile = outputFile;
    }

    public void setAmbiguityHandling(String ambiguityHandling) {
        this.ambiguityHandling = ambiguityHandling;
    }

    public void setMaxAmbiguityFraction(double max_ambiguity_fraction) {
        this.max_ambiguity_fraction = max_ambiguity_fraction;
    }

    public void setCores(int cores) {
        this.cores = cores;
    }

    public void setInputAsPairs(boolean input_as_pairs) {
        this.input_as_pairs = input_as_pairs;
    }

    final static int[][] resolutions = {
        // A,C,G,T
        {1, 0, 0, 0},  // A             -> A (0) (Adenine)
        {0, 1, 0, 0},  // C             -> C (1) (Cytosine)
        {0, 0, 1, 0},  // G             -> G (2) (Guanine)
        {0, 0, 0, 1},  // T             -> T (3) (Thymine)
        {0, 0, 0, 1},  // T             -> U (4) (Uracil)
        {1, 0, 1, 0},  // A | G         -> R (5) (Either Purine)
        {0, 1, 0, 1},  // C | T         -> Y (6) (Either Pyrimidine)
        {0, 1, 1, 0},  // C | G         -> S (7)
        {1, 0, 0, 1},  // A | T         -> W (8)
        {0, 0, 1, 1},  // G | T         -> K (9)
        {1, 1, 0, 0},  // A | C         -> M (10)
        {0, 1, 1, 1},  // C | G | T     -> B (11) (Not Adenine)
        {1, 0, 1, 1},  // A | G | T     -> D (12) (Not Cytosine)
        {1, 1, 0, 1},  // A | C | T     -> H (13) (Not Guanine)
        {1, 1, 1, 0},  // A | C | G     -> V (14) (Not Thymine)
        {1, 1, 1, 1},  // A | C | G | T -> N (15)
        {1, 1, 1, 1},  // A | C | G | T -> ? (16)
        {0, 0, 0, 0},  // GAP           -> - (17)
    };

    final static double[] resolutionsCount = {
        1.0,  // A
        1.0,  // C
        1.0,  // G
        1.0,  // T
        1.0,  // U
        1.0 / 2.0,  // R
        1.0 / 2.0,  // Y
        1.0 / 2.0,  // S
        1.0 / 2.0,  // W
        1.0 / 2.0,  // K
        1.0 / 2.0,  // M
        1.0 / 3.0,  // B
        1.0 / 3.0,  // D
        1.0 / 3.0,  // H
        1.0 / 3.0,  // V
        1.0 / 4.0,  // N
        1.0 / 4.0,  // ?
        0.0,  // GAP
    };

    public void tn93Fasta() {
        PrintWriter f = null;
        try {
            if (!use_stdout) System.out.println("Reading input file...");
            ArrayList<Seq> seqs;
            if (input_as_pairs) {
                ArrayList<Pair<Pair<String, Seq>, Pair<String, Seq>>> pairs = read_seq_pairs();
                tn93_pairs(pairs);
                return;
            }
            else if (use_stdin) {
                seqs = read_fasta_stdin();
            } else {
                seqs = read_fasta(inputFile);
            }
            if (!use_stdout) System.out.println("Calculating distances...");
            tn93(seqs);
        }
        catch(FileNotFoundException e) {
            e.printStackTrace();
        }
        finally {
            if(f != null) f.close();
        }
    }

    public void tn93(ArrayList<Seq> seqs) {
        if (cores > 1)
            tn93_parallel(seqs);
        else 
            tn93_sequential(seqs);
    }
    
    public void tn93_pairs(ArrayList<Pair<Pair<String, Seq>, Pair<String, Seq>>> pairs) {
        if (cores > 1)
            tn93_parallel_pairs(pairs);
        else 
            tn93_sequential_pairs(pairs);
    }

    public void tn93_parallel(ArrayList<Seq> seqs){
        if (cores >= seqs.size()) {
            cores = seqs.size()-1;
        }

        if (!use_stdout) System.out.println("Creating thread pool with " + cores + " threads...");
        ExecutorService executor = Executors.newFixedThreadPool(cores);
        List<Future<Triplet<Integer, Integer, Double>>> futures = new ArrayList<>();
        AtomicReference<PrintWriter> writerRef = new AtomicReference<>();


        if (use_stdout) 
            System.out.println("Source,Target,Distance");
        else {
            try {
                writerRef.set(new PrintWriter(new BufferedWriter(new FileWriter(outputFile))));
            } catch (IOException e) {
                e.printStackTrace();
            }
            writerRef.get().println("Source,Target,Distance");
        }
            

        long total_pairs_to_compute = ((long) seqs.size() * seqs.size() - seqs.size())/2;
        long current_pair = 0;
        long startTime = System.nanoTime(), estimatedTime;

        if (!use_stdout) System.out.println("Submitting jobs...");
        for (int i = 1; i < seqs.size(); ++i) {
            if (!use_stdout) System.out.print("Processing " + i + " of " + seqs.size() + " sequences...\r");
            for (int j = 0; j < i; ++ j) {
                final int row = i;
                final int col = j;
                final Seq seq1 = seqs.get(row);
                final Seq seq2 = seqs.get(col);
                
                futures.add(executor.submit( () -> {
                    double d = tn93(seq1, seq2);
                    String report = "";
                    if (d == -0) d = 0;
                    if (d < this.edgeThreshold){
                        if (enumerate_sequences) {
                            report = String.format("%d,%d,%f", row, col, d);
                        }
                        else
                            report = String.format("%s,%s,%f", seq1.getName(), seq2.getName(), d);

                        if (use_stdout)
                            System.out.println(report);
                        else
                            writerRef.get().println(report);
                    }
                    return new Triplet<>(row, col, d);
                }));

                if (futures.size() > 1000 || total_pairs_to_compute < 1000 ) {
                    for (Future<Triplet<Integer, Integer, Double>> future : futures) {
                        try {
                            future.get();
                        }
                        catch (InterruptedException | ExecutionException e) {
                            e.printStackTrace();
                        }
        
                        current_pair++;
                        update_percent_complete(total_pairs_to_compute, current_pair, startTime); 
                    }
                    futures.clear();
                }
            }
        }

        for (Future<Triplet<Integer, Integer, Double>> future : futures) {
            try {
                future.get();         
            }
            catch (InterruptedException | ExecutionException e) {
                e.printStackTrace();
            }

            current_pair++;
            update_percent_complete(total_pairs_to_compute, current_pair, startTime); 
        }

        executor.shutdown();
        if (!use_stdout) {
            writerRef.get().flush();
            writerRef.get().close();
        }

        //output a mapfile between index and name
        if (enumerate_sequences) {
            try {
                PrintWriter mapfile = new PrintWriter(new BufferedWriter(new FileWriter(outputFile + ".map")));
                for (int i = 0; i < seqs.size(); ++i) {
                    mapfile.println(i + "," + seqs.get(i).getName());
                }
                mapfile.flush();
                mapfile.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        setChanged();
        notifyObservers(100);
        return;
    }

    public void tn93_sequential(ArrayList<Seq> seqs) {
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
          
        long total_pairs_to_compute = (seqs.size() * seqs.size() - seqs.size())/2;
        long current_pair = 0;
        long startTime = System.nanoTime(), estimatedTime;

        for (int i = 1; i < seqs.size(); ++i) {
            if (!use_stdout) System.out.print("Processing " + i + " of " + seqs.size() + " sequences...\r");
            for (int j = 0; j < i; ++ j) {
                double distance = tn93(seqs.get(i), seqs.get(j));
                String report = "";

                if (distance == -0) distance = 0;
                if (distance <= edgeThreshold) {
                    if (enumerate_sequences) {
                        report = String.format("%d,%d,%f", i, j, distance);
                        }
                    else
                        report = String.format("%s,%s,%f", seqs.get(i).getName(), seqs.get(j).getName(), distance);
                    
                    if (use_stdout)
                        System.out.println(report);
                    else
                        f.println(report);
                }
                current_pair++;
                update_percent_complete(total_pairs_to_compute, current_pair, startTime);
            }
        }
        if (!use_stdout) {
            f.flush();
            f.close();
        }

        //output a mapfile between index and name
        if (enumerate_sequences) {
            try {
                PrintWriter mapfile = new PrintWriter(new BufferedWriter(new FileWriter(outputFile + ".map")));
                for (int i = 0; i < seqs.size(); ++i) {
                    mapfile.println(i + "," + seqs.get(i).getName());
                }
                mapfile.flush();
                mapfile.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        setChanged();
        notifyObservers(100);
        return;
    }

    public void tn93_parallel_pairs(ArrayList<Pair<Pair<String, Seq>, Pair<String, Seq>>> pairs) {
        if (cores >= pairs.size()) {
            cores = pairs.size()-1;
        }

        if (!use_stdout) System.out.println("Creating thread pool with " + cores + " threads...");
        ExecutorService executor = Executors.newFixedThreadPool(cores);
        List<Future<Triplet<Integer, Integer, Double>>> futures = new ArrayList<>();
        AtomicReference<PrintWriter> writerRef = new AtomicReference<>();


        if (use_stdout) 
            System.out.println("Source,Target,Distance");
        else {
            try {
                writerRef.set(new PrintWriter(new BufferedWriter(new FileWriter(outputFile))));
            } catch (IOException e) {
                e.printStackTrace();
            }
            writerRef.get().println("Source,Target,Distance");
        }

        long total_pairs_to_compute = pairs.size();
        long current_pair = 0;
        long startTime = System.nanoTime(), estimatedTime;

        if (!use_stdout) System.out.println("Submitting jobs...");

        for (Pair<Pair<String, Seq>, Pair<String, Seq>> pair : pairs) {
            final Pair<String, Seq> seq1 = pair.first;
            final Pair<String, Seq> seq2 = pair.second;
            futures.add(executor.submit( () -> {
                double d = tn93(seq1.second, seq2.second);
                if (d == -0) d = 0;
                if (d < this.edgeThreshold)
                    if (use_stdout)
                        System.out.println(String.format("%s,%s,%f", seq1.first, seq2.first, d));
                    else
                        writerRef.get().println(String.format("%s,%s,%f", seq1.first, seq2.first, d));
                return new Triplet<>(0, 0, d);
            }));

            if (futures.size() > 1000 || total_pairs_to_compute < 1000 ) {
                for (Future<Triplet<Integer, Integer, Double>> future : futures) {
                    try {
                        future.get();
                    }
                    catch (InterruptedException | ExecutionException e) {
                        e.printStackTrace();
                    }
    
                    current_pair++;
                    update_percent_complete(total_pairs_to_compute, current_pair, startTime); 
                }
                futures.clear();
            }
        }

        for (Future<Triplet<Integer, Integer, Double>> future : futures) {
            try {
                future.get();         
            }
            catch (InterruptedException | ExecutionException e) {
                e.printStackTrace();
            }

            current_pair++;
            update_percent_complete(total_pairs_to_compute, current_pair, startTime); 
        }

        executor.shutdown();
        if (!use_stdout) {
            writerRef.get().flush();
            writerRef.get().close();
        }
        setChanged();
        notifyObservers(100);
        return;
    }


    public void tn93_sequential_pairs(ArrayList<Pair<Pair<String, Seq>, Pair<String, Seq>>> pairs) {
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
          
        long total_pairs_to_compute = pairs.size();
        long current_pair = 0;
        long startTime = System.nanoTime(), estimatedTime;


        for (Pair<Pair<String, Seq>, Pair<String, Seq>> pair : pairs) {
            final Pair<String, Seq> seq1 = pair.first;
            final Pair<String, Seq> seq2 = pair.second;
            double distance = tn93(seq1.second, seq2.second);
            if (distance == -0) distance = 0;
            if (distance <= edgeThreshold) {
                if (use_stdout)
                    System.out.println(seq1.first + "," + seq2.first + "," + distance);
                else
                    f.println(seq1.first + "," + seq2.first + "," + distance);
            }
            current_pair++;
            update_percent_complete(total_pairs_to_compute, current_pair, startTime);
        }
        if (!use_stdout) {
            f.flush();
            f.close();
        }
        setChanged();
        notifyObservers(100);
        return;
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


    private boolean exceeds_ambig_threshold(Seq s1, Seq s2, double max_ambiguity_fraction) {
        int ambigs_count = 0;
        int scan_length = Math.min(s1.getSeq_enc().length, s2.getSeq_enc().length);
        int total_non_gap = 0;
        for (int i = 0; i < scan_length; ++i) {
            int c1 = s1.getSeq_enc()[i];
            int c2 = s2.getSeq_enc()[i];
            boolean c1_is_ambig = c1 > 4 && c1 != 17;
            boolean c2_is_ambig = c2 > 4 && c2 != 17;
            if (c1_is_ambig || c2_is_ambig) {
                ambigs_count++;
            }
            else if (!c1_is_ambig && !c2_is_ambig && c1 != 17 && c2 != 17) {
                total_non_gap++;
            }
        }
        return total_non_gap * max_ambiguity_fraction <= ambigs_count;
    }

    private double tn93(Seq s1, Seq s2) {
        double[][] nucl_pair_counts = countPairwiseNucl(s1, s2);
        double[] nucl_counts = getNuclCounts(nucl_pair_counts);

        double total_non_gap = 2.0 / Arrays.stream(nucl_counts).sum();

        double[] nucl_freq = new double[4];
        double auxd = 1.0 / Arrays.stream(nucl_counts).sum();
        for(int i=0; i<4; ++i) 
            nucl_freq[i] = nucl_counts[i] * auxd;

        double dist = 0.0;
        double AG_counts = nucl_pair_counts[Seq.A][Seq.G] + nucl_pair_counts[Seq.G][Seq.A];
        double AG = AG_counts * total_non_gap;
        double CT_counts = nucl_pair_counts[Seq.C][Seq.T] + nucl_pair_counts[Seq.T][Seq.C];
        double CT = CT_counts * total_non_gap;
        double matching = (nucl_pair_counts[Seq.A][Seq.A] + nucl_pair_counts[Seq.C][Seq.C] + nucl_pair_counts[Seq.T][Seq.T] + nucl_pair_counts[Seq.G][Seq.G]) * total_non_gap;
        double tv = 1 - (AG + CT + matching);

        boolean useK2P = false;
        for(int i=0; i<4; ++i) if (nucl_freq[i] == 0.0) useK2P = true;

        if (useK2P) {
            AG = 1 - 2 * (AG + CT) - tv;
            CT = 1 - 2 * tv;
            if (AG > 0 && CT > 0) {
                dist = -0.5 * log(AG) - 0.25 * log(CT);
            } else {
                dist = 1.0;
            }
        } else {
            double fR = nucl_freq[Seq.A] + nucl_freq[Seq.G];
            double fY = nucl_freq[Seq.C] + nucl_freq[Seq.T];
            double K1 = 2 * nucl_freq[Seq.A] * nucl_freq[Seq.G] / fR;
            double K2 = 2 * nucl_freq[Seq.C] * nucl_freq[Seq.T] / fY;
            double K3 = 2 * ( fR * fY - nucl_freq[Seq.A] * nucl_freq[Seq.G] * fY / fR - nucl_freq[Seq.C] * nucl_freq[Seq.T] * fR / fY);
            AG = 1 - AG / K1 - 0.5 * tv / fR;
            CT = 1 - CT / K2 - 0.5 * tv / fY;
            tv = 1 - 0.5 * tv / fY / fR;
            dist = -K1 * log(AG) - K2 * log(CT) - K3 * log(tv);
        }
        return dist;
    }


    private double[][] countPairwiseNucl(Seq s1, Seq s2) {
        double nucl_pair_counts[][] = new double[4][4];
        if ("resolve".equals(ambiguityHandling)) {
            if (max_ambiguity_fraction != -1)
                if (exceeds_ambig_threshold(s1, s2, max_ambiguity_fraction))
                    return countNucl_average(s1.getSeq_enc(), s2.getSeq_enc(), nucl_pair_counts);
            return countNucl_resolve(s1.getSeq_enc(), s2.getSeq_enc(), nucl_pair_counts);
        }
        else if ("average".equals(ambiguityHandling))
            return countNucl_average(s1.getSeq_enc(), s2.getSeq_enc(), nucl_pair_counts);
        else if ("gapmm".equals(ambiguityHandling))
            return countNucl_gapmm(s1.getSeq_enc(), s2.getSeq_enc(), nucl_pair_counts);
        else if ("skip".equals(ambiguityHandling))
            return countNucl_skip(s1.getSeq_enc(), s2.getSeq_enc(), nucl_pair_counts);
        else
        {
            if (max_ambiguity_fraction != -1)
                if (exceeds_ambig_threshold(s1, s2, max_ambiguity_fraction))
                    return countNucl_average(s1.getSeq_enc(), s2.getSeq_enc(), nucl_pair_counts);
            return countNucl_resolve(s1.getSeq_enc(), s2.getSeq_enc(), nucl_pair_counts);
        }
    }


    private static double[] getNuclCounts(double[][] nucl_pair_counts) {
        double[] nucl_freq = new double[4];
        for(int i=0; i<4; ++i) {
            for(int j=0; j<4; ++j) {
                nucl_freq[i] += nucl_pair_counts[i][j];
                nucl_freq[j] += nucl_pair_counts[i][j];
            }
        }
        return nucl_freq;
    }


    private double[][] countNucl_resolve(int[] s1, int[] s2, double[][]nucl_pair_counts) {
        int L = Math.min(s1.length, s2.length);
        for(int i=0; i<L; i++) {
            int c1 = s1[i];
            int c2 = s2[i];
            if (c1 == 17 && c2 == 17) continue;                 // both are gaps; continue
            if (c1 < 4 && c2 < 4) {                             // Neither is ambiguous
                nucl_pair_counts[c1][c2]++;
            }
            else {
                if (c1 < 4) {                                   // if c1 is not ambiguous, c2 is
                    if (resolutionsCount[c2] > 0) {             // if c2 is not a gap
                        if (resolutions[c2][c1] == 1) {         // if c2 can resolve to c1               
                            nucl_pair_counts[c1][c1]++;          // Resolve c2 to c1 
                            continue;
                        }
                        for (int j=0; j<4; j++) {               // else: average over all possible resolutions
                            if (resolutions[c2][j] == 1) {
                                nucl_pair_counts[c1][j] += resolutionsCount[c2];
                            }
                        }
                    }
                }
                else if (c2 < 4){                               // c2 is not ambiguous, c1 is
                    if (resolutionsCount[c1] > 0) {
                        if (resolutions[c1][c2] == 1) {
                            nucl_pair_counts[c2][c2]++;
                            continue;
                        }
                        for (int j=0; j<4; j++) {
                            if (resolutions[c1][j] == 1) {
                                nucl_pair_counts[j][c2] += resolutionsCount[c1];
                            }
                        }
                    }
                } 
                else {                                          // Both c1 and c2 are ambiguous
                    double norm = resolutionsCount[c1] * resolutionsCount[c2]; 
                    if (norm > 0.0) { // if both are not gaps
                        int matched = 0;
                        boolean[] positive_match = new boolean[4];
                        for (int j=0; j<4; j++) { 
                            if (resolutions[c1][j] == 1 && resolutions[c2][j] == 1) {
                                positive_match[j] = true;
                                matched++;
                            }
                        }
                        if (matched > 0) {
                            double norm2 = 1.0/matched;
                            for (int k=0; k<4; k++) {
                                if (positive_match[k]) {
                                    nucl_pair_counts[k][k] += norm2;
                                }
                            }
                            continue;
                        }
                        for (int j=0; j<4; j++) {
                            if (resolutions[c1][j] == 1) {
                                for (int k=0; k<4; k++) {
                                    if (resolutions[c2][k] == 1) {
                                        nucl_pair_counts[j][k] += norm;
                                    }
                                }
                            }
                        }
                    }
                } 
            }
        }
        return nucl_pair_counts;
    }


    private double[][] countNucl_average(int[] s1, int[] s2, double[][] nucl_pair_counts) {
        int L = Math.min(s1.length, s2.length);
        for(int i=0; i<L; i++) {
            int c1 = s1[i];
            int c2 = s2[i];

            if (c1 == 17 || c2 == 17) continue;   

            if (c1 < 4 && c2 < 4) {
                nucl_pair_counts[c1][c2]++;
            }
            else {
                if (c1 < 4) {
                    if (resolutionsCount[c2] > 0) {
                        for (int j=0; j<4; j++) {
                            if (resolutions[c2][j] == 1) {
                                nucl_pair_counts[c1][j] += resolutionsCount[c2];
                            }
                        }
                    }
                }
                else if (c2 < 4){
                    if (resolutionsCount[c1] > 0) {
                        for (int j=0; j<4; j++) {
                            if (resolutions[c1][j] == 1) {
                                nucl_pair_counts[j][c2] += resolutionsCount[c1];
                            }
                        }
                    }
                } 
                else {
                    double norm = resolutionsCount[c1] * resolutionsCount[c2]; 
                    if (norm > 0.0) {
                        for (int j=0; j<4; j++) {
                            if (resolutions[c1][j]==1) {
                                for (int k=0; k<4; k++) {
                                    if (resolutions[c2][k]==1) {
                                        nucl_pair_counts[j][k] += norm;
                                    }
                                }
                            }
                        }
                    }
                } 
            }
        }
        return nucl_pair_counts;
    }


    private static double[][] countNucl_gapmm(int[] s1, int[] s2, double[][] nucl_pair_counts) {
        int L = Math.min(s1.length, s2.length);
        for (int i = 0; i < L; i++) {
            int c1 = s1[i];
            int c2 = s2[i];

            if (c1 == 17 && c2 == 17) continue;

            if (c1 < 4 && c2 < 4) 
                nucl_pair_counts[c1][c2]++;
            else {
                if (c1 == 17 || c2 == 17) {
                    if (c1 == 17) 
                        c1 = 15;
                    else 
                        c2 = 15;
                }
                if (c1 < 4) {
                    if (resolutionsCount[c2] > 0) {
                        for (int j=0; j<4; j++) {
                            if (resolutions[c2][j] == 1) {
                                nucl_pair_counts[c1][j] += resolutionsCount[c2];
                            }
                        }
                    }
                }
                else if (c2 < 4) {
                    if (resolutionsCount[c1] > 0) {
                        for (int j=0; j<4; j++) {
                            if (resolutions[c1][j] == 1) {
                                nucl_pair_counts[j][c2] += resolutionsCount[c1];
                            }
                        }
                    }
                }
                else {
                    double norm = resolutionsCount[c1] * resolutionsCount[c2]; 
                    if (norm > 0.0) {
                        for (int j=0; j<4; j++) {
                            if (resolutions[c1][j]==1) {
                                for (int k=0; k<4; k++) {
                                    if (resolutions[c2][k]==1) {
                                        nucl_pair_counts[j][k] += norm;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return nucl_pair_counts;
    }


    private static double[][] countNucl_skip(int[] s1, int[] s2, double[][] nucl_pair_counts) {
        int L = Math.min(s1.length, s2.length);

        for (int i = 0; i < L; i++) {
            int c1 = s1[i];
            int c2 = s2[i];

            if (c1 < 4 && c2 < 4) {
                nucl_pair_counts[c1][c2]++;
            }
        }
        return nucl_pair_counts;
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


    private static ArrayList<Seq> read_fasta(File inputFile) throws FileNotFoundException {
        Scanner sc = new Scanner(inputFile);
        ArrayList<Seq> a = read_seqs(sc);
        sc.close();
        return a;
    }


    private static ArrayList<Seq> read_fasta_stdin() {
        Scanner sc = new Scanner(System.in);
        ArrayList<Seq> a = read_seqs(sc);
        sc.close();
        return a;
    }


    private static ArrayList<Pair<Pair<String, Seq>, Pair<String, Seq>>> read_seq_pairs() {
        // This supports reading pairs of sequences from stdin
        // Expected format:
        // seq1_name, seq1, seq2_name, seq2\n

        Scanner sc = new Scanner(System.in);
        ArrayList<Pair<Pair<String, Seq>, Pair<String, Seq>>> pairs = new ArrayList<>();

        while(sc.hasNextLine()) {
            String line = sc.nextLine().trim();
            if(line.length() == 0) continue;
            String[] fields = line.split(",");
            if (fields.length != 4) {
                System.err.println("Error: expected 4 fields per line, got " + fields.length);
                System.exit(1);
            }
            String name1 = fields[0];
            String seq1 = fields[1];
            String name2 = fields[2];
            String seq2 = fields[3];

            pairs.add(new Pair<Pair<String, Seq>, Pair<String, Seq>>(
                new Pair<String, Seq>(name1, new Seq(name1, seq1)),
                new Pair<String, Seq>(name2, new Seq(name2, seq2))
            ));
        }
        sc.close();
        return pairs;
    }
}
