package SNP;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicReference;
import java.util.Observable;

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

    public void snpFasta() {
        System.out.println("Running SNPFASTA");
        PrintWriter f = null;
        try {
            System.out.println("Reading input file...");
            ArrayList<Seq> seqs = read_fasta(inputFile);
            System.out.println("Calculating distances...");
            snp(seqs);
        }
        catch(FileNotFoundException e) {
            e.printStackTrace();
        }
        finally {
            if(f != null) f.close();
        }
    }


    public void snp(ArrayList<Seq> seqs){
        if (this.cores >= seqs.size()) {
            this.cores = seqs.size()-1;
        }

        System.out.println("Creating thread pool with " + this.cores + " threads...");
        ExecutorService executor = Executors.newFixedThreadPool(this.cores);
        List<Future<Triplet<Integer, Integer, Integer>>> futures = new ArrayList<>();
        AtomicReference<PrintWriter> f = new AtomicReference<>();
        try {
            f.set(new PrintWriter(outputFile));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        f.get().println("Source,Target,Distance");
        long pairs_count = (seqs.size() * (seqs.size() - 1)) / 2;
        long current_pair = 0;
        long startTime = System.nanoTime(), estimatedTime;

        for (int i = 1; i < seqs.size(); ++i) {
            System.out.print("Processing " + i + " of " + seqs.size() + " sequences...\r");
            for (int j = 0; j < i; ++ j) {
                final int row = i;
                final int col = j;
                final Seq seq1 = seqs.get(row);
                final Seq seq2 = seqs.get(col);
                final int L = seq1.getSeq().length();

                futures.add(executor.submit( () -> {
                    int d = snp(seq1, seq2);
                    if ((float) d / L <= this.edgeThreshold) {
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
                        if (pairs_count < 100 || current_pair % (pairs_count / 100) == 0) {
                            estimatedTime = System.nanoTime() - startTime;
                            int percCompleted = (int) (current_pair*100/pairs_count);
                            System.out.print(String.format("%d%% completed in ", percCompleted));
                            System.out.print(TimeUnit.SECONDS.convert(estimatedTime, TimeUnit.NANOSECONDS));
                            System.out.println(" sec                                ");
                            setChanged();
                            notifyObservers(percCompleted);
                        }
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
            if (pairs_count < 100 || current_pair % (pairs_count / 100) == 0) {
                estimatedTime = System.nanoTime() - startTime;
                int percCompleted = (int) (current_pair*100/pairs_count);
                System.out.print(String.format("%d%% completed in ", percCompleted));
                System.out.print(TimeUnit.SECONDS.convert(estimatedTime, TimeUnit.NANOSECONDS));
                System.out.println(" sec                                ");
                setChanged();
                notifyObservers(percCompleted);
            }
        }

        executor.shutdown();
        f.get().flush();
        f.get().close();
        setChanged();
        notifyObservers(100);
        return;
    }

    public int snp(Seq s1, Seq s2){
        int d = 0;
        String seq1 = s1.getSeq();
        String seq2 = s2.getSeq();
        for (int i = 0; i < seq1.length(); ++i) {
            if (seq1.charAt(i) != seq2.charAt(i) && seq1.charAt(i) != '-' && seq2.charAt(i) != '-') {
                ++d;
            }
        }
        return d;
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
}
