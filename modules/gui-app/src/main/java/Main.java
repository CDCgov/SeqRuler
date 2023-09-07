import TN93.TN93;
import SNP.SNP;
import picocli.CommandLine;

import javax.swing.*;
import javax.swing.event.DocumentListener;
import javax.swing.event.DocumentEvent;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.Observable;
import java.util.Observer;
import javax.swing.filechooser.FileNameExtensionFilter;

import static javax.swing.JOptionPane.showMessageDialog;

@CommandLine.Command(name = "SeqRuler", mixinStandardHelpOptions = true, version = "3.0")
public class Main implements Runnable{
    @CommandLine.Option(names={"-i", "--inFile"}, description="input file with sequences",
            paramLabel = "FILE")
    private File inputFile;
    @CommandLine.Option(names={"-o", "--outFile"},
            description="output file with distances",
            paramLabel = "FILE")
    private File outputFile;
    @CommandLine.Option(names={"-s", "--stdin"}, description="read fasta from stdin. Alternative to reading from a file (-i)", defaultValue = "false")
    private boolean use_stdin;
    @CommandLine.Option(names={"-S", "--stdout"}, description="write distances to stdout. Alternative to writing to a file (-o)", defaultValue = "false")
    private boolean use_stdout;
    @CommandLine.Option(names={"-p", "--pairs"}, description="read pairs of sequences from stdin, calculate distance for each pair. format \"name1, seq1, name2, seq2\\n\"", defaultValue = "false")
    private boolean input_as_pairs;
    @CommandLine.Option(names={"-d", "--distance-method"},
            description="distance metric to use. One of [TN93, SNP]. Default: TN93", defaultValue = "TN93")
    private String distanceMethod;
    @CommandLine.Option(names={"-r", "--run-server"}, description="run jetty server", defaultValue = "false")
    private boolean is_server;
    @CommandLine.Option(names={"-t", "--edge-threshold"},
            description="edges above the threshold are not reported in output", defaultValue = "1.0")
    private String edgeThresholdString;
    @CommandLine.Option(names={"-a", "--ambiguity", "--ambiguities"},
            description="How to handle ambiguous nucleotides. One of [resolve, average, gapmm, skip]", defaultValue = "resolve")
    private String ambiguityHandling;
    @CommandLine.Option(names={"-f", "--fraction"},
            description="Maximum allowable fraction of ambiguities allowed for 'resolve' mode. If exceeded, use 'average' mode.", defaultValue = "1.0")
    private float max_ambiguity_fraction;
    @CommandLine.Option(names={"-c", "--cores"},
            description="Number of cores to use for parallel processing.", defaultValue = "1")
    private int cores;
    @CommandLine.Option(names={"-g", "--ignore-terminal-gaps"},
            description="Ignore terminal gaps at beginning and end of sequences when calculating distances. [SNP only] Default: true")
    private boolean ignoreTerminalGaps=true; //snp only
    @CommandLine.Option(names={"-G", "--ignore-all-gaps"},
            description="Ignore all gaps when calculating distances. [SNP only] Default: false")
    private boolean ignoreAllGaps=false; //snp only



    public void run() {
        if(is_server) {
            try {
                run_server();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }
        if(inputFile == null && use_stdin == false && input_as_pairs == false || outputFile == null && use_stdout == false) {
            javax.swing.SwingUtilities.invokeLater(new Runnable() {
                public void run() {
                    createAndShowGUI();
                }
            });
        }
        else {
            if(distanceMethod == null) distanceMethod = "TN93";
            distanceMethod = distanceMethod.toUpperCase();
            if("TN93".equals(distanceMethod)) {
                TN93 tn93 = new TN93();
                tn93.setEdgeThreshold(Float.parseFloat(edgeThresholdString));
    
                if (use_stdin) 
                    tn93.setUseStdin(true);
                    if (input_as_pairs)
                        tn93.setInputAsPairs(true);
                else
                    tn93.setInputFile(inputFile);

                if (use_stdout)
                    tn93.setUseStdout(true);
                else
                    tn93.setOutputFile(outputFile);
                System.out.println(ambiguityHandling);
                tn93.setAmbiguityHandling(ambiguityHandling);
                tn93.setMaxAmbiguityFraction(max_ambiguity_fraction);
                tn93.setCores(cores);
                tn93.tn93Fasta();
            }
            else if("SNP".equals(distanceMethod)) {
                SNP snp = new SNP();

                if (use_stdin) 
                    snp.setUseStdin(true);
                    if (input_as_pairs)
                        snp.setInputAsPairs(true);
                else
                    snp.setInputFile(inputFile);

                if (use_stdout)
                    snp.setUseStdout(true);
                else
                    snp.setOutputFile(outputFile);

                snp.setEdgeThreshold(Float.parseFloat(edgeThresholdString));
                snp.setCores(cores);
                snp.setIgnoreTerminalGaps(ignoreTerminalGaps);
                snp.snpFasta();
            }
        }
    }
    private static void run_server() throws InterruptedException {
        System.out.println("To stop the server press Ctrl-C");
        Runtime.getRuntime().addShutdownHook(new Thread() {
            @Override
            public void run() {
                System.out.println("Stopping the server...");
            }
        });
        //TODO: Run Server
        while(true) Thread.sleep(1000);
    }

    public static void main(String[] args) {
        CommandLine.run(new Main(), System.out, args);
    }

    private static void createAndShowGUI() {
        //Create and set up the window.
        JFrame frame = new JFrame("SeqRuler");
        JPanel mainPane = new JPanel();
        mainPane.setLayout(new BoxLayout(mainPane, BoxLayout.Y_AXIS));
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        //Create and set up the content pane.
        TN93_Panel panel = new TN93_Panel(frame);
        mainPane.add(panel);
        panel.setOpaque(true); //content panes must be opaque
        frame.setContentPane(mainPane);

        //Display the window.
        frame.pack();
        frame.setVisible(true);
    }
}

class TN93_Panel extends JPanel implements ActionListener, Observer {
    private float edgeThreshold;
    private JToggleButton tn93ModeBut, snpModeBut;
    private JButton inBut, outBut, runBut;
    private JTextField fastaTextField, edgeListTextField, edgeThresholdField;
    private JProgressBar progress;
    private ButtonGroup modeGroup;
    private ButtonGroup ambiguityHandlingGroup;
    private JRadioButton resolveBut, averageBut, gapmmBut, skipBut;
    private JTextField maxAmbiguityFractionField;
    private JLabel maxEdgeThresholdLabel, maxAmbiguityFractionLabel;
    private JCheckBox useMaxCoresCheckbox;
    private JTextField numCoresField;
    private JCheckBox ignoreTerminalGapsCheckbox;
    private JCheckBox ignoreAllGapsCheckbox;


    private File fastaFile, edgeListFile;
    private TN93 tn93;
    private SNP snp;
    private JFrame frame;


    TN93_Panel(JFrame frame) {
        this.frame = frame;
        tn93 = new TN93();
        tn93.addObserver(this);
        snp = new SNP();
        snp.addObserver(this);

        modeGroup = new ButtonGroup();
        tn93ModeBut = new JToggleButton("TN93 Distance");
        snpModeBut = new JToggleButton("SNP Distance");
        modeGroup.add(tn93ModeBut);
        modeGroup.add(snpModeBut);
        tn93ModeBut.setSelected(true);


        inBut = new JButton("Load Fasta");
        outBut = new JButton("Save as: Edge List CSV File");
        runBut = new JButton("Run TN93");

        ambiguityHandlingGroup = new ButtonGroup();
        resolveBut = new JRadioButton("Resolve");
        averageBut = new JRadioButton("Average");
        gapmmBut = new JRadioButton("Gapmm");
        skipBut = new JRadioButton("Skip");

        ambiguityHandlingGroup.add(averageBut);
        ambiguityHandlingGroup.add(resolveBut);
        ambiguityHandlingGroup.add(gapmmBut);
        ambiguityHandlingGroup.add(skipBut);

        averageBut.setSelected(true);

        fastaTextField = new JTextField(20);
        edgeListTextField = new JTextField(20);
        edgeThresholdField = new JTextField("0.015");
        progress = new JProgressBar(0, 100);
        progress.setStringPainted(true);

        inBut.setActionCommand("loadFasta");
        outBut.setActionCommand("specifyEdgeListFile");
        runBut.setActionCommand("runTN93");

        inBut.addActionListener(this);
        outBut.addActionListener(this);
        runBut.addActionListener(this);

        setLayout(new BorderLayout());


        maxEdgeThresholdLabel = new JLabel("Maximum distance threshold:", JLabel.CENTER);
        maxEdgeThresholdLabel.setToolTipText("Distances greater than this value will not be included in output file");
        
        JPanel fastaPanel = new JPanel();
        fastaPanel.setLayout(new GridLayout(5, 2));
        fastaPanel.add(tn93ModeBut);
        fastaPanel.add(snpModeBut);
        fastaPanel.add(maxEdgeThresholdLabel);
        fastaPanel.add(edgeThresholdField);
        fastaPanel.add(fastaTextField);
        fastaPanel.add(inBut);
        fastaPanel.add(edgeListTextField);
        fastaPanel.add(outBut);
        fastaPanel.add(progress);
        fastaPanel.add(runBut);
        add(fastaPanel, BorderLayout.NORTH);


        JPanel ambigsPanel = new JPanel();   
        ambigsPanel.setLayout(new GridLayout(1, 1));
        ambigsPanel.setBorder(BorderFactory.createTitledBorder("Ambiguity Handling"));

        JPanel radioButtonsPanel = new JPanel();
        radioButtonsPanel.setLayout(new GridLayout(1, 4));
        resolveBut.setToolTipText("Count any resolutions that match as a perfect match");
        averageBut.setToolTipText("Average all possible resolutions");
        gapmmBut.setToolTipText("Count character-gap sites as 4-way mismatches, otherwise Average");
        skipBut.setToolTipText("Skip all sites with gaps or ambiguities");
        radioButtonsPanel.add(averageBut);
        radioButtonsPanel.add(resolveBut);
        radioButtonsPanel.add(gapmmBut);
        radioButtonsPanel.add(skipBut);
        ambigsPanel.add(radioButtonsPanel, BorderLayout.NORTH);

        JPanel maxAmbigsPanel = new JPanel();
        maxAmbigsPanel.setLayout(new GridLayout(1, 2));
        maxAmbiguityFractionLabel = new JLabel("Maximum ambiguity fraction: ", JLabel.CENTER);
        maxAmbiguityFractionLabel.setToolTipText("Sequences with a higher proportion of ambiguities than this value will instead be averaged");
        maxAmbigsPanel.add(maxAmbiguityFractionLabel);
        maxAmbiguityFractionField = new JTextField("0.05");
        maxAmbigsPanel.add(maxAmbiguityFractionField, BorderLayout.SOUTH);
        
        add(ambigsPanel, BorderLayout.CENTER);

        JPanel gapHandlingPanel = new JPanel();
        gapHandlingPanel.setLayout(new GridLayout(1, 2));
        gapHandlingPanel.setBorder(BorderFactory.createTitledBorder("Gap Handling"));
        ignoreTerminalGapsCheckbox = new JCheckBox("Ignore terminal gaps", true);
        ignoreTerminalGapsCheckbox.setToolTipText("Gaps at the beginning and end of sequences don't contribue to distance");
        gapHandlingPanel.add(ignoreTerminalGapsCheckbox);
        ignoreAllGapsCheckbox = new JCheckBox("Ignore all gaps", false);
        ignoreAllGapsCheckbox.setToolTipText("Gaps don't contribute to distance");
        gapHandlingPanel.add(ignoreAllGapsCheckbox);

        resolveBut.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                ambigsPanel.add(maxAmbigsPanel);
                ambigsPanel.setLayout(new GridLayout(2, 1));
                frame.pack();
            }
        });

        ActionListener hideMaxAmbiguityListener = new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                ambigsPanel.setLayout(new GridLayout(1, 1));
                ambigsPanel.remove(maxAmbigsPanel);
                frame.pack();
            }
        };

        averageBut.addActionListener(hideMaxAmbiguityListener);
        gapmmBut.addActionListener(hideMaxAmbiguityListener);
        skipBut.addActionListener(hideMaxAmbiguityListener);

        

        tn93ModeBut.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                runBut.setText("Run TN93");
                runBut.setActionCommand("runTN93");
                remove(gapHandlingPanel);
                add(ambigsPanel, BorderLayout.CENTER);
                revalidate();
                repaint();
                frame.pack();

            }
        });

        snpModeBut.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                remove(ambigsPanel);
                add(gapHandlingPanel, BorderLayout.CENTER);
                runBut.setText("Run SNP");
                runBut.setActionCommand("runSNP");
                revalidate();
                repaint();
                frame.pack();
            }
        });



        ignoreAllGapsCheckbox.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(ignoreAllGapsCheckbox.isSelected()) {
                    ignoreTerminalGapsCheckbox.setSelected(true);
                    ignoreTerminalGapsCheckbox.setEnabled(false);
                } else {
                    ignoreTerminalGapsCheckbox.setEnabled(true);
                }
            }
        });


        fastaTextField.getDocument().addDocumentListener(new DocumentListener() {
            @Override
            public void insertUpdate(DocumentEvent e) {
                fastaFile = new File(fastaFile.getParent(), fastaTextField.getText());
            }

            @Override
            public void removeUpdate(DocumentEvent e) {
                fastaFile = new File(fastaFile.getParent(), fastaTextField.getText());
            }

            @Override
            public void changedUpdate(DocumentEvent e) {
                fastaFile = new File(fastaFile.getParent(), fastaTextField.getText());
            }
        });

        edgeListTextField.getDocument().addDocumentListener(new DocumentListener() {
            @Override
            public void insertUpdate(DocumentEvent e) {
                edgeListFile = new File(edgeListFile.getParent(), edgeListTextField.getText());
            }

            @Override
            public void removeUpdate(DocumentEvent e) {
                edgeListFile = new File(edgeListFile.getParent(), edgeListTextField.getText());
            }

            @Override
            public void changedUpdate(DocumentEvent e) {
                edgeListFile = new File(edgeListFile.getParent(), edgeListTextField.getText());
            }
        });

        JPanel coresPanel = new JPanel();
        coresPanel.setLayout(new GridLayout(1, 2));
        useMaxCoresCheckbox = new JCheckBox("Use all available cores", true);
        numCoresField = new JTextField(Integer.toString(Runtime.getRuntime().availableProcessors()));
        numCoresField.setEnabled(false);

        useMaxCoresCheckbox.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if (useMaxCoresCheckbox.isSelected()) {
                    numCoresField.setEnabled(false);
                    numCoresField.setText(Integer.toString(Runtime.getRuntime().availableProcessors()));
                } else {
                    numCoresField.setEnabled(true);
                }
            }
        });

        coresPanel.add(new JLabel("Number of cores to use:", JLabel.CENTER));
        coresPanel.add(numCoresField);
        coresPanel.add(useMaxCoresCheckbox);

        coresPanel.setBorder(BorderFactory.createTitledBorder("Processing"));
        add(coresPanel, BorderLayout.SOUTH);
    }

    public void actionPerformed(ActionEvent e) {
        if("loadFasta".equals(e.getActionCommand())) {
            FileDialog fileDialog = new FileDialog(new Frame(), "Open Fasta file", FileDialog.LOAD);
            fileDialog.setFilenameFilter((dir, name) -> name.endsWith(".fa") || name.endsWith(".fas") || name.endsWith(".fasta"));
            fileDialog.setVisible(true);
        
            String selectedFileDirectory = fileDialog.getDirectory();
            String selectedFileName = fileDialog.getFile();
        
            if (selectedFileDirectory != null && selectedFileName != null) {
                fastaFile = new File(selectedFileDirectory, selectedFileName);
                fastaTextField.setText(fastaFile.getName());
            }
        }
        
        else if ("specifyEdgeListFile".equals(e.getActionCommand())) {
            FileDialog fileDialog = new FileDialog(new Frame(), "Save as: Edge List CSV File", FileDialog.SAVE);
            fileDialog.setFilenameFilter((dir, name) -> name.endsWith(".csv"));
            fileDialog.setVisible(true);
        
            String selectedFileDirectory = fileDialog.getDirectory();
            String selectedFileName = fileDialog.getFile();
        
            if (selectedFileDirectory != null && selectedFileName != null) {
                edgeListFile = new File(selectedFileDirectory, selectedFileName);
                edgeListTextField.setText(edgeListFile.getName());
            }
        }
        else if("runTN93".equals(e.getActionCommand())) {
            if(fastaFile == null) {
                showMessageDialog(null, "Specify an input Fasta file!");
                return;
            }
            if(edgeListFile == null) {
                showMessageDialog(null, "Specify an output file using Save as!");
                return;
            }
            if(!parseEdgeThreshold()) {
                showMessageDialog(null, "Threshold should be number!");
            }
            inactivatePanel();
            tn93.setEdgeThreshold(edgeThreshold);
            tn93.setInputFile(fastaFile);
            tn93.setOutputFile(edgeListFile);
            
            if(resolveBut.isSelected()) {
                tn93.setAmbiguityHandling("resolve");
                tn93.setMaxAmbiguityFraction(Double.parseDouble(maxAmbiguityFractionField.getText()));
            } else if(averageBut.isSelected()) 
                tn93.setAmbiguityHandling("average");
            else if(gapmmBut.isSelected()) 
                tn93.setAmbiguityHandling("gapmm");
            else if(skipBut.isSelected()) 
                tn93.setAmbiguityHandling("skip");
            else 
                tn93.setAmbiguityHandling("resolve");

            if(useMaxCoresCheckbox.isSelected()) {
                tn93.setCores(Runtime.getRuntime().availableProcessors());
            } else {
                try {
                    tn93.setCores(Integer.parseInt(numCoresField.getText()));
                } catch (NumberFormatException ex) {
                    showMessageDialog(null, "Number of cores should be number!");
                }
            }
            if(ignoreTerminalGapsCheckbox.isSelected()) {
                snp.setIgnoreTerminalGaps(true);
            } else {
                snp.setIgnoreTerminalGaps(false);
            }

            if(ignoreAllGapsCheckbox.isSelected()) {
                snp.setIgnoreAllGaps(true);
            } else {
                snp.setIgnoreAllGaps(false);
            }
        
            SwingWorker<Void, Void> worker = new SwingWorker<Void, Void>() {
                @Override
                protected Void doInBackground() {
                    tn93.tn93Fasta();
                    return null;
                }
                @Override
                protected void done() {
                    activatePanel();
                }
            };
            worker.execute();
        }
        else if("runSNP".equals(e.getActionCommand())) {
            if(fastaFile == null) {
                showMessageDialog(null, "Specify an input Fasta file!");
                return;
            }
            if(edgeListFile == null) {
                showMessageDialog(null, "Specify an output file using Save as!");
                return;
            }
            if(!parseEdgeThreshold()) {
                showMessageDialog(null, "Threshold should be number!");
            }
            inactivatePanel();
            snp.setInputFile(fastaFile);
            snp.setOutputFile(edgeListFile);
            snp.setEdgeThreshold(edgeThreshold);

            if(useMaxCoresCheckbox.isSelected()) {
                snp.setCores(Runtime.getRuntime().availableProcessors());
            } else {
                try {
                    snp.setCores(Integer.parseInt(numCoresField.getText()));
                } catch (NumberFormatException ex) {
                    showMessageDialog(null, "Number of cores should be number!");
                }
            }

            if(ignoreTerminalGapsCheckbox.isSelected()) {
                snp.setIgnoreTerminalGaps(true);
            } else {
                snp.setIgnoreTerminalGaps(false);
            }

            SwingWorker<Void,Void> worker = new SwingWorker<Void,Void>() {
                @Override
                protected Void doInBackground() {
                    snp.snpFasta();
                    return null;
                }
                @Override
                protected void done() {
                    activatePanel();
                }
            };
            worker.execute();
        }
    }
    private void inactivatePanel() {
        runBut.setEnabled(false);
        inBut.setEnabled(false);
        outBut.setEnabled(false);
        progress.setValue(0);
    }
    private void activatePanel() {
        runBut.setEnabled(true);
        inBut.setEnabled(true);
        outBut.setEnabled(true);
    }

    private boolean parseEdgeThreshold() {
        try {
            edgeThreshold = Float.parseFloat(edgeThresholdField.getText());
            return true;
        } catch (NumberFormatException e) {
            return false;
        }
    }
    public void update(Observable obj, Object arg) {
        this.progress.setValue((Integer) arg);
    }

}

class StatusPanel extends JPanel {
    protected JTextArea ta;

    StatusPanel() {
        ta = new JTextArea(10, 40);
        add(new JScrollPane(ta), BorderLayout.PAGE_START);
    }
}
