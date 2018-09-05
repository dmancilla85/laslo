/**
 *
 */
package com.laslo.core;

import com.tools.io.InputSequence;
import static java.lang.System.out;
import java.io.FileWriter;
import java.io.IOException;
import java.io.File;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Iterator;
import java.util.List;
import com.opencsv.CSVWriter;
import com.tools.UShuffle;
import com.tools.io.FASTACorrector;
import com.tools.io.SourceFile;
import static com.tools.io.SourceFile.*;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Locale;
import java.util.Map;
import java.util.ResourceBundle;
import static java.util.ResourceBundle.getBundle;
import java.util.Set;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.biojava.nbio.core.sequence.DNASequence;
import static org.biojava.nbio.core.sequence.io.FastaReaderHelper.readFastaDNASequence;
import static org.biojava.nbio.core.sequence.io.GenbankReaderHelper.readGenbankDNASequence;


/**
 * @author David
 *
 */
public class LoopCatcher {

    protected String pathOut;
    protected String pathIn;
    protected ArrayList<String> loopPatterns;
    protected InputSequence inputType;
    protected int minLength;
    protected int maxLength;
    protected int maxWooble;
    protected ResourceBundle bundle;
    protected int maxMismatch;
    protected File[] fileList;
    protected boolean extendedMode;
    private boolean makeRandoms;
    private int numberOfRandoms;
    private int kLetRandoms;
    protected File actualFile;
    protected String additionalSequence;

    public int getkLetRandoms() {
        return kLetRandoms;
    }

    public void setkLetRandoms(int kLetRandoms) {
        this.kLetRandoms = kLetRandoms;
    }

    public String getAdditionalSequence() {
        return additionalSequence;
    }

    public void setAdditionalSequence(String additionalSequence) {
        this.additionalSequence = additionalSequence;
    }

    // Create a new instance of RNAFold4J
    //protected RNAFoldAPI rfa;
    public LoopCatcher(String pathOut, String pathIn,
            ArrayList<String> loopPatterns, String additionalSequence,
            InputSequence inputType,
            int minLength, int maxLength,
            int maxWooble, int maxMismatch, Locale locale,
            int kLetRandoms) {
        this.pathOut = pathOut;
        this.loopPatterns = loopPatterns;
        this.inputType = inputType;
        this.minLength = minLength;
        this.maxLength = maxLength;
        this.maxWooble = maxWooble;
        this.maxMismatch = maxMismatch;
        //this.rfa = new RNAFoldAPI();
        this.fileList = null;
        this.extendedMode = false;
        this.makeRandoms = false;
        this.numberOfRandoms = 0;
        this.additionalSequence = additionalSequence;
        this.kLetRandoms = kLetRandoms;
        this.bundle = getBundle("resources/Bundle", locale);
    }

    /**
     *
     */
    public LoopCatcher() {
        this("", "", new ArrayList<>(), BiologyPatterns.PUM1,
                InputSequence.ENSEMBL, //NOI18N
                4, 16, 2, 0, new Locale("es", "AR"), 2);
    }

    /**
     *
     * @return
     */
    public File[] getFileList() {
        return fileList;
    }

    /**
     *
     * @param fileList
     */
    public void setFileList(File[] fileList) {
        this.fileList = fileList;
    }

    /**
     *
     * @param mode
     */
    public void setIsExtendedMode(boolean mode) {
        this.extendedMode = mode;
    }

    /**
     *
     * @return
     */
    public boolean getIsExtendedMode() {
        return this.extendedMode;
    }

    /**
     *
     * @return
     */
    public String getPathOut() {
        return pathOut;
    }

    /**
     *
     * @param pathOut
     */
    public void setPathOut(String pathOut) {
        this.pathOut = pathOut;
    }

    /**
     *
     * @return
     */
    public String getPathIn() {
        return pathIn;
    }

    /**
     *
     * @param pathIn
     */
    public void setPathIn(String pathIn) {

        File aux = new File(pathIn);

        if (aux.isDirectory()) {
            this.pathIn = pathIn;
        } else {
            this.pathIn = aux.getParent();

        }
    }

    /**
     *
     * @return
     */
    public ArrayList<String> getLoopPatterns() {
        return loopPatterns;
    }

    /**
     *
     * @param loopPatterns
     */
    public void setLoopPatterns(ArrayList<String> loopPatterns) {
        this.loopPatterns = loopPatterns;
    }

    public void setBundle(ResourceBundle bundle) {
        this.bundle = bundle;
    }

    /**
     *
     * @return
     */
    public InputSequence getInputType() {
        return inputType;
    }

    /**
     *
     * @param inputType
     */
    public void setInputType(InputSequence inputType) {
        this.inputType = inputType;
    }

    /**
     *
     * @return
     */
    public int getMinLength() {
        return minLength;
    }

    /**
     *
     * @param minLength
     */
    public void setMinLength(int minLength) {
        this.minLength = minLength;
    }

    /**
     *
     * @return
     */
    public int getMaxLength() {
        return maxLength;
    }

    public void setMaxLength(int maxLength) {
        this.maxLength = maxLength;
    }

    public int getMaxWooble() {
        return maxWooble;
    }

    public void setMaxWooble(int maxWooble) {
        this.maxWooble = maxWooble;
    }

    /**
     *
     * @return
     */
    public int getMaxMismatch() {
        return maxMismatch;
    }

    /**
     *
     * @param maxMismatch
     */
    public void setMaxMismatch(int maxMismatch) {
        this.maxMismatch = maxMismatch;
    }

    /**
     *
     * @param sFileName
     * @param stemResearch
     * @param headerList
     */
    public void writeCSV(String sFileName, List<StemLoop> stemResearch,
            String headerList) {

        CSVWriter writer;

        if (stemResearch.isEmpty() && headerList == null) {
            return;
        }

        try {
            writer = new CSVWriter(new FileWriter(sFileName), ';',
                    CSVWriter.DEFAULT_QUOTE_CHARACTER,
                    CSVWriter.DEFAULT_ESCAPE_CHARACTER,
                    CSVWriter.DEFAULT_LINE_END);

            if (headerList != null) {
                String[] header;
                header = headerList.split(";"); //NOI18N
                writer.writeNext(header);
            }

            for (int i = 0; i < stemResearch.size(); i++) {
                String[] data;
                data = stemResearch.get(i).toRowCSV().split(";"); //NOI18N
                writer.writeNext(data);
                //allData.add(data);
            }

            writer.close();

        } catch (IOException e) {
            // TODO Auto-generated catch block

        }
    }

    /**
     *
     * @param sequence
     * @param pattern
     * @return
     */
    /*public List<Integer> getPatternLocations(String sequence, String pattern) {
        List<Integer> locations = new ArrayList<>();
        String regExp = toRegularExpression(pattern);

        Pattern p = Pattern.compile(regExp);
        Matcher loopFinder = p.matcher(sequence);

        while (loopFinder.find()) {
            locations.add(loopFinder.start());
        }

        return locations;
    }*/

    // Union
    /**
     *
     * @param a
     * @param b
     * @return
     */
    public File[] unionFiles(File[] a, File[] b) {
        Set<File> set;
        set = new HashSet<>(Arrays.asList(a));
        set.addAll(Arrays.asList(b));
        return set.toArray(new File[set.size()]);
    }

    /**
     *
     * @return
     */
    public boolean beginSearch() {
        // To check the elapsed time
        Calendar ini, fin;
        ini = Calendar.getInstance();
        out.println(
                java.text.MessageFormat.format(bundle.getString("START_TIME"), new Object[]{Calendar.getInstance().getTime()}));
        out.flush();

        if (fileList == null) {
            return false;
        }

        if (this.makeRandoms) {
            for (File currentFile : fileList) {
                if (currentFile.isFile()
                        && (currentFile.toString().endsWith(FASTA_EXT)
                        || currentFile.toString().endsWith(FASTA_EXT_2)
                        || currentFile.toString().endsWith(GENBANK_EXT))) {

                    //fileName = currentFile.getName();
                    LinkedHashMap<String, DNASequence> fasta = null;

                    try {
                        if (currentFile.toString().endsWith(GENBANK_EXT)) {
                            fasta = readGenbankDNASequence(currentFile, false);
                        } else {
                            fasta = readFastaDNASequence(currentFile, false);
                        }

                    } catch (IOException ex) {
                        Logger.getLogger(LoopCatcher.class.getName()).log(Level.SEVERE, null, ex);
                    } catch (Exception ex) {
                        Logger.getLogger(LoopCatcher.class.getName()).log(Level.SEVERE, null, ex);
                    }

                    UShuffle.makeShuffleSequences(pathOut, currentFile.getName(),
                            fasta, numberOfRandoms, kLetRandoms);
                }
            }

            File folder;
            folder = new File(pathIn + UShuffle.getRandomDir());
            File[] randomFiles;
            randomFiles = folder.listFiles();
            fileList = unionFiles(fileList, randomFiles);
        }

        // I. File level (hacer un hilo)
        for (File currentFile : fileList) {
            if (currentFile.isFile() && (currentFile.toString().endsWith(FASTA_EXT)
                    || currentFile.toString().endsWith(FASTA_EXT_2)
                    || currentFile.toString().endsWith(GENBANK_EXT))) {

                this.actualFile = currentFile;
                processFile();
            }
        }

        fin = Calendar.getInstance();
        out.println(java.text.MessageFormat.format(bundle.getString("TOTAL_TIME"), new Object[]{((fin.getTimeInMillis() - ini.getTimeInMillis()) / 1000)}) + " s.");

        out.flush();

        // Memory cleaning
        fileList = null;
        return true;
    }

    public void processFile() {

        CSVWriter writer;
        String fileName, fileOut;
        final int MAX_HILOS = 100;
        int secuencias = 0, nHilos = MAX_HILOS, i, count = 0;
        ExecutorService pool;
        CountDownLatch latch;

        try {
            fileName = actualFile.getName();
            out.println(java.text.MessageFormat.format(bundle
                    .getString("FILE_PRINT"), new Object[]{fileName})); //$NON-NLS-1$
            out.flush();
            fileName = fileName.replaceFirst("[.][^.]+$", ""); //NOI18N
            fileOut = this.pathOut + "\\" + fileName + CSV_EXT; //NOI18N //NOI18N

            if (new File(fileOut).exists()) {
                try {
                    (new File(fileOut)).delete();

                } catch (Exception io) {
                    out.println(io.getMessage());
                    out.flush();
                    return;
                }
            }

            Calendar.getInstance();

            // Generation of the iterator of {id,sequence}
            LinkedHashMap<String, DNASequence> fasta;

            if (actualFile.getName().endsWith(GENBANK_EXT)) {
                fasta = readGenbankDNASequence(actualFile, false);
            } else {
                fasta = readFastaDNASequence(actualFile, false);
            }

            if (fasta.isEmpty()) {
                out.println(bundle.getString("INVALID_FILE_FORMAT"));
                out.println(bundle.getString("TRYING_TO_FIX"));
                boolean formatFile = FASTACorrector.formatFile(
                        actualFile.getAbsolutePath());
                if (formatFile) {
                    fasta = readFastaDNASequence(actualFile, false);
                } else {
                    out.println(bundle.getString("CANT_PROCESS"));
                    return;
                }
            }

            this.inputType = SourceFile.detectHeader(fasta.entrySet().iterator()
                    .next().getValue().getOriginalHeader());
            int listSize = fasta.size();
            //i = 1;

            writer = new CSVWriter(new FileWriter(fileOut), ';',
                    CSVWriter.DEFAULT_QUOTE_CHARACTER,
                    CSVWriter.DEFAULT_ESCAPE_CHARACTER,
                    CSVWriter.DEFAULT_LINE_END);
            writer.writeNext(StemLoop.getHeader(this.inputType).split(";")); //NOI18N

            // II. Transcript level
            //pool = Executors.newFixedThreadPool(listSize);
            //latch = new CountDownLatch(listSize);

            if(listSize < nHilos)
                nHilos = listSize;
            
            pool = Executors.newFixedThreadPool(nHilos);
            latch = new CountDownLatch(nHilos);
            i = 1;
            
            for (Map.Entry<String, DNASequence> entry : fasta.entrySet()) {
                
                DNASequence element = entry.getValue();               
                Iterator<String> patternItr = loopPatterns.iterator();
                count++;
                
                LoopCatcherThread thread = new LoopCatcherThread(extendedMode, 
                        additionalSequence, maxLength, minLength, element,
                        inputType, patternItr, writer);

                if(i <= nHilos){
                    thread.setLatch(latch);
                    pool.execute(thread);
                    i++;
                } else {
                    i = 1;
                    
                    if(count < nHilos)
                        nHilos = count;
                    
                    out.println("Esperando hijos...");
                    latch.await();
                    out.println("Terminando pool");
                    pool.shutdown();
                    pool = Executors.newFixedThreadPool(nHilos);
                    latch = new CountDownLatch(nHilos);
                    
                }
                        
                
                
                
            }
            //out.println("Esperando hijos...");
            //latch.await();
            //out.println("Terminando pool");
            //pool.shutdown();

            writer.close();
            writer = null;
            fasta = null;

            out.printf(bundle.getString("RESUME"),
                    this.minLength,
                    this.maxLength,
                    secuencias);
            out.flush();
        } catch (FileNotFoundException ex) {
            /*Logger.getLogger(LoopCatcher.class.getName())
                    .log(Level.SEVERE, null, ex);*/
            out.println(bundle.getString("CANT_OPEN_FILE"));
        } catch (IOException ex) {
            Logger.getLogger(LoopCatcher.class.getName())
                    .log(Level.SEVERE, null, ex);
        } catch (Exception ex) {
            Logger.getLogger(LoopCatcher.class.getName())
                    .log(Level.SEVERE, null, ex);
        }

        writer = null;
        fileName = null;
        fileOut = null;

        System.gc();
    }

    public boolean isMakeRandoms() {
        return makeRandoms;
    }

    public void setMakeRandoms(boolean makeRandoms) {
        this.makeRandoms = makeRandoms;
    }

    public int getNumberOfRandoms() {
        return numberOfRandoms;
    }

    public void setNumberOfRandoms(int numberOfRandoms) {
        this.numberOfRandoms = numberOfRandoms;
    }
}
