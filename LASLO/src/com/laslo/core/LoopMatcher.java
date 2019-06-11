/*
 * Copyright (C) 2018 David A. Mancilla
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
package com.laslo.core;

import com.opencsv.CSVWriter;
import com.tools.UShuffle;
import com.tools.io.FASTACorrector;
import com.tools.io.InputSequence;
import com.tools.io.SourceFile;
import static com.tools.io.SourceFile.*;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import static java.lang.Math.round;
import static java.lang.System.err;
import static java.lang.System.out;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import org.biojava.nbio.core.sequence.DNASequence;
import java.util.ResourceBundle;
import static org.biojava.nbio.core.sequence.io.FastaReaderHelper.readFastaDNASequence;
import static org.biojava.nbio.core.sequence.io.GenbankReaderHelper.readGenbankDNASequence;
import javax.swing.JProgressBar;

/**
 * @author David
 *
 */
public class LoopMatcher {

    private String pathOut;
    private String pathIn;
    private ArrayList<String> loopPatterns;
    private InputSequence inputType;
    private int minLength;
    private int maxLength;
    private int maxWooble;
    private ResourceBundle bundle;
    private int maxMismatch;
    private File[] fileList;
    private boolean extendedMode;
    private boolean makeRandoms;
    private int numberOfRandoms;
    private int kLetRandoms;
    private File actualFile;
    private String additionalSequence;
    private boolean searchReverse;
    private int progress;
    private JProgressBar jpBar;

    /**
     *
     * @param pathOut
     * @param pathIn
     * @param loopPatterns
     * @param additionalSequence
     * @param inputType
     * @param minLength
     * @param maxLength
     * @param maxWooble
     * @param maxMismatch
     * @param locale
     * @param kLetRandoms
     * @param searchReverse
     */
    public LoopMatcher(String pathOut, String pathIn,
            ArrayList<String> loopPatterns, String additionalSequence,
            InputSequence inputType, int minLength, int maxLength,
            int maxWooble, int maxMismatch, Locale locale, int kLetRandoms,
            boolean searchReverse) {
        this.pathOut = pathOut;
        this.loopPatterns = loopPatterns;
        this.inputType = inputType;
        this.minLength = minLength;
        this.maxLength = maxLength;
        this.maxWooble = maxWooble;
        this.maxMismatch = maxMismatch;
        this.fileList = null;
        this.extendedMode = false;
        this.makeRandoms = false;
        this.progress = 0;
        this.numberOfRandoms = 0;
        this.additionalSequence = additionalSequence;
        this.kLetRandoms = kLetRandoms;
        this.bundle = ResourceBundle.getBundle("resources/Bundle", locale);
        SequenceAnalizer.setBundle(bundle);
        this.searchReverse = searchReverse;
    }

    /**
     * Constructor
     */
    public LoopMatcher() {
        this("", "", new ArrayList<>(), BiologicPatterns.PUM1,
                InputSequence.ENSEMBL, //NOI18N
                4, 16, 2, 0, new Locale("es", "ES"), 2, false);
        SequenceAnalizer.setBundle(bundle);
    }

    /**
     * 
     * @return 
     */
    public boolean isExtendedMode() {
        return extendedMode;
    }

    public void setProgressBar(JProgressBar bar){
        this.jpBar = bar;
    }
    
    /**
     * 
     * @param extendedMode 
     */
    public void setExtendedMode(boolean extendedMode) {
        this.extendedMode = extendedMode;
    }

    /**
     * 
     * @return 
     */
    public boolean isSearchReverse() {
        return searchReverse;
    }

    /**
     * 
     * @param searchReverse 
     */
    public void setSearchReverse(boolean searchReverse) {
        this.searchReverse = searchReverse;
    }
    
    public int getProgress(){
        return this.progress;
    }

    /**
     * 
     * @return 
     */
    public int getkLetRandoms() {
        return kLetRandoms;
    }

    /**
     * 
     * @param kLetRandoms 
     */
    public void setkLetRandoms(int kLetRandoms) {
        this.kLetRandoms = kLetRandoms;
    }

    /**
     * 
     * @return 
     */
    public String getAdditionalSequence() {
        return additionalSequence;
    }

    /**
     * 
     * @param additionalSequence 
     */
    public void setAdditionalSequence(String additionalSequence) {
        this.additionalSequence = additionalSequence;
    }

    /**
     * 
     * @return 
     */
    @SuppressWarnings("ReturnOfCollectionOrArrayField")
    public File[] getFileList() {
        return fileList;
    }

    /**
     *
     * @param mode
     */
    public void setIsExtendedMode(boolean mode) {
        this.setExtendedMode(mode);
    }

    /**
     *
     * @return
     */
    public boolean getIsExtendedMode() {
        return this.isExtendedMode();
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
    @SuppressWarnings("ReturnOfCollectionOrArrayField")
    public ArrayList<String> getLoopPatterns() {
        return loopPatterns;
    }

    /**
     *
     * @param loopPatterns
     */
    @SuppressWarnings("AssignmentToCollectionOrArrayFieldFromParameter")
    public void setLoopPatterns(ArrayList<String> loopPatterns) {
        this.loopPatterns = loopPatterns;
    }

    /**
     *
     * @param bundle
     */
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

    /**
     *
     * @param maxLength
     */
    public void setMaxLength(int maxLength) {
        this.maxLength = maxLength;
    }

    /**
     *
     * @return
     */
    public int getMaxWooble() {
        return maxWooble;
    }

    /**
     *
     * @param maxWooble
     */
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

            for (StemLoop stemResearch1 : stemResearch) {
                String[] data;
                data = stemResearch1.toRowCSV().split(";");
                writer.writeNext(data);
            }

            writer.close();

        } catch (IOException e) {
            // TODO Auto-generated catch block
            err.println(java.text.MessageFormat.format(
                            getBundle()
                                    .getString("ERROR_EX"), new Object[] {e.getMessage()}));
        }
    }

    /**
     * Union two lists of files
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
    public boolean startReadingFiles() {
        // To check the elapsed time
        Calendar ini, fin;
        boolean isGenBank = false;
        ini = Calendar.getInstance();
        out.println(java.text.MessageFormat.format(getBundle().getString("START_TIME"),
                new Object[]{Calendar.getInstance().getTime()}));
        out.flush();

        if (getFileList() == null) {
            return false;
        }

        if (this.isMakeRandoms()) {
            out.print(getBundle().getString("MAKING_RANDOM_SEQUENCES"));
            for (File currentFile : getFileList()) {
                if (currentFile.isFile()
                        && (currentFile.toString().endsWith(getFASTA_EXT())
                        || currentFile.toString().endsWith(getFASTA_EXT_2())
                        || currentFile.toString().endsWith(getGENBANK_EXT()))) {

                    //fileName = currentFile.getName();
                    LinkedHashMap<String, DNASequence> dnaFile = null;

                    try {
                        if (currentFile.toString().endsWith(getGENBANK_EXT())) {
                            dnaFile = readGenbankDNASequence(currentFile, false);
                            isGenBank = true;
                        } else {
                            dnaFile = readFastaDNASequence(currentFile, false);
                        }

                    } catch (IOException ex) {
                        err.println(java.text.MessageFormat.format(
                            getBundle()
                                    .getString("ERROR_EX"), new Object[] {ex.getMessage()}));
                    } catch (Exception ex) {
                        err.println(java.text.MessageFormat.format(
                            getBundle()
                                    .getString("ERROR_EX"), new Object[] {ex.getMessage()}));
                    }
                    UShuffle.makeShuffleSequences(getPathOut(), currentFile.getName(),
                            dnaFile, getNumberOfRandoms(), getkLetRandoms(), isGenBank);
                }
            }
            out.print(getBundle().getString("DONE"));

            File folder;
            folder = new File(getPathIn() + UShuffle.getRandomDir());
            File[] randomFiles;
            randomFiles = folder.listFiles();
            setFileList(unionFiles(getFileList(), randomFiles));
        }

        // I. File level (hacer un hilo)
        for (File currentFile : getFileList()) {
            if (currentFile.isFile() && 
                    (currentFile.toString().endsWith(getGENBANK_EXT()) || 
                    currentFile.toString().endsWith(getFASTA_EXT()) ||
                    currentFile.toString().endsWith(getFASTA_EXT_2()))) {

                this.setActualFile(currentFile);
                callProcessThreads();
            }
        }

        fin = Calendar.getInstance();
        out.print(java.text.MessageFormat.format(getBundle()
                .getString("TOTAL_TIME"),
                new Object[]{((fin.getTimeInMillis() 
                        - ini.getTimeInMillis()) / 1000)}));

        out.flush();

        // Memory cleaning
        setFileList(null);
        return true;
    }

    /**
     * Process the files selected
     */
    @SuppressWarnings({"ValueOfIncrementOrDecrementUsed", "empty-statement"})
    public void callProcessThreads() {

        CSVWriter writer;
        boolean genbank;
        String fileName, fileOut;
        final int MAX_HILOS = 10;
        int totalSecuencias;
        int nHilos = MAX_HILOS;
        int i;
        int count;
        count = 0;
        ExecutorService pool;
        CountDownLatch latch;
        Calendar ini, fin;
        String mode = "local";
        
        if(this.extendedMode){
            mode = "global";
        }

        try {
            fileName = getActualFile().getName();
            //genbank = false;
            out.print(java.text.MessageFormat.format(getBundle()
                    .getString("FILE_PRINT"), new Object[]{fileName})); 
            out.flush();
            fileName = fileName.replaceFirst("[.][^.]+$", ""); 
            fileOut = this.getPathOut() + "\\" + fileName + "." +
                    mode + "." + this.minLength + "." + this.maxLength +
                    getCSV_EXT();

            if (new File(fileOut).exists()) {
                try {
                    (new File(fileOut)).delete();

                } catch (Exception io) {
                    err.println(io.getMessage());
                    err.flush();
                    return;
                }
            }

            // Generation of the iterator of {id,sequence}
            LinkedHashMap<String, DNASequence> fasta;

            if (getActualFile().getName().endsWith(getGENBANK_EXT())) {
                fasta = readGenbankDNASequence(getActualFile());
                genbank = true;
            } else {
                fasta = readFastaDNASequence(getActualFile(), false);
                genbank = false;
            }

            if (fasta.isEmpty()) {
                err.println(getBundle().getString("INVALID_FILE_FORMAT"));
                err.println(getBundle().getString("TRYING_TO_FIX"));
                boolean formatFile = FASTACorrector
                        .formatFile(getActualFile().getAbsolutePath());
                if (formatFile) {
                    fasta = readFastaDNASequence(getActualFile(), false);
                } else {
                    err.println(getBundle().getString("CANT_PROCESS"));
                    return;
                }
            }

            if (!genbank) {
                this.setInputType(SourceFile.detectHeader(
                        fasta.entrySet().iterator().next().getValue()
                                .getOriginalHeader()));
            } else {
                this.setInputType(InputSequence.GENBANK);
            }

            int listSize = fasta.size();

            writer = new CSVWriter(new FileWriter(fileOut), ';',
                    CSVWriter.DEFAULT_QUOTE_CHARACTER,
                    CSVWriter.DEFAULT_ESCAPE_CHARACTER,
                    CSVWriter.DEFAULT_LINE_END);
            writer.writeNext(StemLoop.getHeader(this.getInputType())
                    .split(";")); 

            // II. Transcript level
            if (listSize < nHilos) {
                nHilos = listSize;
            }

            pool = Executors.newFixedThreadPool(listSize);
            latch = new CountDownLatch(nHilos);
            i = 1;

            ini = Calendar.getInstance();
            totalSecuencias = fasta.entrySet().size();
            this.progress = 0;
                    
            for (Map.Entry<String, DNASequence> entry : fasta.entrySet()) {

                DNASequence element = entry.getValue();
                Iterator<String> patternItr = getLoopPatterns().iterator();
                count++;
                
                LoopMatcherThread thread = new LoopMatcherThread(
                        isExtendedMode(), getAdditionalSequence(), 
                        getMaxLength(), getMinLength(), element, 
                        getInputType(), patternItr, writer, isSearchReverse(),
                        bundle);

                this.progress = (int)round(count/(double)totalSecuencias * 100);
                jpBar.setValue(progress);
                
                if (i++ <= nHilos) {
                    thread.setLatch(latch);
                    pool.execute(thread);
                } else {
                    i = 1;
                    latch.await();
                    pool.shutdown();
                    
                    if (fasta.size() - count < nHilos) {
                        nHilos = fasta.size() - count;
                    }
                    
                    if (nHilos > 0) {
                        pool = Executors.newFixedThreadPool(nHilos);
                        latch = new CountDownLatch(nHilos);
                    }
                }
                
                // Print
                /*if(this.progress % 15 == 0 || this.progress % 6 == 0){
                    out.printf(getBundle()
                            .getString("PROGRESO"), 
                            this.progress);
                }*/
                
            }
            
            jpBar.setValue(100);

            if (latch.getCount() > 0) {
                latch.await();
                pool.shutdown();
            }
            
            if(!pool.isShutdown()){
                pool.shutdown();
            }
            
            while(!pool.isTerminated());

            writer.close();
            fasta.clear();

            fin = Calendar.getInstance();
            //getBundle().getString("SUMMARY")
            out.printf(" = {%d secuencias} - Tiempo: %3.2fs.\n", 
                    totalSecuencias, 
                    (fin.getTimeInMillis() - ini.getTimeInMillis()) / (double)1000);
            
            /*out.print(java.text.MessageFormat
                    .format(getBundle()
                            .getString("SEQUENCES"), new Object[] {totalSecuencias}));
            out.print(java.text.MessageFormat.format(getBundle().getString("TIME"), 
                    new Object[] {(fin.getTimeInMillis() - ini.getTimeInMillis()) / 1000}));
            out.println();*/

        } catch (FileNotFoundException ex) {
            err.println(getBundle().getString("CANT_OPEN_FILE"));
            err.println(java.text.MessageFormat.format(
                            getBundle()
                                    .getString("ERROR_EX"), new Object[] {ex.getMessage()}));
        } catch (Exception ex) {
            err.println(java.text.MessageFormat.format(
                            getBundle()
                                    .getString("ERROR_EX"), new Object[] {ex.getMessage()}));
        } 
    }

    /**
     *
     * @return
     */
    public boolean isMakeRandoms() {
        return makeRandoms;
    }

    /**
     *
     * @param makeRandoms
     */
    public void setMakeRandoms(boolean makeRandoms) {
        this.makeRandoms = makeRandoms;
    }

    /**
     *
     * @return
     */
    public int getNumberOfRandoms() {
        return numberOfRandoms;
    }

    /**
     *
     * @param numberOfRandoms
     */
    public void setNumberOfRandoms(int numberOfRandoms) {
        this.numberOfRandoms = numberOfRandoms;
    }

    /**
     * @return the actualFile
     */
    public File getActualFile() {
        return actualFile;
    }

    /**
     * @return the bundle
     */
    public ResourceBundle getBundle() {
        return bundle;
    }

    /**
     * @param actualFile the actualFile to set
     */
    public void setActualFile(File actualFile) {
        this.actualFile = actualFile;
    }

    /**
     * @param fileList the fileList to set
     */
    @SuppressWarnings("AssignmentToCollectionOrArrayFieldFromParameter")
    public void setFileList(File[] fileList) {
        this.fileList = fileList;
    }
}
