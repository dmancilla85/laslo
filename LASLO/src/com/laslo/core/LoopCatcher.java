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
import java.util.ResourceBundle;
import static java.util.ResourceBundle.getBundle;
import java.util.Set;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
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
    protected boolean searchReverse;

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
    public LoopCatcher(String pathOut, String pathIn,
            ArrayList<String> loopPatterns, String additionalSequence,
            InputSequence inputType,
            int minLength, int maxLength,
            int maxWooble, int maxMismatch,
            Locale locale, int kLetRandoms,
            boolean searchReverse) {
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
        this.searchReverse = searchReverse;
    }

    /**
     *
     */
    public LoopCatcher() {
        this("", "", new ArrayList<>(), BiologyPatterns.PUM1,
                InputSequence.ENSEMBL, //NOI18N
                4, 16, 2, 0, new Locale("es", "AR"), 2, false);
    }
    
    public boolean isExtendedMode() {
        return extendedMode;
    }
    
    public void setExtendedMode(boolean extendedMode) {
        this.extendedMode = extendedMode;
    }
    
    public boolean isSearchReverse() {
        return searchReverse;
    }
    
    public void setSearchReverse(boolean searchReverse) {
        this.searchReverse = searchReverse;
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
            
            for (int i = 0; i < stemResearch.size(); i++) {
                String[] data;
                data = stemResearch.get(i).toRowCSV().split(";"); //NOI18N
                writer.writeNext(data);
                //allData.add(data);
            }
            
            writer.close();
            
        } catch (IOException e) {
            // TODO Auto-generated catch block
            out.println("ERROR: " + e.getMessage());
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
        ini = Calendar.getInstance();
        out.println(
                java.text.MessageFormat.format(bundle.getString("START_TIME"),
                        new Object[]{Calendar.getInstance().getTime()}));
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
                        out.println("ERROR: " + ex.getMessage());
                    } catch (Exception ex) {
                        out.println("ERROR: " + ex.getMessage());
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
                callProcessThreads();
            }
        }
        
        fin = Calendar.getInstance();
        out.print(java.text.MessageFormat.format(bundle.getString("TOTAL_TIME"),
                new Object[]{((fin.getTimeInMillis() - ini.getTimeInMillis()) / 1000)}) + " s.");
        
        out.flush();

        // Memory cleaning
        fileList = null;
        return true;
    }

    /**
     * Process the files selected
     */
    public void callProcessThreads() {
        
        CSVWriter writer;
        boolean genbank;
        String fileName, fileOut;
        final int MAX_HILOS = 10;
        int secuencias, nHilos = MAX_HILOS, i, count = 0;
        ExecutorService pool;
        CountDownLatch latch;
        Calendar ini, fin;
        
        try {
            fileName = actualFile.getName();
            genbank = false;
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
                genbank = true;
            } else {
                fasta = readFastaDNASequence(actualFile,
                        actualFile.length() > (52428800));
                genbank = false;
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
            
            if (!genbank) {
                this.inputType = SourceFile.detectHeader(fasta.entrySet().iterator()
                        .next().getValue().getOriginalHeader());
            } else {
                this.inputType = InputSequence.GENBANK;
            }
            
            int listSize = fasta.size();
            
            writer = new CSVWriter(new FileWriter(fileOut), ';',
                    CSVWriter.DEFAULT_QUOTE_CHARACTER,
                    CSVWriter.DEFAULT_ESCAPE_CHARACTER,
                    CSVWriter.DEFAULT_LINE_END);
            writer.writeNext(StemLoop.getHeader(this.inputType).split(";")); //NOI18N

            // II. Transcript level
            if (listSize < nHilos) {
                nHilos = listSize;
            }
            
            pool = Executors.newFixedThreadPool(listSize);
            latch = new CountDownLatch(listSize);
            i = 0;
            
            ini = Calendar.getInstance();
            secuencias = 0;
            
            for (Map.Entry<String, DNASequence> entry : fasta.entrySet()) {
                
                DNASequence element = entry.getValue();
                Iterator<String> patternItr = loopPatterns.iterator();
                count++;
                secuencias++;
                LoopCatcherThread thread = new LoopCatcherThread(extendedMode,
                        additionalSequence, maxLength, minLength, element,
                        inputType, patternItr, writer, searchReverse);
                
                thread.setLatch(latch);
                pool.execute(thread);
                
            }
            
            if (latch.getCount() > 0) {
                latch.await();
                pool.shutdown();
            }
            
            pool.shutdown();
            
            writer.close();
            writer = null;
            fasta.clear();
            fasta = null;
            
            out.print(" Secuencias: " + secuencias);
            fin = Calendar.getInstance();
            out.print(" Tiempo: " + (fin.getTimeInMillis()
                    - ini.getTimeInMillis()) / 1000 + " s.");
            out.println();
            
        } catch (FileNotFoundException ex) {
            /*Logger.getLogger(LoopCatcher.class.getName())
                    .log(Level.SEVERE, null, ex);*/
            out.println(bundle.getString("CANT_OPEN_FILE"));
        } catch (IOException ex) {
            out.println("ERROR: " + ex.getMessage());
        } catch (Exception ex) {
            out.println("ERROR: " + ex.getMessage());
        }
        
        writer = null;
        fileName = null;
        fileOut = null;
        System.gc();
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
}
