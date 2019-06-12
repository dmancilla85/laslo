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
import com.tools.RNAfold;
import com.tools.io.InputSequence;
import java.util.Iterator;
import static com.laslo.core.SequenceAnalizer.*;
import com.tools.OSValidator;
import static java.lang.System.err;
import static java.lang.System.out;
import java.util.ResourceBundle;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.Semaphore;
import org.biojava.nbio.core.sequence.DNASequence;

/**
 *
 * @author David A. Mancilla
 */
public class LoopMatcherThread implements Runnable {

    private boolean extendedMode;
    private boolean searchReverse;
    private String additionalSequence;
    private int maxLength;
    private int minLength;
    private InputSequence inputType;
    private Iterator<String> patternItr;
    private CSVWriter writer;
    private DNASequence dnaElement;
    private static Semaphore MUTEX = new Semaphore(1);
    private static Semaphore SEM;
    private static boolean started = false;
    private CountDownLatch latch;
    private ResourceBundle bundle;
    
    /**
     * 
     * @param extendedMode
     * @param additionalSequence
     * @param maxLength
     * @param minLength
     * @param dnaElement
     * @param inputType
     * @param patternItr
     * @param writer
     * @param searchReverse 
     * @param bundle
     */
    public LoopMatcherThread(boolean extendedMode, String additionalSequence,
            int maxLength, int minLength, DNASequence dnaElement,
            InputSequence inputType, Iterator<String> patternItr,
            CSVWriter writer, boolean searchReverse, ResourceBundle bundle) {

        int countThreads;
        this.extendedMode = extendedMode;
        this.additionalSequence = additionalSequence;
        this.maxLength = maxLength;
        this.minLength = minLength;
        this.inputType = inputType;
        this.dnaElement = dnaElement;
        this.patternItr = patternItr;
        this.writer = writer;
        this.searchReverse = searchReverse;
        this.bundle = bundle;

        if (!started) {

            countThreads = OSValidator.getNumberOfCPUCores();
            countThreads-=1;
            
            /*if (countThreads > 1) {
                countThreads -= 1;
            }*/

            out.println(java.text.MessageFormat.format(bundle
                    .getString("USING_N_CORES"), new Object[] {countThreads}));
            SEM = new Semaphore(countThreads);
            started = true;
        }
    }

    /**
     * @return the bundle
     */
    public ResourceBundle getBundle() {
        return bundle;
    }
    
    public void setBundle(ResourceBundle bundle){
        this.bundle = bundle;
    }
    
    /**
     *
     */
    @Override
    public void run() {

        RNAfold fold = new RNAfold();

        // sólo para modo dos
        if (isExtendedMode()) {
            fold = new RNAfold(getDnaElement().getRNASequence()
                    .getSequenceAsString());
        }
        // III. Loop level
        while (getPatternItr().hasNext()) {

            String currentPattern = getPatternItr().next().trim().toUpperCase();

            // 1. Stem research
            if (isExtendedMode()) {
                
                try {
                    getSEM().acquire();
                    sequenceExtendedResearch(getDnaElement(), 
                            fold.getStructure(),
                            currentPattern, getWriter(), false, getMaxLength(), 
                            getMinLength(), getInputType(), 
                            getAdditionalSequence());
                } catch (InterruptedException ex) {
                    err.println(java.text.MessageFormat.format(
                            getBundle()
                                    .getString("ERROR_EX"), new Object[] {ex.getMessage()}));
                } finally {
                    getSEM().release();
                }
                
                if (isSearchReverse()) {
                    try {
                        getSEM().acquire();
                        sequenceExtendedResearch(getDnaElement(),
                                fold.getStructure(), 
                                currentPattern, getWriter(), true, getMaxLength(), 
                                getMinLength(), getInputType(), 
                                getAdditionalSequence());
                        
                    } catch (InterruptedException ex) {
                        err.println(java.text.MessageFormat.format(
                            getBundle()
                                    .getString("ERROR_EX"), new Object[] {ex.getMessage()}));
                    } finally {
                        getSEM().release();
                    }
                }
            } else {
                
                SequenceAnalizer.sequenceResearch(getDnaElement(), currentPattern, 
                        getWriter(), false, getMaxLength(), getMinLength(), 
                        getInputType(), getAdditionalSequence());
                
                if (isSearchReverse()) {
                    SequenceAnalizer.sequenceResearch(getDnaElement(), 
                            currentPattern, getWriter(), true, getMaxLength(), 
                            getMinLength(), getInputType(), 
                            getAdditionalSequence());
                }
                
            }
        }

        getLatch().countDown();
    }

    /**
     * @return the additionalSequence
     */
    public String getAdditionalSequence() {
        return additionalSequence;
    }

    /**
     * @return the dnaElement
     */
    public DNASequence getDnaElement() {
        return dnaElement;
    }

    /**
     * @return the inputType
     */
    public InputSequence getInputType() {
        return inputType;
    }

    /**
     * @return the latch
     */
    public CountDownLatch getLatch() {
        return latch;
    }

    /**
     * @return the maxLength
     */
    public int getMaxLength() {
        return maxLength;
    }

    /**
     * @return the minLength
     */
    public int getMinLength() {
        return minLength;
    }

    /**
     * @return the patternItr
     */
    public Iterator<String> getPatternItr() {
        return patternItr;
    }

    /**
     * @return the writer
     */
    public CSVWriter getWriter() {
        return writer;
    }

    /**
     * @return the extendedMode
     */
    public boolean isExtendedMode() {
        return extendedMode;
    }

    /**
     * @return the searchReverse
     */
    public boolean isSearchReverse() {
        return searchReverse;
    }

    /**
     * @param additionalSequence the additionalSequence to set
     */
    public void setAdditionalSequence(String additionalSequence) {
        this.additionalSequence = additionalSequence;
    }

    /**
     * @param dnaElement the dnaElement to set
     */
    public void setDnaElement(DNASequence dnaElement) {
        this.dnaElement = dnaElement;
    }

    /**
     * @param extendedMode the extendedMode to set
     */
    public void setExtendedMode(boolean extendedMode) {
        this.extendedMode = extendedMode;
    }

    /**
     * @param inputType the inputType to set
     */
    public void setInputType(InputSequence inputType) {
        this.inputType = inputType;
    }

    /**
     * @param maxLength the maxLength to set
     */
    public void setMaxLength(int maxLength) {
        this.maxLength = maxLength;
    }

    /**
     * @param minLength the minLength to set
     */
    public void setMinLength(int minLength) {
        this.minLength = minLength;
    }

    /**
     * @param patternItr the patternItr to set
     */
    public void setPatternItr(Iterator<String> patternItr) {
        this.patternItr = patternItr;
    }

    /**
     * @param searchReverse the searchReverse to set
     */
    public void setSearchReverse(boolean searchReverse) {
        this.searchReverse = searchReverse;
    }

    /**
     * @param writer the writer to set
     */
    public void setWriter(CSVWriter writer) {
        this.writer = writer;
    }
    
    /**
     * 
     * @return 
     */
    public static Semaphore getMUTEX() {
        return MUTEX;
    }

    /**
     * @return the SEM
     */
    public static Semaphore getSEM() {
        return SEM;
    }

    /**
     * @return the started
     */
    public static boolean isStarted() {
        return started;
    }

    /**
     * @param aMUTEX the MUTEX to set
     */
    public static void setMUTEX(Semaphore aMUTEX) {
        MUTEX = aMUTEX;
    }

    /**
     * @param aSEM the SEM to set
     */
    public static void setSEM(Semaphore aSEM) {
        SEM = aSEM;
    }

    /**
     * @param aStarted the started to set
     */
    public static void setStarted(boolean aStarted) {
        started = aStarted;
    }
    
    /**
     *
     * @param latch
     */
    public void setLatch(CountDownLatch latch) {
        this.latch = latch;
    }
}
