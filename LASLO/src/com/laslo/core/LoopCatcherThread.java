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
import static java.lang.System.out;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.Semaphore;
import org.biojava.nbio.core.sequence.DNASequence;

/**
 *
 * @author David A. Mancilla
 */
public class LoopCatcherThread implements Runnable {

    protected boolean extendedMode;
    protected boolean searchReverse;
    protected String additionalSequence;
    protected int maxLength;
    protected int minLength;
    protected InputSequence inputType;
    protected Iterator<String> patternItr;
    protected CSVWriter writer;
    protected DNASequence dnaElement;
    protected final static Semaphore MUTEX = new Semaphore(1);
    protected static Semaphore SEM;
    private static boolean started = false;
    private CountDownLatch latch;

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
     */
    public LoopCatcherThread(boolean extendedMode, String additionalSequence,
            int maxLength, int minLength, DNASequence dnaElement,
            InputSequence inputType, Iterator<String> patternItr,
            CSVWriter writer, boolean searchReverse) {

        int count = 0;
        this.extendedMode = extendedMode;
        this.additionalSequence = additionalSequence;
        this.maxLength = maxLength;
        this.minLength = minLength;
        this.inputType = inputType;
        this.dnaElement = dnaElement;
        this.patternItr = patternItr;
        this.writer = writer;
        this.searchReverse = searchReverse;

        if (!started) {

            count = OSValidator.getNumberOfCPUCores();

            if (count > 1) {
                count = count - 1;
            }

            out.println("[Using " + count + " CPU cores.]");
            SEM = new Semaphore(count);
            started = true;
        }
    }

    /**
     *
     * @param latch
     */
    public void setLatch(CountDownLatch latch) {
        this.latch = latch;
    }

    /**
     *
     */
    @Override
    public void run() {

        RNAfold fold = new RNAfold();

        // sólo para modo dos
        if (extendedMode) {
            fold = new RNAfold(dnaElement.getRNASequence()
                    .getSequenceAsString());
        }
        // III. Loop level
        while (patternItr.hasNext()) {

            String currentPattern = patternItr.next().trim().toUpperCase();

            // 1. Stem research
            if (extendedMode) {

                try {
                    SEM.acquire();
                    sequenceExtendedResearch(
                            dnaElement, 
                            fold.getStructure(),
                            currentPattern, writer, 
                            false, maxLength, minLength,
                            inputType, additionalSequence);
                } catch (InterruptedException ex) {
                    out.println("ERROR: " + ex.getMessage());
                } finally {
                    SEM.release();
                } 

                if (searchReverse) {
                    try {
                        SEM.acquire();
                        sequenceExtendedResearch(
                                dnaElement,
                                fold.getStructure(), 
                                currentPattern, writer,
                                true, maxLength, minLength, 
                                inputType, additionalSequence);
                    } catch (InterruptedException ex) {
                        out.println("ERROR: " + ex.getMessage());
                    } finally {
                        SEM.release();
                    }
                }
            } else {

                SequenceAnalizer.sequenceResearch(dnaElement, currentPattern, 
                        writer, false, maxLength, minLength, inputType, 
                        additionalSequence);

                if (searchReverse) {
                    SequenceAnalizer.sequenceResearch(dnaElement, currentPattern, 
                        writer, true, maxLength, minLength, inputType, 
                        additionalSequence);
                }

            }
        }

        latch.countDown();

    }

    
}
