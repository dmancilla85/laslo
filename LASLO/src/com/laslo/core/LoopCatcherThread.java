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
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import static com.laslo.core.SequenceAnalizer.*;
import com.tools.OSValidator;
import static java.lang.System.out;
import java.util.Map;
import java.util.concurrent.CountDownLatch;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.features.Qualifier;
import org.biojava.nbio.core.sequence.features.TextFeature;

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
    //private final static Semaphore MUTEX = new Semaphore(1);
    //private static Semaphore SEM;
    //private static boolean started = false;
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

        /*if (!started) {

            count = OSValidator.getNumberOfCPUCores();

            if (count > 1) {
                count = count - 1;
            }

            out.println("[Using " + count + " CPU cores.]");
            //SEM = new Semaphore(count);
            started = true;
        }*/
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

                /*try {
                    SEM.acquire();*/
                    sequenceExtendedResearch(dnaElement, fold.getStructure(),
                            currentPattern, writer, false);
                /*} catch (InterruptedException ex) {
                    out.println("ERROR: " + ex.getMessage());
                } finally {
                    SEM.release();
                }*/

                if (searchReverse) {
                    /*try {
                        SEM.acquire();*/
                        sequenceExtendedResearch(dnaElement,
                                fold.getStructure(), currentPattern, writer,
                                true);
                    /*} catch (InterruptedException ex) {
                        out.println("ERROR: " + ex.getMessage());
                    } finally {
                        SEM.release();
                    }*/
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

    /**
     *
     * @param fastaSeq
     * @param hairpinSeq
     * @param stemLoopPattern
     * @param writer
     * @param invert
     * @return
     */
    public int sequenceExtendedResearch(DNASequence fastaSeq, String hairpinSeq,
            String stemLoopPattern, CSVWriter writer, boolean invert) {
        List<StemLoop> slrList = new ArrayList<>();
        StemLoop slr;
        slr = null;
        int size;
        size = 0;
        int length = this.maxLength, posAux, k = 1;
        int loopPos = 0, loopLength = stemLoopPattern.length();
        boolean isValidHairpin;
        String gene = "", synonym = "", note = "";

        String rnaSequence = fastaSeq.getSequenceAsString().replace('T', 'U');
        int sequenceLength = rnaSequence.length();

        if (invert) {
            stemLoopPattern = reverseIt(stemLoopPattern);
        }

        // Convert the original loop pattern to a regular expression
        String regExp = toRegularExpression(stemLoopPattern);
        Pattern p = Pattern.compile(regExp);
        Matcher loopFinder = p.matcher(rnaSequence);

        String rnaLoop = "", rnaSeq = "", hairpinModel = ""; //NOI18N

        // As exists loop matches
        while (loopFinder.find()) {
            loopPos = loopFinder.start();
            isValidHairpin = true;
            rnaLoop = rnaSequence.substring(loopPos, loopFinder.end());
            slr = new StemLoop(this.inputType);

            if (this.inputType != InputSequence.GENBANK) {
                slr.setTags(fastaSeq.getOriginalHeader());
            } else {

                Map qual = ((TextFeature) fastaSeq.getFeaturesByType("gene")
                        .toArray()[0]).getQualifiers();

                if (!qual.isEmpty()) {
                    gene = ((Qualifier) ((ArrayList) (qual.get("gene"))).get(0))
                            .getValue();
                    synonym = ((Qualifier) ((ArrayList) (qual.get("gene_synonym"))).get(0))
                            .getValue();
                    note = ((Qualifier) ((ArrayList) (qual.get("note"))).get(0))
                            .getValue();
                }

                slr.setTags(gene, synonym, fastaSeq.getAccession().getID(),
                        note,
                        ((TextFeature) fastaSeq.getFeaturesByType("CDS")
                                .toArray()[0]).getSource());
            }

            try {
                if ((loopPos - length) > 0
                        && (loopPos + loopLength + length) < sequenceLength) {

                    rnaSeq = rnaSequence.substring(loopPos - length,
                            loopPos + loopLength + length);
                    hairpinModel = hairpinSeq.substring(loopPos - length,
                            loopPos + loopLength + length);

                    isValidHairpin = isRNAPair(rnaSequence.charAt(loopPos - 1),
                            rnaSequence.charAt(loopPos + rnaLoop.length()));

                    if (isValidHairpin) {
                        hairpinModel = SequenceAnalizer.isValidHairpin(
                                maxLength, minLength,hairpinModel, loopLength,
                                loopPos, rnaSeq);
                        isValidHairpin = hairpinModel.length() > 0;
                    }

                } else {
                    isValidHairpin = false;
                }

            } catch (Exception e) {
                out.println("ERROR: " + e.getMessage());
            }

            if (isValidHairpin) {
                int extIzq = hairpinModel.lastIndexOf("(") + 1; //NOI18N
                int extDer = hairpinModel.length() - extIzq - loopLength;
                rnaSeq = rnaSequence.substring(loopPos - extIzq, loopPos
                        + loopLength + extDer);

                posAux = rnaSequence.indexOf(rnaSeq);

                // Fill the fields
                try {
                    slr.setAdditional5Seq(rnaSequence
                            .substring(posAux - k, posAux));
                    slr.setAdditional3Seq(rnaSequence.substring(posAux
                            + rnaSeq.length(), posAux + rnaSeq.length() + k));
                } catch (IndexOutOfBoundsException e) {
                    slr.setAdditional3Seq("");
                    slr.setAdditional5Seq("");
                }

                slr.setReverse(invert);

                if (this.inputType == InputSequence.GENBANK) {
                    slr.setLocation(loopPos - extIzq);
                }

                slr.setRnaHairpinSequence(rnaSeq);
                slr.setLoop(rnaLoop);
                slr.setStartsAt(loopPos - extIzq);
                slr.setStructure(hairpinModel);
                slr.setSequenceLength(sequenceLength);
                slr.checkPairments();
                slr.setMfe(new RNAfold(rnaSeq).getMfe());
                slr.setNLoop(extIzq);
                slr.setPercent_AG();

                if (!invert) {
                    slr.setLoopPattern(stemLoopPattern);
                } else {
                    slr.setLoopPattern(reverseIt(stemLoopPattern));
                }

                slr.setEndsAt(loopFinder.end() + extDer);
                slr.setPercA_sequence(
                        (rnaSequence.length() - rnaSequence.replace("A", "")
                        .length()) / (float) rnaSequence.length());
                slr.setPercG_sequence(
                        (rnaSequence.length() - rnaSequence.replace("G", "")
                        .length()) / (float) rnaSequence.length());
                slr.setPercC_sequence(
                        (rnaSequence.length() - rnaSequence.replace("C", "")
                        .length()) / (float) rnaSequence.length());
                slr.setPercU_sequence(
                        (rnaSequence.length() - rnaSequence.replace("U", "")
                        .length()) / (float) rnaSequence.length());

                if (this.additionalSequence.length() > 0) {
                    slr.setAdditionalSeqLocations(
                            getPatternLocations(rnaSequence,
                                    this.additionalSequence));
                }

                slr.setRelativePos((double) slr.getStartsAt()
                        / (double) rnaSequence.length());

                slrList.add(slr);
            }
        }

        slr = null;
        Iterator<StemLoop> itr = slrList.iterator();

        while (itr.hasNext()) {
            StemLoop element = itr.next();
            /*try {
                MUTEX.acquire();*/
                writer.writeNext(element.toRowCSV().split(";")); //NOI18N
            /*} catch (InterruptedException ex) {
                out.println("ERROR: " + ex.getMessage());
            } finally {
                MUTEX.release();
            }*/
        }

        size = slrList.size();
        slrList.clear();
        slrList = null;

        return size;
    }
}
