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
import static com.laslo.core.PairmentAnalizer.*;
import com.tools.OSValidator;
import static java.lang.System.out;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.Semaphore;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.apache.commons.lang.StringUtils;
import org.biojava.nbio.core.sequence.DNASequence;
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
    protected DNASequence dnaElement;
    protected InputSequence inputType;
    protected Iterator<String> patternItr;
    protected CSVWriter writer;
    private final static Semaphore MUTEX = new Semaphore(1);
    private static Semaphore SEM;
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
        this.dnaElement = dnaElement;
        this.inputType = inputType;
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
        if (extendedMode) 
            fold = new RNAfold(dnaElement.getRNASequence()
                    .getSequenceAsString());

        // III. Loop level
        while (patternItr.hasNext()) {

            String currentPattern = patternItr.next().trim().toUpperCase();

            // 1. Stem research
            if (extendedMode) {

                try {
                    SEM.acquire();
                    sequenceExtendedResearch(dnaElement, fold.getStructure(),
                            currentPattern, writer, false);
                } catch (InterruptedException ex) {
                    out.println("ERROR: " + ex.getMessage());
                } finally {
                    SEM.release();
                }

                if (searchReverse) {
                    try {
                        SEM.acquire();
                        sequenceExtendedResearch(dnaElement,
                                fold.getStructure(), currentPattern, writer,
                                true);
                    } catch (InterruptedException ex) {
                        out.println("ERROR: " + ex.getMessage());
                    } finally {
                        SEM.release();
                    }
                }

            } else {
                sequenceResearch(dnaElement, currentPattern, writer, false);

                if (searchReverse) {
                    sequenceResearch(dnaElement, currentPattern, writer, true);
                }
            }
        }

        latch.countDown();

    }

    /**
     *
     * @param sequence
     * @param pattern
     * @return
     */
    public List<Integer> getPatternLocations(String sequence, String pattern) {
        List<Integer> locations = new ArrayList<>();
        String regExp = toRegularExpression(pattern);

        Pattern p = Pattern.compile(regExp);
        Matcher loopFinder = p.matcher(sequence);

        while (loopFinder.find()) {
            locations.add(loopFinder.start());
        }

        return locations;
    }

    /**
     *
     * @param fastaSeq
     * @param stemLoopPattern
     * @param writer
     * @param invert
     * @return
     */
    public int sequenceResearch(DNASequence fastaSeq, String stemLoopPattern,
            CSVWriter writer, boolean invert) {

        List<StemLoop> slrList = new ArrayList<>();
        StemLoop slr;
        slr = null;
        int size, posAux, k = 1;
        size = 0;

        int loopPos = 0, loopLength = stemLoopPattern.length();
        boolean isValidHairpin;
        String rnaSequence = fastaSeq.getRNASequence().getSequenceAsString();
        final int sequenceLength = rnaSequence.length();

        if (invert) {
            stemLoopPattern = reverseIt(stemLoopPattern);
        }

        RNAfold fold = new RNAfold();

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

            int length = this.maxLength;
            slr = new StemLoop(this.inputType);

            if (this.inputType != InputSequence.GENBANK) {
                slr.setTags(fastaSeq.getOriginalHeader());
            } else {
                slr.setTags(fastaSeq.getAccession().getID(),
                        fastaSeq.getAccession().getVersion().toString(),
                        fastaSeq.getDescription(),
                        ((TextFeature) fastaSeq.getFeaturesByType("CDS")
                                .toArray()[0]).getSource());
            }

            try {
                if ((loopPos - length) > 0
                        && (loopPos + loopLength + length) < sequenceLength) {
                    rnaSeq = rnaSequence.substring(loopPos - length,
                            loopPos + loopLength + length);

                    isValidHairpin = isComplementaryRNA(rnaSeq
                            .charAt(length - 1), rnaSeq.charAt(length
                            + loopLength))
                            || isComplementaryRNAWooble(rnaSeq
                                    .charAt(length - 1), rnaSeq.charAt(length
                                    + loopLength));

                    if (isValidHairpin) {

                        try {
                            SEM.acquire();
                            fold = new RNAfold(rnaSeq);
                        } catch (InterruptedException ex) {
                            out.println("ERROR: " + ex.getMessage());
                        } finally {
                            SEM.release();
                        }

                        hairpinModel = fold.getStructure();

                        if (rnaSeq.length() != hairpinModel.length()) {
                            out.println("Error NO COINCIDEN:" + rnaSeq + " - "
                                    + hairpinModel);
                        }

                        if (fold.getMfe() == 0.0) {
                            isValidHairpin = false;
                        } else {
                            hairpinModel = isValidHairpin(hairpinModel,
                                    loopLength, loopPos, rnaSeq);
                            isValidHairpin = hairpinModel.length() > 0;
                        }
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
                                    additionalSequence));
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
            try {
                MUTEX.acquire();
                writer.writeNext(element.toRowCSV().split(";")); //NOI18N
            } catch (InterruptedException ex) {
                out.println("ERROR: " + ex.getMessage());
            } finally {
                MUTEX.release();
            }
        }

        size = slrList.size();
        slrList.clear();
        slrList = null;

        return size;
    }

    /**
     *
     * @param hairpin
     * @param loopLength
     * @param loopPos
     * @param seq
     * @return
     */
    @SuppressWarnings("empty-statement")
    public String isValidHairpin(String hairpin, int loopLength, int loopPos,
            String seq) {
        boolean ret, begin = true;
        String extremoIzq = "", extremoDer = ""; //NOI18N
        int stemLength = (hairpin.length() - loopLength) / 2;
        int newLength = stemLength;

        // Pasos
        // 0. Verificar que estèn los 3 simbolos
        ret = hairpin.indexOf('.') >= 0
                && hairpin.indexOf('(') >= 0
                && hairpin.indexOf(')') >= 0;

        // 1. Verificar que en la posición del loop no haya brackets
        ret = ret && hairpin.substring(stemLength, stemLength + loopLength)
                .replaceAll("\\.", "").length() == 0; //NOI18N
        // 1'. Verificar que haya un cierre
        if (ret) {
            ret = hairpin.charAt(stemLength - 1) == '('
                    && hairpin.charAt(stemLength + loopLength) == ')';
        }

        // 2. Analizar lado izq. y derecho del stem desde extremo 5',
        if (ret) {
            // 2.a 	Si encuentra un . o un ) en el ext izquierdo, remover 
            // 		base (stemLength-1)
            // Repetir 2.a hasta que complete el recorrido. 
            extremoIzq = hairpin.substring(0, stemLength);

            try {
                for (int i = 0; i < stemLength; i++) {
                    if (extremoIzq.charAt(i) == '.'
                            || extremoIzq.charAt(i) == ')') {
                        if (begin) {
                            newLength--;
                            begin = true;
                        } else {
                            if (extremoIzq.charAt(i) == ')') {
                                ret = false;
                            }
                        }
                    } else {
                        begin = false;
                    }
                }

                hairpin = hairpin.substring(stemLength - newLength,
                        hairpin.length() - (stemLength - newLength));
                stemLength = newLength;

                // 2.b 	Si encuentra un . o un ( en el ext derecho, remover base 
                //		(stemLength-1)
                // Repetir 2.b hasta que complete el recorrido. 	
                extremoDer = hairpin.substring(stemLength + loopLength,
                        hairpin.length());
                begin = true;

                for (int i = stemLength - 1; i >= 0; i--) {
                    if (extremoDer.charAt(i) == '.'
                            || extremoDer.charAt(i) == '(') {
                        if (begin) {
                            newLength--;
                            begin = true;
                        } else {
                            if (extremoDer.charAt(i) == '(') {
                                ret = false;
                                //newLength--;
                            }
                        }
                    } else {
                        begin = false;
                    }
                }

                hairpin = hairpin.substring(0, stemLength + loopLength);
                extremoDer = extremoDer.substring(0, newLength - 1);
                hairpin = hairpin + extremoDer;

                int auxR = StringUtils.countMatches(hairpin, ")");
                int auxL = StringUtils.countMatches(hairpin, "(");
                int i = hairpin.lastIndexOf("(");
                while (auxL > 0 && auxR > 0) {
                    if (hairpin.charAt(i--) == '(') {
                        auxR--;
                        auxL--;
                    }
                }

                hairpin = hairpin.substring(i + 1, hairpin.length());

                int k;
                for (k = hairpin.length() - 1; hairpin.charAt(k) == '.'; k--);

                hairpin = hairpin.substring(0, k + 1);
            } catch (IndexOutOfBoundsException e) {
                out.println("Error: " + e.getMessage());
            }
        }

        if (ret) {
            ret = StringUtils.countMatches(hairpin, "(") >= minLength
                    && StringUtils.countMatches(hairpin, ")") >= minLength;
        }
        // Si la longitud es menor a la esperada, rechazar.
        if (!ret) {
            return ""; //NOI18N
        } else {
        }
        return hairpin;
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

        String rnaSequence = fastaSeq.getRNASequence().getSequenceAsString();
        final int sequenceLength = rnaSequence.length();

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
                slr.setTags(fastaSeq.getAccession().getID(),
                        fastaSeq.getAccession().getVersion().toString(),
                        fastaSeq.getDescription(),
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

                    isValidHairpin = isRNAPair(rnaSeq.charAt(length - 1),
                            rnaSeq.charAt(length + loopLength));

                    if (isValidHairpin) {
                        hairpinModel = isValidHairpin(hairpinModel, loopLength,
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
            try {
                MUTEX.acquire();
                writer.writeNext(element.toRowCSV().split(";")); //NOI18N
            } catch (InterruptedException ex) {
                out.println("ERROR: " + ex.getMessage());
            } finally {
                MUTEX.release();
            }
        }

        size = slrList.size();
        slrList.clear();
        slrList = null;

        return size;
    }
}
