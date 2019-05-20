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

import static com.laslo.core.LoopMatcherThread.MUTEX;
import static com.laslo.core.LoopMatcherThread.SEM;
import com.opencsv.CSVWriter;
import com.tools.RNAfold;
import com.tools.io.InputSequence;
import static java.lang.System.out;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.apache.commons.lang.StringUtils;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.features.Qualifier;
import org.biojava.nbio.core.sequence.features.TextFeature;

/**
 * @author David
 *
 */
public class SequenceAnalizer {

    private static boolean omitRNAFold = false;

    /**
     *
     * @param fastaSeq
     * @param stemLoopPattern
     * @param writer
     * @param invert
     * @param maxLength
     * @param minLength
     * @param inputType
     * @param additionalSeq
     * @return
     */
    public static synchronized int sequenceResearch(DNASequence fastaSeq, String stemLoopPattern,
            CSVWriter writer, boolean invert, int maxLength, int minLength,
            InputSequence inputType, String additionalSeq) {

        List<StemLoop> slrList = new ArrayList<>();
        StemLoop slr;
        slr = null;
        int size, posAux, k = 1;
        size = 0;
        String gene = "", note = "", synonym = "", id = "", cds = "";
        int loopPos = 0, loopLength = stemLoopPattern.length();
        boolean isValidHairpin;
        String rnaSequence = fastaSeq.getSequenceAsString().replace('T', 'U');

        int sequenceLength = rnaSequence.length();

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

            int length = maxLength;
            slr = new StemLoop(inputType);

            if (inputType != InputSequence.GENBANK) {
                slr.setTags(fastaSeq.getOriginalHeader());
            } else {

                if (!fastaSeq.getOriginalHeader().contains("@")) {
                    Map qual = ((TextFeature) fastaSeq.getFeaturesByType("gene")
                            .toArray()[0]).getQualifiers();

                    if (!qual.isEmpty()) {
                        gene = ((Qualifier) ((ArrayList) (qual.get("gene"))).get(0))
                                .getValue();
                        synonym = ((Qualifier) ((ArrayList) (qual.get("gene_synonym"))).get(0))
                                .getValue();
                        note = ((Qualifier) ((ArrayList) (qual.get("note"))).get(0))
                                .getValue();
                        id = fastaSeq.getAccession().getID();
                        cds = ((TextFeature) fastaSeq.getFeaturesByType("CDS")
                                .toArray()[0]).getSource();
                    }
                } else {
                    String auxH[] = fastaSeq.getOriginalHeader().split("@");

                    gene = auxH[0];
                    synonym = auxH[1];
                    note = auxH[2];
                    id = auxH[3];
                    cds = auxH[4];
                }

                slr.setTags(gene, synonym, id, note, cds);
            }

            try {
                if ((loopPos - length) > 0
                        && (loopPos + loopLength + length) < sequenceLength) {
                    rnaSeq = rnaSequence.substring(loopPos - length,
                            loopPos + loopLength + length);

                    isValidHairpin = isRNAPair(rnaSequence.charAt(loopPos - 1),
                            rnaSequence.charAt(loopPos + rnaLoop.length()));

                    if (isValidHairpin) {
                        isValidHairpin = isRNAPair(rnaSequence.charAt(loopPos - 2),
                                rnaSequence.charAt(loopPos + rnaLoop.length() + 1));
                    }

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

                        if (fold.getMfe() >= 0.0 && !omitRNAFold) {
                            isValidHairpin = false;
                        } else {
                            hairpinModel = isValidHairpin(minLength,
                                    hairpinModel, loopLength, loopPos, rnaSeq);
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

                if (inputType == InputSequence.GENBANK) {
                    slr.setLocation(loopPos - extIzq);
                }

                slr.setRnaHairpinSequence(rnaSeq);
                slr.setLoop(rnaLoop);
                slr.setStartsAt(loopPos - extIzq);
                slr.setStructure(hairpinModel);
                slr.setSequenceLength(sequenceLength);
                slr.checkPairments();
                slr.checkInternalLoops();
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

                if (additionalSeq.trim().length() > 0) {
                    slr.setAdditionalSeqLocations(
                            getPatternLocations(rnaSequence,
                                    additionalSeq));
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

    @SuppressWarnings("empty-statement")
    public static String isValidHairpin(int minLength,
            String hairpin, int loopLength, int loopPos,
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
                        } else if (extremoIzq.charAt(i) == ')') {
                            ret = false;
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
                        } else if (extremoDer.charAt(i) == '(') {
                            ret = false;
                            //newLength--;
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
     * @param maxLength
     * @param minLength
     * @param inputType
     * @param additionalSequence
     * @return
     */
    public static int sequenceExtendedResearch(
            DNASequence fastaSeq,
            String hairpinSeq,
            String stemLoopPattern,
            CSVWriter writer, boolean invert,
            int maxLength, int minLength,
            InputSequence inputType,
            String additionalSequence) {
        List<StemLoop> slrList = new ArrayList<>();
        StemLoop slr;
        int size = -1;
        int length = maxLength, posAux, k = 1;
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

        String rnaLoop, rnaSeq, hairpinModel = ""; //NOI18N

        // As exists loop matches
        while (loopFinder.find()) {
            loopPos = loopFinder.start();
            isValidHairpin = true;
            rnaLoop = rnaSequence.substring(loopPos, loopFinder.end());
            slr = new StemLoop(inputType);

            if (inputType != InputSequence.GENBANK) {
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
                        hairpinModel = isValidHairpin(minLength,
                                hairpinModel, loopLength, loopPos, rnaSeq);
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

                if (inputType == InputSequence.GENBANK) {
                    slr.setLocation(loopPos - extIzq);
                }

                slr.setRnaHairpinSequence(rnaSeq);
                slr.setLoop(rnaLoop);
                slr.setStartsAt(loopPos - extIzq);
                slr.setStructure(hairpinModel);
                slr.setSequenceLength(sequenceLength);
                slr.checkPairments();
                slr.checkInternalLoops();
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

                if (additionalSequence.length() > 0) {
                    slr.setAdditionalSeqLocations(
                            getPatternLocations(rnaSequence, additionalSequence));
                }

                slr.setRelativePos((double) slr.getStartsAt()
                        / (double) rnaSequence.length());

                slrList.add(slr);
            }
        }

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

        return size;
    }

    /**
     *
     * @param sequence
     * @param pattern
     * @return
     */
    public static List<Integer> getPatternLocations(String sequence, String pattern) {
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
     * @param base1
     * @param base2
     * @return
     */
    public static boolean isComplementaryDNA(char base1, char base2) {
        // ok
        boolean isComplement = false;

        switch (base1) {
            case 'A':
                isComplement = (base2 == 'T');
                break;
            case 'T':
                isComplement = (base2 == 'A');
                break;
            case 'G':
                isComplement = (base2 == 'C');
                break;
            case 'C':
                isComplement = (base2 == 'G');
                break;
        }

        return isComplement;
    }

    /**
     *
     * @param source
     * @return
     */
    public static String reverseIt(String source) {
        int i, len = source.length();
        StringBuilder dest = new StringBuilder(len);

        for (i = (len - 1); i >= 0; i--) {
            dest.append(source.charAt(i));
        }

        return dest.toString();
    }

    public static boolean isRNAPair(char base1, char base2) {

        if (!isComplementaryRNA(base1, base2)) {
            return isComplementaryRNAWooble(base1, base2);
        } else {
            return true;
        }
    }

    /**
     *
     * @param base1
     * @param base2
     * @return
     */
    public static boolean isComplementaryRNA(char base1, char base2) {
        // ok
        boolean isComplement = false;

        /**
         * The nucleic acid codes are: A --> adenosine M --> A C (amino) C -->
         * cytidine S --> G C (strong) G --> guanine W --> A T (weak) T -->
         * thymidine B --> G T C U --> uridine D --> G A T R --> G A (purine) H
         * --> A C T Y --> T C (pyrimidine) V --> G C A K --> G T (keto) N --> A
         * G C T (any) - gap of indeterminate length
         */
        switch (base1) {
            case 'A':
                isComplement = (base2 == 'U') || (base2 == 'W')
                        || (base2 == 'K') || (base2 == 'Y')
                        || (base2 == 'B') || (base2 == 'D')
                        || (base2 == 'H') || (base2 == 'V');
                break;
            case 'U':
                isComplement = (base2 == 'A') || (base2 == 'W')
                        || (base2 == 'M') || (base2 == 'R')
                        || (base2 == 'D')
                        || (base2 == 'H') || (base2 == 'V');
                break;
            case 'G':
                isComplement = (base2 == 'C') || (base2 == 'S')
                        || (base2 == 'M') || (base2 == 'Y')
                        || (base2 == 'B')
                        || (base2 == 'H') || (base2 == 'V');
                break;
            case 'C':
                isComplement = (base2 == 'G') || (base2 == 'S')
                        || (base2 == 'K') || (base2 == 'R')
                        || (base2 == 'B');
                break;
        }

        return isComplement;
    }

    /**
     *
     * @param base1
     * @param base2
     * @return
     */
    public static boolean isMismatch(char base1, char base2) {
        boolean isMismatch;
        isMismatch = !isComplementaryDNA(base1, base2)
                && !isComplementaryRNA(base1, base2)
                && !isComplementaryRNAWooble(base1, base2);

        return isMismatch;
    }

    /**
     * The nucleic acid codes are: A --> adenosine M --> A C (amino) C -->
     * cytidine S --> G C (strong) G --> guanine W --> A T (weak) T -->
     * thymidine B --> G T C U --> uridine D --> G A T R --> G A (purine) H -->
     * A C T Y --> T C (pyrimidine) V --> G C A K --> G T (keto) N --> A G C T
     * (any) - gap of indeterminate length
     *
     * @param base1
     * @param base2
     * @return
     */
    public static boolean isComplementaryRNAWooble(char base1, char base2) {

        boolean isComplement = false;

        switch (base1) {
            case 'U':
                isComplement = (base2 == 'G') || (base2 == 'S')
                        || (base2 == 'K') || (base2 == 'R')
                        || (base2 == 'B') || (base2 == 'I');
                break;

            case 'G':
                isComplement = (base2 == 'U') || (base2 == 'W')
                        || (base2 == 'K') || (base2 == 'Y')
                        || (base2 == 'B') || (base2 == 'D')
                        || (base2 == 'H') || (base2 == 'V');
                break;

            case 'I':
                isComplement = (base2 == 'U') || (base2 == 'A')
                        || (base2 == 'C');
                break;

            case 'A':
                isComplement = (base2 == 'I');
                break;

            case 'C':
                isComplement = (base2 == 'I');
                break;

        }

        return isComplement;
    }

    /**
     *
     * @param dnaSequence
     * @return
     */
    public static String toComplementaryRNA(String dnaSequence) {
        String rnaSequence = ""; //$NON-NLS-1$
        int i;

        for (i = 0; i < dnaSequence.length(); i++) {
            switch (dnaSequence.charAt(i)) {
                case 'A':
                    rnaSequence = rnaSequence.concat("U"); //$NON-NLS-1$
                    break;
                case 'T':
                    rnaSequence = rnaSequence.concat("A"); //$NON-NLS-1$
                    break;
                case 'G':
                    rnaSequence = rnaSequence.concat("C"); //$NON-NLS-1$
                    break;
                case 'C':
                    rnaSequence = rnaSequence.concat("G"); //$NON-NLS-1$
                    break;
                default:
                    rnaSequence = rnaSequence.concat(String.valueOf(
                            dnaSequence.charAt(i))); // test
            }
        }

        return rnaSequence;
    }

    /**
     *
     * @param fastaPattern
     * @return
     */
    public static String toRegularExpression(String fastaPattern) {
        String regExp = null;
        boolean extendedReplace = true;

        regExp = fastaPattern.replaceAll("N", "[AUTGC]"); //$NON-NLS-1$ 

        if (extendedReplace) {
            regExp = regExp.replaceAll("W", "[AUT]"); //$NON-NLS-1$ 
            regExp = regExp.replaceAll("S", "[CG]"); //$NON-NLS-1$ 
            regExp = regExp.replaceAll("M", "[AC]"); //$NON-NLS-1$
            regExp = regExp.replaceAll("K", "[GUT]"); //$NON-NLS-1$ 
            regExp = regExp.replaceAll("R", "[GA]"); //$NON-NLS-1$ 
            regExp = regExp.replaceAll("Y", "[CUT]"); //$NON-NLS-1$ 
            regExp = regExp.replaceAll("B", "[GCUT]"); //$NON-NLS-1$ 
            regExp = regExp.replaceAll("D", "[GUTA]"); //$NON-NLS-1$
            regExp = regExp.replaceAll("H", "[CUTA]"); //$NON-NLS-1$ 
            regExp = regExp.replaceAll("V", "[CGA]"); //$NON-NLS-1$ 
        }
        return regExp;
    }

    /**
     *
     * @param rna
     * @param sequence
     * @param nMin
     * @return
     */
    public static int findSlippageSequence(String rna, String sequence,
            int nMin) {
        /**
         * To search for repeated tracts as CA[N]
         */
        int init = 0, n = 0, position = 0, indexOf;

        indexOf = rna.indexOf(sequence, init);

        while (indexOf > 0) {

            if (init == 0) {
                position = indexOf;
            }

            init = indexOf + sequence.length();
            n++;

            indexOf = rna.indexOf(sequence, init);
        }

        if (n >= nMin) {
            return position;
        } else {
            return 0;
        }
    }

    /**
     *
     * @param loop
     * @return
     */
    public static boolean checkInternalPairments(String loop) {
        char base1, base2;
        boolean ret = true;

        // Check if loop is good (There's no pairment between his bases)
        if (loop.length() > 5) {
            for (int i = 0; i < 1/*loop.length()/2 - 1*/; i++) {

                base1 = loop.charAt(i);
                base2 = loop.charAt(loop.length() - 1 - i);

                if (SequenceAnalizer.isComplementaryRNAWooble(base1, base2)
                        || SequenceAnalizer.isComplementaryRNA(base1, base2)) {
                    ret = false;
                }
            }
        }
        return ret;
    }

}
