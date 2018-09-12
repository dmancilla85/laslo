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
import java.util.Arrays;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.Semaphore;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.apache.commons.lang.StringUtils;
import org.biojava.nbio.core.sequence.DNASequence;

public class LoopCatcherThread implements Runnable {

    protected boolean extendedMode;
    protected String additionalSequence;
    protected int maxLength;
    protected int minLength;
    protected DNASequence dnaElement;
    protected InputSequence inputType;
    protected Iterator<String> patternItr;
    protected CSVWriter writer;
    private final static Semaphore MUTEX = new Semaphore(1);
    //private final static Semaphore SEM = new Semaphore(4);
    private static Semaphore SEM;
    private CountDownLatch latch;

    public void setLatch(CountDownLatch latch) {
        this.latch = latch;
    }

    public static void setCores(){
        SEM = new Semaphore(OSValidator.getNumberOfCPUCores());
    }
    
    public LoopCatcherThread(boolean extendedMode, String additionalSequence,
            int maxLength, int minLength, DNASequence dnaElement,
            InputSequence inputType, Iterator<String> patternItr,
            CSVWriter writer) {
        this.extendedMode = extendedMode;
        this.additionalSequence = additionalSequence;
        this.maxLength = maxLength;
        this.minLength = minLength;
        this.dnaElement = dnaElement;
        this.inputType = inputType;
        this.patternItr = patternItr;
        this.writer = writer;
    }

    @Override
    public void run() {

        String rnaSequence = "";
        RNAfold fold = new RNAfold();

        // sólo para modo dos
        if (extendedMode) {
            rnaSequence = dnaElement.getRNASequence().getSequenceAsString();
            fold = new RNAfold(rnaSequence);
        }

        // III. Loop level
        while (patternItr.hasNext()) {

            String currentPattern = patternItr.next().trim().toUpperCase();

            try {
                SEM.acquire();
            } catch (InterruptedException ex) {
                Logger.getLogger(LoopCatcherThread.class.getName()).log(Level.SEVERE, null, ex);
            } finally {
                SEM.release();
            }
            
            // 1. Stem research
            if (extendedMode) {
                sequenceExtendedResearch(
                        rnaSequence,
                        dnaElement.getOriginalHeader(),
                        fold.getStructure(),
                        currentPattern,
                        writer);
            } else {
                sequenceResearch(dnaElement, currentPattern, writer);
            }
        }
       
        out.println("Threads remaining: " + latch.getCount());
        latch.countDown();
        
    }

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
     * @return
     */
    public int sequenceResearch(DNASequence fastaSeq,
            String stemLoopPattern, CSVWriter writer) {

        List<StemLoop> slrList = new ArrayList<>();
        //List<Integer> mismatchs = new ArrayList<>();
        StemLoop slr;
        slr = null;
        int size, posAux, k = 3;
        size = 0;
        int loopPos = 0, loopLength = stemLoopPattern.length();
        boolean isValidHairpin;
        String rnaSequence = fastaSeq.getRNASequence().getSequenceAsString();
        final int sequenceLength = rnaSequence.length();

        RNAfold fold;

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
            slr.setTags(fastaSeq.getOriginalHeader());
            try {
                if ((loopPos - length) > 0
                        && (loopPos + loopLength + length) < sequenceLength) {
                    rnaSeq = rnaSequence.substring(loopPos - length,
                            loopPos + loopLength + length);

                    isValidHairpin = isComplementaryRNA(rnaSeq.charAt(length-1),
                            rnaSeq.charAt(length + loopLength))
                            || isComplementaryRNAWooble(rnaSeq.charAt(length-1),
                                    rnaSeq.charAt(length + loopLength));

                    if (isValidHairpin) {
                        fold = new RNAfold(rnaSeq);
                        hairpinModel = fold.getStructure();

                        if (rnaSeq.length() != hairpinModel.length()) {
                            out.println("Error NO COINCIDEN:" + rnaSeq + " - "
                                    + hairpinModel);
                        }

                        if (fold.getMfe()
                                //mfe.mfe 
                                == 0.0) {
                            isValidHairpin = false;
                        } else {
                            hairpinModel = isValidHairpin(hairpinModel, loopLength, loopPos, rnaSeq);
                            isValidHairpin = hairpinModel.length() > 0;
                        }
                    }

                } else {
                    isValidHairpin = false;
                }

            } catch (Exception e) {
                out.println(Arrays.toString(e.getStackTrace()));
            }

            if (isValidHairpin) {
                int extIzq = hairpinModel.lastIndexOf("(") + 1; //NOI18N
                int extDer = hairpinModel.length() - extIzq - loopLength;
                rnaSeq = rnaSequence.substring(loopPos - extIzq, loopPos
                        + loopLength + extDer);

                posAux = rnaSequence.indexOf(rnaSeq);

                // Fill the fields
                slr.setRnaHairpinSequence(rnaSeq);
                slr.setLoop(rnaLoop);
                slr.setStartsAt(loopPos - extIzq);
                slr.setStructure(hairpinModel);
                slr.setSequenceLength(sequenceLength);
                try {
                    slr.setAdditional5Seq(rnaSequence
                            .substring(posAux - k, posAux));
                    slr.setAdditional3Seq(rnaSequence.substring(posAux
                            + rnaSeq.length(), posAux + rnaSeq.length() + k));
                } catch (IndexOutOfBoundsException e) {
                    slr.setAdditional3Seq("");
                    slr.setAdditional5Seq("");
                }
                slr.checkPairments();
                slr.setMfe(new RNAfold(rnaSeq).getMfe()
                //rfa.getMFE(rnaSeq.getBytes()).mfe
                );
                slr.setPredecessorLoop(extIzq);
                slr.setNLoop(extIzq);
                slr.setPercent_AG();
                slr.setLoopPattern(stemLoopPattern);
                slr.setEndsAt(loopFinder.end() + extDer);
                slr.setPercA_sequence(
                        (rnaSequence.length() - rnaSequence.replace("A", "") //NOI18N
                        .length()) / (float) rnaSequence.length());
                slr.setPercG_sequence(
                        (rnaSequence.length() - rnaSequence.replace("G", "") //NOI18N
                        .length()) / (float) rnaSequence.length());
                slr.setPercC_sequence(
                        (rnaSequence.length() - rnaSequence.replace("C", "") //NOI18N
                        .length()) / (float) rnaSequence.length());
                slr.setPercU_sequence(
                        (rnaSequence.length() - rnaSequence.replace("U", "") //NOI18N
                        .length()) / (float) rnaSequence.length());

                if(this.additionalSequence.length() > 0)
                slr.setAdditionalSeqLocations(getPatternLocations(rnaSequence,
                        this.additionalSequence));

                slr.setRelativePos((double) slr.getStartsAt() / (double) rnaSequence.length());

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
                Logger.getLogger(LoopCatcherThread.class.getName()).log(Level.SEVERE, null, ex);
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
    public String isValidHairpin(String hairpin, int loopLength, int loopPos, String seq) {
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
            // 2.a 	Si encuentra un . o un ) en el extremo izquierdo, remover 
            // 		base (stemLength-1)
            // Repetir 2.a hasta que complete el recorrido. 
            extremoIzq = hairpin.substring(0, stemLength);

            try {
                for (int i = 0; i < stemLength; i++) {
                    if (extremoIzq.charAt(i) == '.' || extremoIzq.charAt(i) == ')') {
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

                // 2.b 	Si encuentra un . o un ( en el extremo derecho, remover base 
                //		(stemLength-1)
                // Repetir 2.b hasta que complete el recorrido. 	
                extremoDer = hairpin.substring(stemLength + loopLength,
                        hairpin.length());
                begin = true;

                for (int i = stemLength - 1; i >= 0; i--) {
                    if (extremoDer.charAt(i) == '.' || extremoDer.charAt(i) == '(') {
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
                    && //NOI18N
                    StringUtils.countMatches(hairpin, ")") >= minLength; //NOI18N
        }
        // Si la longitud es menor a la esperada, rechazar.
        if (!ret) {
            return ""; //NOI18N
        } else {
        }
        return hairpin;
    }

    public int sequenceExtendedResearch(String rnaSequence, String header,
            String hairpinSeq, String stemLoopPattern, CSVWriter writer) {
        List<StemLoop> slrList = new ArrayList<>();
        StemLoop slr;
        slr = null;
        int size;
        size = 0;
        int loopPos = 0, loopLength = stemLoopPattern.length();
        boolean isValidHairpin;

        final int sequenceLength = rnaSequence.length();

        // Convert the original loop pattern to a regular expression
        String regExp = toRegularExpression(stemLoopPattern);
        Pattern p = Pattern.compile(regExp);
        Matcher loopFinder = p.matcher(rnaSequence);

        String rnaLoop = "", rnaSeq = "", hairpinModel = ""; //NOI18N

        // As exists loop matches
        while (loopFinder.find()) {
            loopPos = loopFinder.start();
            //mismatchs.clear();
            isValidHairpin = true;
            rnaLoop = rnaSequence.substring(loopPos, loopFinder.end());
            int length = this.maxLength, posAux, k = 3;
            slr = new StemLoop(this.inputType);
            slr.setTags(header);
            try {
                if ((loopPos - length) > 0
                        && (loopPos + loopLength + length) < sequenceLength) {

                    rnaSeq = rnaSequence.substring(loopPos - length,
                            loopPos + loopLength + length);
                    hairpinModel = hairpinSeq.substring(loopPos - length,
                            loopPos + loopLength + length);
                    
                    isValidHairpin = isComplementaryRNA(rnaSeq.charAt(length-1),
                            rnaSeq.charAt(length + loopLength))
                            || isComplementaryRNAWooble(rnaSeq.charAt(length-1),
                                    rnaSeq.charAt(length + loopLength));

                    if (rnaSeq.length() != hairpinModel.length()) {
                        out.println("Error NO COINCIDEN:" + rnaSeq + " - "
                                + hairpinModel);
                    }

                    if(isValidHairpin){
                        hairpinModel = isValidHairpin(hairpinModel, loopLength, 
                                loopPos, rnaSeq);
                        isValidHairpin = hairpinModel.length() > 0;
                    }

                } else {
                    isValidHairpin = false;
                }

            } catch (Exception e) {
                out.println(Arrays.toString(e.getStackTrace()));
            }

            if (isValidHairpin) {
                int extIzq = hairpinModel.lastIndexOf("(") + 1; //NOI18N
                int extDer = hairpinModel.length() - extIzq - loopLength;
                rnaSeq = rnaSequence.substring(loopPos - extIzq, loopPos
                        + loopLength + extDer);

                posAux = rnaSequence.indexOf(rnaSeq);

                // Fill the fields
                slr.setRnaHairpinSequence(rnaSeq);
                slr.setLoop(rnaLoop);
                slr.setStartsAt(loopPos - extIzq);
                slr.setStructure(hairpinModel);
                slr.setSequenceLength(sequenceLength);
                try {
                    slr.setAdditional5Seq(rnaSequence
                            .substring(posAux - k, posAux));
                    slr.setAdditional3Seq(rnaSequence.substring(posAux
                            + rnaSeq.length(), posAux + rnaSeq.length() + k));
                } catch (IndexOutOfBoundsException e) {
                    slr.setAdditional3Seq("");
                    slr.setAdditional5Seq("");
                }
                slr.checkPairments();
                slr.setMfe(new RNAfold(rnaSeq).getMfe());
                slr.setPredecessorLoop(extIzq);
                slr.setNLoop(extIzq);
                slr.setPercent_AG();
                slr.setLoopPattern(stemLoopPattern);
                slr.setEndsAt(loopFinder.end() + extDer);
                slr.setPercA_sequence(
                        (rnaSequence.length() - rnaSequence.replace("A", "") //NOI18N
                        .length()) / (float) rnaSequence.length());
                slr.setPercG_sequence(
                        (rnaSequence.length() - rnaSequence.replace("G", "") //NOI18N
                        .length()) / (float) rnaSequence.length());
                slr.setPercC_sequence(
                        (rnaSequence.length() - rnaSequence.replace("C", "") //NOI18N
                        .length()) / (float) rnaSequence.length());
                slr.setPercU_sequence(
                        (rnaSequence.length() - rnaSequence.replace("U", "") //NOI18N
                        .length()) / (float) rnaSequence.length());

                // Other signals in cDNA
                slr.setAdditionalSeqLocations(getPatternLocations(rnaSequence,
                        this.additionalSequence));
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
                Logger.getLogger(LoopCatcherThread.class.getName()).log(Level.SEVERE, null, ex);
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
