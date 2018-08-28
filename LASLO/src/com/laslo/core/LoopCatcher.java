/**
 *
 */
package com.laslo.core;

import com.tools.ShuffleSeq;
import static com.laslo.core.PairmentAnalizer.*;
import com.tools.fasta.InputSequence;
import static java.lang.System.out;
import java.io.FileWriter;
import java.io.IOException;
import java.io.File;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import com.opencsv.CSVWriter;
import com.tools.RNAfold;
import com.tools.fasta.FASTACorrector;
import static com.tools.fasta.FastaID.CSV_EXT;
import static com.tools.fasta.FastaID.FASTA_EXT;
import static com.tools.fasta.FastaID.FASTA_EXT_2;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Locale;
import java.util.Map;
import static java.util.ResourceBundle.getBundle;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.lang.StringUtils;
import org.biojava.nbio.core.sequence.DNASequence;
import static org.biojava.nbio.core.sequence.io.FastaReaderHelper.readFastaDNASequence;
import rnafold4j.MFEData;
import rnafold4j.RNAFoldAPI;

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
    protected Locale currentLocale;
    protected int maxMismatch;
    protected File[] fileList;
    protected boolean extendedMode;
    private boolean makeRandoms;
    private int numberOfRandoms;
    protected File actualFile;

    // Create a new instance of RNAFold4J
    protected RNAFoldAPI rfa;

    public LoopCatcher(String pathOut, String pathIn,
            ArrayList<String> loopPatterns, InputSequence inputType,
            int minLength, int maxLength,
            int maxWooble, int maxMismatch, Locale locale) {
        this.pathOut = pathOut;
        this.loopPatterns = loopPatterns;
        this.inputType = inputType;
        this.minLength = minLength;
        this.maxLength = maxLength;
        this.maxWooble = maxWooble;
        this.maxMismatch = maxMismatch;
        this.rfa = new RNAFoldAPI();
        this.fileList = null;
        this.extendedMode = false;
        this.makeRandoms = false;
        this.numberOfRandoms = 0;
        this.currentLocale = locale;
    }

    /**
     *
     */
    public LoopCatcher() {
        this("", "", new ArrayList<>(), InputSequence.ENSEMBL, //NOI18N
                4, 16, 2, 0, new Locale("es", "AR"));
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

    public void setCurrentLocale(Locale locale){
        this.currentLocale  = locale;
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
    public void writeCSV(String sFileName, List<StemLoop> stemResearch, String headerList) {
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
     * @param fastaSeq
     * @param stemLoopPattern
     * @param stemsDetected
     * @param writer
     * @return
     */
    public int sequenceResearch(DNASequence fastaSeq, String stemLoopPattern,
            ArrayList<StemPosition> stemsDetected, CSVWriter writer) {

        List<StemLoop> slrList = new ArrayList<>();
        StemLoop slr;
        slr = null;
        List<Integer> mismatchs = new ArrayList<>();
        int size;
        String hairpinModel;
        size = 0;

        int loopPos = 0, // Position of the stem loop
                mismatchCount, woobleCount,
                // Number of G=U pairs. Typically 1, max 4
                lengthIR, // Length of the inverted repeat sequences
                loopLength = stemLoopPattern.length(), pos1,
                pos2;
        boolean existsInvertedRepeat;
        existsInvertedRepeat = false;
        boolean checkWooblePair, isValidHairpin;
        // Upper all characters
        // Make the retro transcription
        String rnaSequence = fastaSeq.getSequenceAsString().toUpperCase().replace('T', 'U');
        final int sequenceLength = rnaSequence.length();

        // Convert the original loop pattern to a regular expression
        String regExp = toRegularExpression(stemLoopPattern);
        Pattern p = Pattern.compile(regExp);
        Matcher loopFinder = p.matcher(rnaSequence);

        String rnaLoop, rnaSeq;

        // As exists loop matches
        while (loopFinder.find()) {
            existsInvertedRepeat = true;
            loopPos = loopFinder.start();
            lengthIR = 0;
            mismatchCount = 0;
            woobleCount = 0;
            mismatchs.clear();
            isValidHairpin = true;

            rnaLoop = rnaSequence.substring(loopFinder.start(), loopFinder.end());

            if (!checkInternalPairments(rnaLoop)) {
                existsInvertedRepeat = false;
            }

            while (existsInvertedRepeat) {
                // To check that the search doesn't exceed the bounds of the
                // sequence
                checkWooblePair = (woobleCount < this.maxWooble);
                pos1 = loopPos - (lengthIR + 1);
                pos2 = loopPos + loopLength + lengthIR;

                if (pos1 >= 0 && pos2 < sequenceLength) {
                    try {

                        if (isComplementaryRNA(rnaSequence.
                                charAt(pos1), rnaSequence.charAt(pos2))
                                || (isComplementaryRNAWooble(rnaSequence
                                        .charAt(pos1),
                                        rnaSequence.charAt(pos2))
                                && checkWooblePair)
                                || mismatchCount < this.maxMismatch) {
                            existsInvertedRepeat = true;

                            if (!isComplementaryRNA(
                                    rnaSequence.charAt(pos1), rnaSequence
                                    .charAt(pos2))) {
                                // A mismatch shouldn't be considered if is at
                                // the border

                                if (isComplementaryRNAWooble(rnaSequence
                                        .charAt(pos1),
                                        rnaSequence.charAt(pos2))
                                        && woobleCount < this.maxWooble) {
                                    woobleCount++;

                                } else if (lengthIR > 0 && mismatchCount
                                        < this.maxMismatch) {
                                    mismatchCount++;
                                    mismatchs.add(lengthIR + 1);
                                    mismatchs.add(loopPos + loopLength
                                            + lengthIR);
                                } else {
                                    existsInvertedRepeat = false;
                                }
                            }

                            lengthIR++;
                        } else {
                            existsInvertedRepeat = false;
                        }

                    } catch (IndexOutOfBoundsException e) {
                        out.println("Error: " + rnaSequence + " - i: "
                                + (loopPos - lengthIR) + " - i2: " //$NON-NLS-3$ //$NON-NLS-2$ //NOI18N
                                // //$NON-NLS-3$
                                + (loopPos + loopLength + lengthIR)); // $NON-NLS-1$
                        out.flush();

                    }
                } else {
                    existsInvertedRepeat = false;
                }
            }

            rnaSeq = rnaSequence.substring(loopPos - lengthIR, loopPos
                    + loopLength + lengthIR);

            // Remove mismatchs from borders
            boolean cut = false;

            if (mismatchCount > 0) {
                for (int i = 0; i < (rnaSeq.length() - loopLength) / 2 && !cut
                        && mismatchCount > 0; i++) {
                    if (isMismatch(rnaSeq.charAt(i),
                            rnaSeq.charAt(rnaSeq.length() - 1 - i))) {
                        lengthIR--;
                        mismatchCount--;

                    } else {
                        cut = true;
                    }
                }
            }

            // Control the percent of mismatches/totalLength
            if (this.minLength <= lengthIR && lengthIR <= this.maxLength) {
                if (mismatchCount > 0 && lengthIR < 5) {
                    isValidHairpin = false;
                }

            } else {
                isValidHairpin = false;
            }

            slr = new StemLoop(this.inputType);
            slr.setTags(fastaSeq.getOriginalHeader());
            slr.setStartsAt(loopPos - lengthIR + 1);

            /*
            * if( stemsDetected.contains(new
            * StemPosition(slr.getTranscriptID(), slr.getStartsAt())) ){
            * isValidHairpin = false; slr = null; } else
            * stemsDetected.add(new StemPosition(slr.getTranscriptID(),
            * slr.getStartsAt()));
             */
            rnaSeq = rnaSequence.substring(loopPos - lengthIR, loopPos
                    + loopLength + lengthIR);

            MFEData mfe = rfa.getMFE(rnaSeq.getBytes());

            if (isValidHairpin) {
                slr.setRnaHairpinSequence(rnaSeq);

                slr.setLoop(rnaSequence.substring(loopPos, loopPos
                        + loopLength));
                hairpinModel = slr.drawHairpinViennaStructure();

                if (mfe.mfe == 0 || !hairpinModel
                        .equals(new String(mfe.structure).trim())) {
                    isValidHairpin = false;
                }
            }

            if (isValidHairpin) {
                // Fill the fields
                slr.setSequenceLength(sequenceLength);
                slr.setMfe(mfe.mfe);
                slr.setStructure(new String(mfe.structure));
                slr.setPredecessorLoop(lengthIR);
                slr.setNLoop(lengthIR);
                slr.setPercent_AG();
                slr.setLoopPattern(stemLoopPattern);
                slr.setMismatches(mismatchCount);
                slr.setPercent_GU(woobleCount);
                slr.setPercent_AU();
                slr.setPercent_CG();
                slr.setStartsAt(loopPos - lengthIR + 1);
                slr.setEndsAt(loopFinder.end() + lengthIR + 1);
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
                // ----------------------
                // PUM1 Protein (pumilio)
                slr.setPumilioLocations(getPatternLocations(rnaSequence,
                        BiologyPatterns.PUMILIO));
                // Polyadenilations Sites
                // Look for the Kozak pattern
                // Special TEST buscar otros posibles stems
                slr.setRelativePos((double) slr.getStartsAt()
                        / (double) rnaSequence.length());

                slrList.add(slr);
            }

        }

        slr = null;
        mismatchs = null;

        Iterator<StemLoop> itr = slrList.iterator();

        while (itr.hasNext()) {
            StemLoop element = itr.next();
            writer.writeNext(element.toRowCSV().split(";")); //NOI18N
        }

        size = slrList.size();
        slrList.clear();
        slrList = null;

        return size;
    }

    /**
     *
     * @param fastaSeq
     * @param stemLoopPattern
     * @param stemsDetected
     * @param writer
     * @return
     */
    public int sequenceExtendedResearch(DNASequence fastaSeq,
            String stemLoopPattern, ArrayList<StemPosition> stemsDetected,
            CSVWriter writer) {

        List<StemLoop> slrList = new ArrayList<>();
        List<Integer> mismatchs = new ArrayList<>();
        StemLoop slr;
        slr = null;
        int size;
        size = 0;
        int loopPos = 0,
                loopLength = stemLoopPattern.length();
        boolean isValidHairpin;
        String rnaSequence = fastaSeq.getSequenceAsString().
                toUpperCase().replace('T', 'U');
        final int sequenceLength = rnaSequence.length();
        MFEData mfe = null;
        //RNAfold fold = null;

        // Convert the original loop pattern to a regular expression
        String regExp = toRegularExpression(stemLoopPattern);
        Pattern p = Pattern.compile(regExp);
        Matcher loopFinder = p.matcher(rnaSequence);

        String rnaLoop = "", rnaSeq = "", hairpinModel = ""; //NOI18N

        // As exists loop matches
        while (loopFinder.find()) {
            loopPos = loopFinder.start();
            mismatchs.clear();
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

                    mfe = rfa.getMFE(rnaSeq.getBytes());
                    //fold = new RNAfold(rnaSeq);

                    hairpinModel = new String(mfe.structure);
                    //fold.getStructure();

                    if (/*fold.getMfe()*/mfe.mfe == 0.0) {
                        isValidHairpin = false;
                    } else {
                        hairpinModel = isValidHairpin(hairpinModel, loopLength);
                        isValidHairpin = hairpinModel.length() > 0;
                    }

                } else {
                    isValidHairpin = false;
                }

            } catch (Exception e) {
                /*out.println("ERROR!: Length" + length + " - Hairpin: "
                        + hairpinModel + " - seq: " + rnaSeq + ". Exception: "
                        + e.getMessage());*/
                out.println(Arrays.toString(e.getStackTrace()));
            }

            if (isValidHairpin) {
                int extIzq = hairpinModel.lastIndexOf("(") + 1; //NOI18N
                int extDer = hairpinModel.length() - extIzq - loopLength;
                rnaSeq = rnaSequence.substring(loopPos - extIzq, loopPos
                        + loopLength + extDer);

                // Fill the fields
                slr.setRnaHairpinSequence(rnaSeq);
                slr.setLoop(rnaLoop);
                slr.setStartsAt(loopPos - extIzq);
                slr.setStructure(hairpinModel);
                slr.setSequenceLength(sequenceLength);
                slr.checkPairments();
                slr.setMfe(/*fold.getMfe()*/mfe.mfe);
                slr.setPredecessorLoop(extIzq);
                slr.setNLoop(extIzq);
                slr.setPercent_AG();
                slr.setLoopPattern(stemLoopPattern);
                //slr.setPercent_AU();
                //slr.setPercent_CG();
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
                // ----------------------
                // PUM1 Protein (pumilio)
                slr.setPumilioLocations(getPatternLocations(rnaSequence,
                        BiologyPatterns.PUMILIO));
                // Polyadenilations Sites
                // Look for the Kozak pattern
                // Special TEST buscar otros posibles stems
                slr.setRelativePos((double) slr.getStartsAt() / (double) rnaSequence.length());

                slrList.add(slr);
            }

        }

        slr = null;
        mismatchs = null;

        Iterator<StemLoop> itr = slrList.iterator();

        while (itr.hasNext()) {
            StemLoop element = itr.next();
            writer.writeNext(element.toRowCSV().split(";")); //NOI18N
        }

        size = slrList.size();
        slrList.clear();
        slrList = null;

        return size;
    }
    
    @SuppressWarnings("empty-statement")
    public String isValidHairpin(String hairpin,
            int loopLength) {
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
        ret = ret && hairpin.charAt(stemLength - 1) == '('
                && hairpin.charAt(stemLength + loopLength) == ')';

        // 2. Analizar lado izq. y derecho del stem desde extremo 5',
        if (ret) {
            // 2.a Si encuentra un . o un ) en el extremo izquierdo, remover base (stemLength-1)
            // repetir 2.a hasta que complete el recorrido. 
           extremoIzq = hairpin.substring(0, stemLength);
           
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

            // 2.b Si encuentra un . o un ( en el extremo derecho, remover base (stemLength-1)
            // repetir 2.b hasta que complete el recorrido. 	
            extremoDer = hairpin.substring(stemLength + loopLength, hairpin.length());
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
            while(auxL > 0 && auxR > 0){
                if(hairpin.charAt(i--) == '('){
                    auxR--;
                    auxL--;
                }
            }
            
            hairpin = hairpin.substring(i + 1, hairpin.length());
            
            int k;
            for(k=hairpin.length()-1; hairpin.charAt(k) == '.';k--);
            
            hairpin = hairpin.substring(0, k + 1);
                
        }
        
        if(ret)
            ret = StringUtils.countMatches(hairpin, "(") >= minLength && //NOI18N
                StringUtils.countMatches(hairpin, ")") >= minLength; //NOI18N
                
        // Si la longitud es menor a la esperada, rechazar.
        if (!ret) {
            return ""; //NOI18N
        } else {
        }
        return hairpin;
    }

    public int searchMode2(String rnaSequence, String header,
            RNAfold fold,
            String stemLoopPattern, ArrayList<StemPosition> stemsDetected,
            CSVWriter writer) {
        List<StemLoop> slrList = new ArrayList<>();
        List<Integer> mismatchs = new ArrayList<>();
        StemLoop slr;
        slr = null;
        int size;
        size = 0;
        int loopPos = 0,
                loopLength = stemLoopPattern.length();
        boolean isValidHairpin;
   
        final int sequenceLength = rnaSequence.length();
        String  fullStr = ""; //NOI18N

        // Convert the original loop pattern to a regular expression
        String regExp = toRegularExpression(stemLoopPattern);
        Pattern p = Pattern.compile(regExp);
        Matcher loopFinder = p.matcher(rnaSequence);

        String rnaLoop = "", rnaSeq = "", hairpinModel = ""; //NOI18N
        fullStr = fold.getStructure();
        
        // As exists loop matches
        while (loopFinder.find()) {
            loopPos = loopFinder.start();
            mismatchs.clear();
            isValidHairpin = true;
            rnaLoop = rnaSequence.substring(loopPos, loopFinder.end());
            int length = this.maxLength;
            slr = new StemLoop(this.inputType);
            slr.setTags(header);
            try {
                if ((loopPos - length) > 0
                        && (loopPos + loopLength + length) < sequenceLength) {
                    
                    rnaSeq = rnaSequence.substring(loopPos - length, loopPos
                            + loopLength + length);

                    hairpinModel = fullStr.substring(loopPos - length, loopPos
                            + loopLength + length);;
                    //fold.getStructure();

                    if (fold.getMfe()== 0.0) {
                        isValidHairpin = false;
                    } else {
                        hairpinModel = isValidHairpin(hairpinModel, loopLength);
                        isValidHairpin = hairpinModel.length() > 0;
                    }

                } else {
                    isValidHairpin = false;
                }

            } catch (Exception e) {
                /*out.println("ERROR!: Length" + length + " - Hairpin: "
                        + hairpinModel + " - seq: " + rnaSeq + ". Exception: "
                        + e.getMessage());*/
                out.println(Arrays.toString(e.getStackTrace()));
            }

            if (isValidHairpin) {
                length = (hairpinModel.length() - loopLength) / 2;
                rnaSeq = rnaSequence.substring(loopPos - length, loopPos
                        + loopLength + length);

                // Fill the fields
                slr.setRnaHairpinSequence(rnaSeq);
                slr.setLoop(rnaLoop);
                slr.setStartsAt(loopPos - length);
                slr.setStructure(hairpinModel);
                slr.setSequenceLength(sequenceLength);
                slr.checkPairments();
                slr.setMfe(fold.getMfe());
                slr.setPredecessorLoop(length);
                slr.setNLoop(length);
                slr.setPercent_AG();
                slr.setLoopPattern(stemLoopPattern);
                slr.setPercent_AU();
                slr.setPercent_CG();
                slr.setStartsAt(loopPos - length + 1);
                slr.setEndsAt(loopFinder.end() + length + 1);
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
                // ----------------------
                // PUM1 Protein (pumilio)
                slr.setPumilioLocations(getPatternLocations(rnaSequence,
                        BiologyPatterns.PUMILIO));
                // Polyadenilations Sites
                // Look for the Kozak pattern
                // Special TEST buscar otros posibles stems
                slr.setRelativePos((double) slr.getStartsAt() / (double) rnaSequence.length());

                slrList.add(slr);
            }

        }

        slr = null;
        mismatchs = null;

        Iterator<StemLoop> itr = slrList.iterator();

        while (itr.hasNext()) {
            StemLoop element = itr.next();
            writer.writeNext(element.toRowCSV().split(";")); //NOI18N
        }

        size = slrList.size();
        slrList.clear();
        slrList = null;

        return size;
    }

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
        String fileName;

        ini = Calendar.getInstance();
        out.println(
                java.text.MessageFormat.format(getBundle("resources/Bundle", currentLocale).
                        getString("START_TIME"), new Object[] {Calendar.getInstance().getTime()}));
        out.flush();

        if (fileList == null) {
            return false;
        }

        if (this.makeRandoms) {
            for (File currentFile : fileList) {
                if (currentFile.isFile() && (currentFile.toString().endsWith(FASTA_EXT)
                        || currentFile.toString().endsWith(FASTA_EXT_2))) {

                    fileName = currentFile.getName();

                    if (fileName.contains(FASTA_EXT)) {
                        fileName = fileName.replaceAll(FASTA_EXT, ""); //NOI18N
                        ShuffleSeq.makeShuffleSequences(pathOut, fileName, FASTA_EXT, numberOfRandoms);
                    } else {
                        fileName = fileName.replaceAll(FASTA_EXT_2, ""); //NOI18N
                        ShuffleSeq.makeShuffleSequences(pathOut, fileName, FASTA_EXT, numberOfRandoms);
                    }

                }
            }

            File folder;
            folder = new File(pathIn + ShuffleSeq.getRandomDir());
            File[] randomFiles;
            randomFiles = folder.listFiles();
            fileList = unionFiles(fileList, randomFiles);
        }

        // I. File level (hacer un hilo)
        for (File currentFile : fileList) {
            if (currentFile.isFile() && (currentFile.toString().endsWith(FASTA_EXT)
                    || currentFile.toString().endsWith(FASTA_EXT_2))) {

                this.actualFile = currentFile;
                run();
            }
        }

        fin = Calendar.getInstance();
        out.println(java.text.MessageFormat.format(getBundle("resources/Bundle", currentLocale).getString("TOTAL_TIME"), 
                new Object[] {((fin.getTimeInMillis() - ini.getTimeInMillis())/ 1000)} ) + " s."); 
        out.flush();

        // Memory cleaning
        fileList = null;
        return true;
    }

    public void run() {

        ArrayList<StemPosition> currentStems = null;
        CSVWriter writer;
        String fileName, fileOut;
        //String auxParam;
        String rnaSequence = ""; //NOI18N
        RNAfold fold = null;
        int i, secuencias = 0;

        try {
            fileName = actualFile.getName();
            out.println(java.text.MessageFormat.format(getBundle("resources/Bundle", currentLocale)
                    .getString("FILE_PRINT"), 
                    new Object[] {fileName})); //$NON-NLS-1$
            out.flush();
            fileName = fileName.replaceAll(FASTA_EXT, ""); //NOI18N
            fileName = fileName.replaceAll(FASTA_EXT_2, ""); //NOI18N
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
            LinkedHashMap<String, DNASequence> fasta = readFastaDNASequence(
                    actualFile, false);

            if (fasta.isEmpty()) {
                out.println(getBundle("resources/Bundle", currentLocale).getString("INVALID_FILE_FORMAT"));
                out.println(getBundle("resources/Bundle", currentLocale).getString("TRYING_TO_FIX"));
                boolean formatFile = FASTACorrector.formatFile(
                        actualFile.getAbsolutePath());
                if (formatFile) {
                    fasta = readFastaDNASequence(actualFile, false);
                } else {
                    out.println(getBundle("resources/Bundle", currentLocale).getString("CANT_PROCESS"));
                    return;
                }
            }
            int listSize = fasta.size();
            i = 1;

            writer = new CSVWriter(new FileWriter(fileOut), ';',
                    CSVWriter.DEFAULT_QUOTE_CHARACTER,
                    CSVWriter.DEFAULT_ESCAPE_CHARACTER,
                    CSVWriter.DEFAULT_LINE_END);
            writer.writeNext(StemLoop.getHeader(this.inputType).split(";")); //NOI18N

            // II. Transcript level
            for (Map.Entry<String, DNASequence> entry : fasta.entrySet()) {

                DNASequence element = entry.getValue();
                StemLoop aux = new StemLoop(this.inputType);
                aux.setTags(element.getOriginalHeader());
                /*out.println("Analizando " + aux.getGeneSymbol() + " - "
                    + aux.getGeneID() + " transcripto " + aux.getTranscriptID() 
                    + "...");*/
                // Here begins one cycle
                Iterator<String> patternItr = loopPatterns.iterator();

                // sólo para modo dos
                
                if(extendedMode){
                    rnaSequence = element.getSequenceAsString().
                        toUpperCase().replace('T', 'U');
                    fold = new RNAfold(rnaSequence);
                }
                //MFEData fold = rfa.getMFE(rnaSequence.getBytes());
                secuencias++;
                
                // III. Loop level
                while (patternItr.hasNext()) {

                    String currentPattern = patternItr.next();
                    currentPattern = currentPattern.trim().toUpperCase();
                    
                    // 1. Stem research
                    if (extendedMode) {
                        searchMode2(rnaSequence,
                                element.getOriginalHeader(),
                                fold,
                                currentPattern,
                                currentStems, writer);
                        //sequenceExtendedResearch(element, currentPattern, currentStems, writer);
                    } else {
                        sequenceExtendedResearch(element, currentPattern, currentStems, writer);
                        //sequenceResearch(element, currentPattern, currentStems, writer);
                    }
                }

                if (i++ % 100 == 0) {
                    out.printf(getBundle("resources/Bundle", currentLocale).getString("PERCENT_PROCESSED_SEQS"),
                            i, listSize, //$NON-NLS-1$
                            (i / (float) listSize) * 100);
                    out.flush();
                }

            }

            writer.close();
            writer = null;
            fasta = null;

            // FASTA with the results
            //auxParam = "Stem loops length between " + this.minLength
            //         + " and " + this.maxLength + ".\n";
            //auxParam += "Wooble pairs allowed up to " + this.maxWooble
            //        + ".\n";
            //auxParam += "Mismatch allowed up to " + this.maxMismatch
            //        + ".\n";
            //auxParam += "Total sequences analized: " + secuencias
            //        + ".\n";
            currentStems = null;
            out.printf(getBundle("resources/Bundle", currentLocale).getString("RESUME"),
                            this.minLength, 
                            this.maxLength, 
                            secuencias);
            out.flush();
        } catch (FileNotFoundException ex) {
            Logger.getLogger(LoopCatcher.class.getName()).log(Level.SEVERE, null, ex);
            out.println(getBundle("resources/Bundle", currentLocale).getString("CANT_OPEN_FILE"));
        } catch (IOException ex) {
            Logger.getLogger(LoopCatcher.class.getName()).log(Level.SEVERE, null, ex);
        }
        currentStems = null;
        writer = null;
        fileName = null;
        fileOut = null;
        //auxParam = null;
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
