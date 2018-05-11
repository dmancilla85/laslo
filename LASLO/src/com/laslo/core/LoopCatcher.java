/**
 *
 */
package com.laslo.core;

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
import com.tools.joewandy.Sequence;
import com.tools.io.FastaReader;
import com.tools.io.FastaReader.FastaFileContent;
import com.opencsv.CSVWriter;
import java.util.logging.Level;
import java.util.logging.Logger;
import rnafold4j.MFEData;
import rnafold4j.RNAFoldAPI;

/**
 * @author David
 *
 */
public class LoopCatcher implements Runnable {

    protected String pathOut;
    protected ArrayList<String> loopPatterns;
    protected InputSequence inputType;
    protected int minLength;
    protected int maxLength;
    protected int maxWooble;
    protected int maxMismatch;
    protected File[] fileList;
    protected boolean extendedMode;
    protected static final String LOG_EXT = ".log";
    protected static final String CSV_EXT = ".csv";
    protected static final String FASTA_EXT = ".fasta";
    protected static final String FASTA_EXT_2 = ".fa";
    // Create a new instance of RNAFold4J
    protected RNAFoldAPI rfa;

    public LoopCatcher(String pathOut, String pathIn,
            ArrayList<String> loopPatterns, InputSequence inputType,
            int minLength, int maxLength,
            int maxWooble, int maxMismatch) {
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
    }

    public LoopCatcher() {
        this.pathOut = "";
        this.loopPatterns = new ArrayList<>();
        this.inputType = InputSequence.ENSEMBL;
        this.minLength = 4;
        this.maxLength = 16;
        this.maxWooble = 2;
        this.maxMismatch = 0;
        this.rfa = new RNAFoldAPI();
    }

    public File[] getFileList() {
        return fileList;
    }

    public void setFileList(File[] fileList) {
        this.fileList = fileList;
    }

    public void setIsExtendedMode(boolean mode) {
        this.extendedMode = mode;
    }

    public boolean getIsExtendedMode() {
        return this.extendedMode;
    }

    public String getPathOut() {
        return pathOut;
    }

    public void setPathOut(String pathOut) {
        this.pathOut = pathOut;
    }

    public ArrayList<String> getLoopPatterns() {
        return loopPatterns;
    }

    public void setLoopPatterns(ArrayList<String> loopPatterns) {
        this.loopPatterns = loopPatterns;
    }

    public InputSequence getInputType() {
        return inputType;
    }

    public void setInputType(InputSequence inputType) {
        this.inputType = inputType;
    }

    public int getMinLength() {
        return minLength;
    }

    public void setMinLength(int minLength) {
        this.minLength = minLength;
    }

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

    public int getMaxMismatch() {
        return maxMismatch;
    }

    public void setMaxMismatch(int maxMismatch) {
        this.maxMismatch = maxMismatch;
    }

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
                header = headerList.split(";");
                writer.writeNext(header);
            }

            for (int i = 0; i < stemResearch.size(); i++) {
                String[] data;
                data = stemResearch.get(i).toRowCSV().split(";");
                writer.writeNext(data);
                //allData.add(data);
            }

            writer.close();

        } catch (IOException e) {
            // TODO Auto-generated catch block

        }
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
     * @param fastaSeq
     * @param stemLoopPattern
     * @param stemsDetected
     * @param writer
     * @return
     */
    public int sequenceResearch(Sequence fastaSeq, String stemLoopPattern,
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
                loopLength = stemLoopPattern.length(),
                sequenceLength = fastaSeq.getSequenceString().length(), pos1,
                pos2;
        boolean existsInvertedRepeat;
        existsInvertedRepeat = false;
        boolean checkWooblePair, isValidHairpin;
        // Upper all characters
        String dnaSequence = fastaSeq.getSequenceString().toUpperCase();

        // Make the retro transcription
        String rnaSequence = dnaSequence.replace('T', 'U');

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
                                || (isComplementaryRNAWooble(rnaSequence.charAt(pos1),
                                        rnaSequence.charAt(pos2)) && checkWooblePair)
                                || mismatchCount < this.maxMismatch) {
                            existsInvertedRepeat = true;

                            if (!isComplementaryRNA(
                                    rnaSequence.charAt(pos1), rnaSequence.charAt(pos2))) {
                                // A mismatch shouldn't be considered if is at
                                // the border

                                if (isComplementaryRNAWooble(rnaSequence.charAt(pos1),
                                        rnaSequence.charAt(pos2))
                                        && woobleCount < this.maxWooble) {
                                    woobleCount++;

                                } else if (lengthIR > 0 && mismatchCount
                                        < this.maxMismatch) {
                                    mismatchCount++;
                                    mismatchs.add(lengthIR + 1);
                                    mismatchs.add(loopPos + loopLength + lengthIR);
                                } else {
                                    existsInvertedRepeat = false;
                                }
                            }

                            lengthIR++;
                        } else {
                            existsInvertedRepeat = false;
                        }

                    } catch (IndexOutOfBoundsException e) {
                        out.println("Error: " + dnaSequence + " - i: "
                                + (loopPos - lengthIR) + " - i2: " //$NON-NLS-3$ //$NON-NLS-2$
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
            slr.setTags(fastaSeq.getId());
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

                slr.setLoop(rnaSequence.substring(loopPos, loopPos + loopLength));
                hairpinModel = slr.drawHairpinViennaStructure();

                if (mfe.mfe == 0 || !hairpinModel.equals(new String(mfe.structure).trim())) {
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
            writer.writeNext(element.toRowCSV().split(";"));
        }

        size = slrList.size();
        slrList.clear();
        slrList = null;

        return size;
    }

    public int sequenceExtendedResearch(Sequence fastaSeq, String stemLoopPattern,
            ArrayList<StemPosition> stemsDetected, CSVWriter writer) {

        List<StemLoop> slrList = new ArrayList<>();
        StemLoop slr;
        slr = null;
        List<Integer> mismatchs = new ArrayList<>();
        int size, i;
        char aux1, aux2;
        size = 0;

        int loopPos = 0,
                loopLength = stemLoopPattern.length(),
                sequenceLength = fastaSeq.getSequenceString().length();

        boolean isValidHairpin;
        String dnaSequence = fastaSeq.getSequenceString().toUpperCase();
        String pairs = "";
        String rnaSequence = dnaSequence.replace('T', 'U');
        MFEData mfe = null;

        // Convert the original loop pattern to a regular expression
        String regExp = toRegularExpression(stemLoopPattern);
        Pattern p = Pattern.compile(regExp);
        Matcher loopFinder = p.matcher(rnaSequence);

        String rnaLoop = "", rnaSeq = "", hairpinModel = "";

        // As exists loop matches
        while (loopFinder.find()) {
            loopPos = loopFinder.start();
            mismatchs.clear();
            isValidHairpin = true;
            rnaLoop = rnaSequence.substring(loopPos, loopFinder.end());
            int length = this.maxLength;
            slr = new StemLoop(this.inputType);
            slr.setTags(fastaSeq.getId());
            try {
                if ((loopPos - length) > 0
                        && (loopPos + loopLength + length) < sequenceLength) {
                    rnaSeq = rnaSequence.substring(loopPos - length, loopPos
                            + loopLength + length);

                    mfe = rfa.getMFE(rnaSeq.getBytes());
                    hairpinModel = new String(mfe.structure).trim();

                    if (mfe.mfe == 0) {
                        isValidHairpin = false;
                    }

                    if (isValidHairpin) {
                        boolean cut = true;

                        for (i = 0; i < length && isValidHairpin && cut; i++) {
                            if (hairpinModel.charAt(hairpinModel.length() - 1 - i) == '.'
                                    && hairpinModel.charAt(i) == '.') {

                                if (length - i - 1 > this.minLength) {
                                    isValidHairpin = false;
                                }
                            } else {
                                cut = false;
                            }
                        }

                        if (isValidHairpin && i > 0) {
                            hairpinModel = hairpinModel.substring(i, hairpinModel.length() - i);
                            rnaSeq = rnaSeq.substring(i, rnaSeq.length() - i);
                            length = length - i;
                        }

                    }

                    if (isValidHairpin) {
                        slr.setRnaHairpinSequence(rnaSeq);
                        slr.setLoop(rnaLoop);
                        slr.setStartsAt(loopPos - length);
                        //slr.checkPairments();
                        /*isValidHairpin = hairpinModel.length() 
                            - hairpinModel.replace("(", "").length() >= length * 0.5;*/
                    }

                    if (isValidHairpin) {
                        for (i = length; i < length + loopLength && isValidHairpin; i++) {
                            aux1 = hairpinModel.charAt(i);

                            if (aux1 != '.') {
                                isValidHairpin = false;
                            }
                        }
                    }

                    if (isValidHairpin) {
                        isValidHairpin = hairpinModel.charAt(0) == '('
                                || hairpinModel.charAt(hairpinModel.length() - 1) == ')';
                    }
                    if (isValidHairpin) {
                        aux1 = hairpinModel.charAt(length - 1);
                        aux2 = hairpinModel.charAt(length + loopLength);
                        isValidHairpin = aux1 == '(' && aux2 == ')';
                    }

                } else {
                    isValidHairpin = false;
                }

            } catch (Exception e) {
                out.println("ERROR!: Length" + length + " - Hairpin: " + hairpinModel + " - seq: " + rnaSeq);
            }

            if (isValidHairpin) {
                // Fill the fields
                slr.setSequenceLength(sequenceLength);
                slr.setMfe(mfe.mfe);
                slr.setStructure(hairpinModel);
                slr.setPredecessorLoop(length);
                slr.setNLoop(length);
                slr.setPercent_AG();
                slr.setLoopPattern(stemLoopPattern);
                slr.setPercent_AU();
                slr.setPercent_CG();
                slr.setStartsAt(loopPos - length + 1);
                slr.setEndsAt(loopFinder.end() + length + 1);
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
            writer.writeNext(element.toRowCSV().split(";"));
        }

        size = slrList.size();
        slrList.clear();
        slrList = null;

        return size;
    }

    @Override
    public void run() {
        ArrayList<StemPosition> currentStems = null;
        // Auxiliar values
        int i;
        FastaReader fasta = null;
        List<Sequence> list = null;
        Iterator<Sequence> itr = null;
        CSVWriter writer;
        String auxParam;
        boolean delete;
        // To check the elapsed time
        Calendar ini, fin;
        String fileOut, fileName;

        ini = Calendar.getInstance();
        out.println("Process started at " + Calendar.getInstance().getTime());
        out.flush();

        if (fileList == null) {
            return;
        }

        // I. File level
        for (File listOfFile : fileList) {
            if (listOfFile.isFile() && (listOfFile.toString().endsWith(FASTA_EXT)
                    || listOfFile.toString().endsWith(FASTA_EXT_2))) {
                try {
                    fileName = listOfFile.getName();
                    out.println("File " + fileName); //$NON-NLS-1$
                    out.flush();
                    fileName = fileName.replaceAll(FASTA_EXT, "");
                    fileName = fileName.replaceAll(FASTA_EXT_2, "");
                    fileOut = this.pathOut + "\\" + fileName + CSV_EXT;

                    if (new File(fileOut).exists()) {
                        try {
                            delete = (new File(fileOut)).delete();
                        } catch (Exception io) {
                            out.println(io.getMessage());
                            out.flush();
                        }
                    }

                    Calendar.getInstance();
                    // Generation of the iterator of {id,sequence}
                    fasta = new FastaReader(FastaFileContent.DNA_SEQUENCE, listOfFile.getAbsolutePath());
                    list = fasta.processFile();
                    if (list.isEmpty()) {
                        out.println("No hay archivos para procesar.");
                        System.exit(0);
                    }
                    itr = list.iterator();
                    int listSize = list.size();
                    i = 1;

                    writer = new CSVWriter(new FileWriter(fileOut), ';',
                            CSVWriter.DEFAULT_QUOTE_CHARACTER,
                            CSVWriter.DEFAULT_ESCAPE_CHARACTER,
                            CSVWriter.DEFAULT_LINE_END);
                    writer.writeNext(StemLoop.getHeader(this.inputType).split(";"));

                    // II. Transcript level
                    while (itr.hasNext()) {
                        Sequence element = itr.next();
                        StemLoop aux = new StemLoop(this.inputType);
                        aux.setTags(element.getId());

                        // Here begins one cycle
                        Iterator<String> patternItr = loopPatterns.iterator();

                        while (patternItr.hasNext()) {

                            String currentPattern = patternItr.next();
                            currentPattern = currentPattern.trim();

                            // 1. Stem research
                            if (extendedMode) {
                                sequenceExtendedResearch(element, currentPattern,
                                        currentStems, writer);
                            } else {
                                sequenceResearch(element, currentPattern, currentStems, writer);
                            }
                        }

                        if (i++ % 3000 == 0) {
                            out.printf("%d out of %d. %2.2f%% processed...%n",
                                    i, listSize, //$NON-NLS-1$
                                    (i / (float) listSize) * 100);
                            out.flush();
                        }

                    }

                    writer.close();
                    writer = null;
                    // FASTA with the results
                    auxParam = "Stem loops length between " + this.minLength
                            + " and " + this.maxLength + ".\n";
                    auxParam += "Wooble pairs allowed up to " + this.maxWooble
                            + ".\n";
                    auxParam += "Mismatch allowed up to " + this.maxMismatch
                            + ".\n";
                    currentStems = null;
                    out.println(auxParam);
                    out.flush();
                } catch (IOException ex) {
                    Logger.getLogger(LoopCatcher.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        }

        fin = Calendar.getInstance();
        out.println("Tiempo: " + ((fin.getTimeInMillis() - ini.getTimeInMillis())
                / 1000) + " segundos."); //$NON-NLS-2$ //$NON-NLS-2$
        out.flush();

        // Memory cleaning
        list = null;
        itr = null;
        fasta = null;
        fileList = null;
    }
}
