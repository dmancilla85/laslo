/**
 *
 */
package com.laslo.core;

import com.tools.ShuffleSeq;
import static com.laslo.core.PairmentAnalizer.*;
import com.tools.io.InputSequence;
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
import com.tools.UShuffle;
import com.tools.io.FASTACorrector;
import com.tools.io.SourceFile;
import static com.tools.io.SourceFile.CSV_EXT;
import static com.tools.io.SourceFile.FASTA_EXT;
import static com.tools.io.SourceFile.FASTA_EXT_2;
import static com.tools.io.SourceFile.GENBANK_EXT;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Locale;
import java.util.Map;
import java.util.ResourceBundle;
import static java.util.ResourceBundle.getBundle;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.lang.StringUtils;
import org.biojava.nbio.core.sequence.DNASequence;
import static org.biojava.nbio.core.sequence.io.FastaReaderHelper.readFastaDNASequence;
import static org.biojava.nbio.core.sequence.io.GenbankReaderHelper.readGenbankDNASequence;
//import rnafold4j.MFEData;
//import rnafold4j.RNAFoldAPI;

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

    public int getkLetRandoms() {
        return kLetRandoms;
    }

    public void setkLetRandoms(int kLetRandoms) {
        this.kLetRandoms = kLetRandoms;
    }

    public String getAdditionalSequence() {
        return additionalSequence;
    }

    public void setAdditionalSequence(String additionalSequence) {
        this.additionalSequence = additionalSequence;
    }

    // Create a new instance of RNAFold4J
    //protected RNAFoldAPI rfa;
    public LoopCatcher(String pathOut, String pathIn,
            ArrayList<String> loopPatterns, String additionalSequence,
            InputSequence inputType,
            int minLength, int maxLength,
            int maxWooble, int maxMismatch, Locale locale) {
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
        this.bundle = getBundle("resources/Bundle", locale);
    }

    /**
     *
     */
    public LoopCatcher() {
        this("", "", new ArrayList<>(), BiologyPatterns.PUM1,
                InputSequence.ENSEMBL, //NOI18N
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
     *
     * @param fastaSeq
     * @param stemLoopPattern
     * @param writer
     * @return
     */
    public int sequenceResearch(DNASequence fastaSeq,
            String stemLoopPattern,  CSVWriter writer) {

        List<StemLoop> slrList = new ArrayList<>();
        //List<Integer> mismatchs = new ArrayList<>();
        StemLoop slr;
        slr = null;
        int size, posAux, k=3;
        size = 0;
        int loopPos = 0, loopLength = stemLoopPattern.length();
        boolean isValidHairpin;
        String rnaSequence = fastaSeq.getRNASequence().getSequenceAsString();
        final int sequenceLength = rnaSequence.length();
        //MFEData mfe = null;
        RNAfold fold = null;

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
            int length = this.maxLength;
            slr = new StemLoop(this.inputType);
            slr.setTags(fastaSeq.getOriginalHeader());
            try {
                if ((loopPos - length) > 0
                        && (loopPos + loopLength + length) < sequenceLength) {
                    rnaSeq = rnaSequence.substring(loopPos - length,
                            loopPos + loopLength + length);

                    //mfe = rfa.getMFE(rnaSeq.getBytes());
                    fold = new RNAfold(rnaSeq);

                    hairpinModel = //new String(mfe.structure);
                            fold.getStructure();

                    if (rnaSeq.length() != hairpinModel.length()) {
                        out.println("Error NO COINCIDEN:" + rnaSeq + " - "
                                + hairpinModel);
                    }

                    if (fold.getMfe()/*mfe.mfe*/ == 0.0) {
                        isValidHairpin = false;
                    } else {
                        hairpinModel = isValidHairpin(hairpinModel, loopLength, loopPos, rnaSeq);
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

                posAux = rnaSequence.indexOf(rnaSeq);
                
                // Fill the fields
                slr.setRnaHairpinSequence(rnaSeq);
                slr.setLoop(rnaLoop);
                slr.setStartsAt(loopPos - extIzq);
                slr.setStructure(hairpinModel);
                slr.setSequenceLength(sequenceLength);
                try{
                slr.setAdditional5Seq( rnaSequence
                        .substring(posAux - 1 -k, posAux - 1));
                slr.setAdditional5Seq( rnaSequence.substring(posAux 
                        + rnaSeq.length() + 1,posAux + rnaSeq.length() + 1 +k) );
                } catch(IndexOutOfBoundsException e){
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
                // ----------------------
                // PUM1 Protein (pumilio)
                slr.setAdditionalSeqLocations(getPatternLocations(rnaSequence,
                        this.additionalSequence));
                // Polyadenilations Sites
                // Look for the Kozak pattern
                // Special TEST buscar otros posibles stems
                slr.setRelativePos((double) slr.getStartsAt()
                        / (double) rnaSequence.length());

                slrList.add(slr);
            }

        }

        slr = null;
        //mismatchs = null;

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
                    

                    if (rnaSeq.length() != hairpinModel.length()) {
                        out.println("Error NO COINCIDEN:" + rnaSeq + " - "
                                + hairpinModel);
                    }
                    
                    hairpinModel = isValidHairpin(hairpinModel, loopLength, loopPos, rnaSeq);
                    isValidHairpin = hairpinModel.length() > 0;

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
                try{
                slr.setAdditional5Seq( rnaSequence
                        .substring(posAux - 1 -k, posAux - 1));
                slr.setAdditional5Seq( rnaSequence.substring(posAux 
                        + rnaSeq.length() + 1,posAux + rnaSeq.length() + 1 +k) );
                } catch(IndexOutOfBoundsException e){
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
        ini = Calendar.getInstance();
        out.println(
                java.text.MessageFormat.format(bundle.getString("START_TIME"), new Object[]{Calendar.getInstance().getTime()}));
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
                        if(currentFile.toString().endsWith(GENBANK_EXT))
                            fasta = readGenbankDNASequence(currentFile, false);
                        else
                            fasta = readFastaDNASequence(currentFile, false);
                        
                    } catch (IOException ex ) {
                        Logger.getLogger(LoopCatcher.class.getName()).log(Level.SEVERE, null, ex);
                    } catch (Exception ex) {
                        Logger.getLogger(LoopCatcher.class.getName()).log(Level.SEVERE, null, ex);
                    }

                    UShuffle.makeShuffleSequences(pathOut, currentFile.getName(),
                            fasta, numberOfRandoms, kLetRandoms);
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
                    || currentFile.toString().endsWith(FASTA_EXT_2)
                    || currentFile.toString().endsWith(GENBANK_EXT))) {

                this.actualFile = currentFile;
                run();
            }
        }

        fin = Calendar.getInstance();
        out.println(java.text.MessageFormat.format(bundle.getString("TOTAL_TIME"), new Object[]{((fin.getTimeInMillis() - ini.getTimeInMillis()) / 1000)}) + " s.");

        out.flush();

        // Memory cleaning
        fileList = null;
        return true;
    }

    public void run() {

        //ArrayList<StemPosition> currentStems = null;
        CSVWriter writer;
        String fileName, fileOut;
        //String auxParam;
        String rnaSequence = ""; //NOI18N
        RNAfold fold = new RNAfold();
        int i, secuencias = 0;

        try {
            fileName = actualFile.getName();
            out.println(java.text.MessageFormat.format(bundle.getString("FILE_PRINT"), new Object[]{fileName})); //$NON-NLS-1$
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
            
            if(actualFile.getName().endsWith(GENBANK_EXT))
                fasta = readGenbankDNASequence(actualFile, false);
            else
                fasta = readFastaDNASequence(actualFile, false);

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

            this.inputType = SourceFile.detectHeader(fasta.entrySet().iterator()
                    .next().getValue().getOriginalHeader());
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

                // Here begins one cycle
                Iterator<String> patternItr = loopPatterns.iterator();

                // sólo para modo dos
                if (extendedMode) {
                    rnaSequence = element.getRNASequence().getSequenceAsString();
                    fold = new RNAfold(rnaSequence);
                }
                secuencias++;

                // III. Loop level
                while (patternItr.hasNext()) {

                    String currentPattern = patternItr.next();
                    currentPattern = currentPattern.trim().toUpperCase();

                    // 1. Stem research
                    if (extendedMode) {
                        sequenceExtendedResearch(rnaSequence,
                                element.getOriginalHeader(),
                                fold.getStructure(),
                                currentPattern, writer);
                    } else {
                        sequenceResearch(element, currentPattern, writer);
                    }
                }

                if (i++ % 100 == 0) {
                    out.printf(bundle.getString("PERCENT_PROCESSED_SEQS"),
                            i, listSize, //$NON-NLS-1$
                            (i / (float) listSize) * 100);
                    out.flush();
                }

            }

            writer.close();
            writer = null;
            fasta = null;

            out.printf(bundle.getString("RESUME"),
                    this.minLength,
                    this.maxLength,
                    secuencias);
            out.flush();
        } catch (FileNotFoundException ex) {
            Logger.getLogger(LoopCatcher.class.getName())
                    .log(Level.SEVERE, null, ex);
            out.println(bundle.getString("CANT_OPEN_FILE"));
        } catch (IOException ex) {
            Logger.getLogger(LoopCatcher.class.getName())
                    .log(Level.SEVERE, null, ex);
        } catch (Exception ex) {
            Logger.getLogger(LoopCatcher.class.getName())
                    .log(Level.SEVERE, null, ex);
        }
        
        writer = null;
        fileName = null;
        fileOut = null;
        
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
