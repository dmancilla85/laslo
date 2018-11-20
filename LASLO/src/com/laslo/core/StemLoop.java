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

import static com.laslo.core.SequenceAnalizer.reverseIt;
import com.tools.io.InputSequence;
import com.tools.io.BioMartFastaID;
import com.tools.io.GenericID;
import com.tools.io.SourceFile;
import com.tools.io.FlyBaseFastaID;
import com.tools.io.EnsemblFastaID;
import com.tools.io.GenBankID;
import static java.lang.System.out;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;
import org.apache.commons.lang.StringUtils;

class PosComparator implements Comparator<Integer>{
    
    @Override
    public int compare(Integer a, Integer b){
        return a > b ? -1 :(a < b ? 1 : 0);
    }
}

/**
 * 13/12/2016
 *
 * @author David A. Mancilla
 */
public class StemLoop {

    protected SourceFile id_fasta;
    protected InputSequence mode;
    protected String rnaHairpinSequence;
    protected String loop;
    protected String hairpinStructure;
    protected String additional5Seq;
    protected String additional3Seq;
    protected String loopPattern;
    protected String viennaStructure;
    protected int sequenceLength;
    protected int startsAt;
    protected int endsAt;
    protected int bulge;
    protected int internalLoops;
    protected char predecessor2Loop;
    protected char predecessorLoop;
    protected char n2Loop;
    protected char n5Loop;
    protected char n6Loop;
    protected char n7Loop;
    protected char n8Loop;
    protected double percent_AG;
    protected double percent_GU;
    protected double percent_CG;
    protected double percent_AU;
    protected double percA_sequence;
    protected double percG_sequence;
    protected double percC_sequence;
    protected double percU_sequence;
    protected double relativePos;
    protected boolean reversed;
    protected double mfe;
    protected List<Integer> additionalSeqLocations;

    /**
     *
     * @return
     */
    public String getAdditional5Seq() {
        return additional5Seq;
    }

    /**
     *
     * @param additional5Seq
     */
    public void setAdditional5Seq(String additional5Seq) {
        this.additional5Seq = additional5Seq;
    }

    /**
     *
     * @return
     */
    public String getAdditional3Seq() {
        return additional3Seq;
    }

    /**
     *
     * @param additional3Seq
     */
    public void setAdditional3Seq(String additional3Seq) {
        this.additional3Seq = additional3Seq;
    }

    /**
     *
     * @param number
     * @return
     */
    private String getFormattedNumber(Double number, int digits) {
        NumberFormat numberFormatter
                = NumberFormat.getNumberInstance(Locale.getDefault());
        numberFormatter.setMinimumFractionDigits(digits);
        return numberFormatter.format(number);
    }
    
    private String getFormattedNumber(Double number) {
        return getFormattedNumber(number,3);
    }

    /**
     *
     * @param mode
     */
    public StemLoop(InputSequence mode) {
        super();
        this.mode = mode;

        switch (mode) {
            case ENSEMBL:
                this.id_fasta = new EnsemblFastaID();
                break;

            case FLYBASE:
                this.id_fasta = new FlyBaseFastaID();
                break;

            case BIOMART:
                this.id_fasta = new BioMartFastaID();
                break;

            case GENERIC:
                this.id_fasta = new GenericID();
                break;
                
            case GENBANK:
                this.id_fasta = new GenBankID();
                break;

        }

        //this.stemLoopId = 0;
        this.loop = ""; //$NON-NLS-1$
        this.rnaHairpinSequence = null;
        this.sequenceLength = 0;
        this.percent_AG = 0;
        this.percent_GU = 0;
        this.percent_CG = 0;
        this.percent_AU = 0;
        this.percA_sequence = 0;
        this.percC_sequence = 0;
        this.percG_sequence = 0;
        this.percU_sequence = 0;
        this.loopPattern = "";
        this.endsAt = 0;
        this.startsAt = 0;
        this.bulge = 0;
        this.reversed = false;
        this.internalLoops = 0;
        this.relativePos = 0.0;
        this.predecessorLoop = 0;
        this.predecessor2Loop = 0;
        this.n2Loop = this.n5Loop = this.n6Loop
                = this.n7Loop = this.n8Loop = 0;

        this.mfe = (float) 0.0;
        this.viennaStructure = "";
        this.additional5Seq = "";
        this.additional3Seq = "";
        this.additionalSeqLocations = new ArrayList<>();
    }

    /**
     *
     * @return
     */
    public int getMismatches() {
        return bulge;
    }

    /**
     *
     * @param invert
     */
    public void setReverse(boolean invert) {
        this.reversed = invert;
    }

    /**
     *
     * @return
     */
    public boolean getReversed() {
        return reversed;
    }

    /**
     *
     * @param mismatches
     */
    public void setMismatches(int mismatches) {
        this.bulge = mismatches;
    }

    /**
     *
     * @param pattern
     */
    public void setLoopPattern(String pattern) {
        this.loopPattern = pattern;
    }

    /**
     *
     * @return
     */
    public double getRelativePos() {
        return this.relativePos;
    }

    /**
     *
     * @param relativePos
     */
    public void setRelativePos(double relativePos) {
        this.relativePos = relativePos;
    }

    /**
     *
     * @param mode
     * @return
     */
    public static String getHeader(InputSequence mode) {

        String header = "";

        switch (mode) {
            case ENSEMBL:
                header = EnsemblFastaID.getHeader();
                break;

            case FLYBASE:
                header = FlyBaseFastaID.getHeader();
                break;

            case BIOMART:
                header = BioMartFastaID.getHeader();
                break;

            case GENBANK:
                header = GenBankID.getHeader();
                break;

            case GENERIC:
                header = GenericID.getHeader();
                break;
        }

        return header
                + "LoopPattern" + SourceFile.ROW_DELIMITER
                + "LoopID" + SourceFile.ROW_DELIMITER
                + "TerminalPair" + SourceFile.ROW_DELIMITER
                + "Sense" + SourceFile.ROW_DELIMITER
                + "N-2" + SourceFile.ROW_DELIMITER
                + "N-1" + SourceFile.ROW_DELIMITER
                + "N2" + SourceFile.ROW_DELIMITER
                + "N5" + SourceFile.ROW_DELIMITER
                + "N6" + SourceFile.ROW_DELIMITER
                + "N7" + SourceFile.ROW_DELIMITER
                + "N8" + SourceFile.ROW_DELIMITER
                + "Loop" + SourceFile.ROW_DELIMITER
                + "StemLoopSequence" + SourceFile.ROW_DELIMITER
                + "Additional5Seq" + SourceFile.ROW_DELIMITER
                + "Additional3Seq" + SourceFile.ROW_DELIMITER
                + "PredictedStructure" + SourceFile.ROW_DELIMITER
                + "ViennaBracketStr" + SourceFile.ROW_DELIMITER
                + "Pairments" + SourceFile.ROW_DELIMITER
                + "WooblePairs" + SourceFile.ROW_DELIMITER
                + "Bulges" + SourceFile.ROW_DELIMITER
                + "InternalLoops" + SourceFile.ROW_DELIMITER
                + "SequenceLength" + SourceFile.ROW_DELIMITER
                + "StartsAt" + SourceFile.ROW_DELIMITER
                + "EndsAt" + SourceFile.ROW_DELIMITER
                + "A_PercentSequence" + SourceFile.ROW_DELIMITER
                + "C_PercentSequence" + SourceFile.ROW_DELIMITER
                + "G_PercentSequence" + SourceFile.ROW_DELIMITER
                + "U_PercentSequence" + SourceFile.ROW_DELIMITER
                + "AU_PercentPairs" + SourceFile.ROW_DELIMITER
                + "CG_PercentPairs" + SourceFile.ROW_DELIMITER
                + "GU_PercentPairs" + SourceFile.ROW_DELIMITER
                + "PurinePercentStem" + SourceFile.ROW_DELIMITER
                + "RnaFoldMFE" + SourceFile.ROW_DELIMITER
                + "RelativePosition" + SourceFile.ROW_DELIMITER
                + "AdditionalSeqMatches" + SourceFile.ROW_DELIMITER
                + "AdditionalSeqPositions" + SourceFile.ROW_DELIMITER;
    }

    /**
     *
     * @return
     */
    public int getEndsAt() {
        return endsAt;
    }

    public String isReverse(){
        if(!this.reversed)
            return "+";
        else
            return "-";
    }
    
    /**
     *
     * @return
     */
    public String getGeneID() {
        return this.id_fasta.getGeneID();
    }

    /**
     *
     * @return
     */
    public String getGeneSymbol() {
        if (this.id_fasta.getGeneSymbol() != null) {
            return this.id_fasta.getGeneSymbol();
        } else {
            return "";
        }
    }

    /**
     *
     * @return
     */
    public String getGUPairs() {
        long pairs = Math.round(this.percent_GU * this.getPairments());

        return pairs + "";
    }

    /**
     *
     * @return
     */
    public String getHairpinStructure() {
        return hairpinStructure; //drawHairpinStructure();
    }

    /**
     *
     * @return
     */
    public String getLoop() {
        
        String theLoop = this.loop;
        
        if (reversed) {
            theLoop = reverseIt(theLoop);
        }
        
        return theLoop;
    }

    /**
     *
     * @return
     */
    public String getLoopID() {
        String theLoop, terminal;

        theLoop = this.loop;
        terminal = this.getTerminalPair();

        if (reversed) {
            theLoop = reverseIt(theLoop);
            terminal = reverseIt(terminal);
        }

        return this.predecessorLoop + theLoop
                + "(" + terminal + ")|" + this.getPairments();
    }

    /**
     *
     * @return
     */
    public String getLoopPattern() {

        return this.loopPattern;
    }

    /**
     *
     * @return
     */
    public InputSequence getMode() {
        return mode;
    }

    /**
     *
     * @return
     */
    public char getN2Loop() {
        return n2Loop;
    }

    /**
     *
     * @return
     */
    public char getN5Loop() {
        return n5Loop;
    }

    /**
     *
     * @return
     */
    public char getN6Loop() {
        return n6Loop;
    }

    /**
     *
     * @return
     */
    public char getN7Loop() {
        return n7Loop;
    }

    /**
     *
     * @return
     */
    public char getN8Loop() {
        return n8Loop;
    }

    /**
     *
     * @return
     */
    public double getPercA_sequence() {
        return percA_sequence;
    }

    /**
     *
     * @return
     */
    public double getPercC_sequence() {
        return percC_sequence;
    }

    /**
     *
     * @return
     */
    public double getPercent_AG() {
        return percent_AG;
    }

    /**
     *
     * @return
     */
    public double getPercent_AU() {
        return percent_AU;
    }

    /**
     *
     * @return
     */
    public double getPercent_CG() {
        return percent_CG;
    }

    /**
     *
     * @return
     */
    public double getPercent_GU() {
        return percent_GU;
    }

    /**
     *
     * @return
     */
    public double getPercG_sequence() {
        return percG_sequence;
    }

    /**
     *
     * @return
     */
    public double getPercU_sequence() {
        return percU_sequence;
    }

    /**
     *
     * @return
     */
    public char getPredecessorLoop() {
        return predecessorLoop;
    }

    /**
     *
     * @return
     */
    public String getAdditionalSequenceCount() {

        Integer count = 0;

        if (this.additionalSeqLocations != null) {
            count = this.additionalSeqLocations.size();
        }

        return count.toString();
    }

    /**
     *
     * @return
     */
    public String getAdditionalSequenceLocations() {

        String locations = ""; //$NON-NLS-1$

        if (this.additionalSeqLocations == null) {
            return "";
        }

        Iterator<Integer> itr = this.additionalSeqLocations.iterator();

        while (itr.hasNext()) {

            Integer element = itr.next();

            if (itr.hasNext()) {
                locations = locations.concat(element + ","); //$NON-NLS-1$
            } else {
                locations = locations.concat(Integer.toString(element));
            }
        }

        return locations;
    }

    /**
     *
     * @return
     */
    public String getRnaHairpinSequence() {
        return rnaHairpinSequence;
    }

    /**
     *
     * @return
     */
    public int getSequenceLength() {
        return sequenceLength;
    }

    /**
     *
     * @return
     */
    public int getStartsAt() {
        return startsAt;
    }

    /**
     *
     * @return
     */
    public int getPairments() {
        int count = StringUtils.countMatches(viennaStructure, "(");
        return count;
    }

    /**
     *
     * @return
     */
    public String getTerminalPair() {
        Character a = rnaHairpinSequence.charAt(viennaStructure.lastIndexOf("("));
        Character b = rnaHairpinSequence.charAt(viennaStructure.indexOf(")"));

        return a.toString() + b.toString();
    }

    /**
     *
     * @return
     */
    public String getTranscriptID() {
        return this.id_fasta.getTranscriptID();
    }

    /**
     *
     * @param endsAt
     */
    public void setEndsAt(int endsAt) {
        this.endsAt = endsAt;
    }

    /**
     *
     * @param loopPattern
     */
    public void setLoop(String loopPattern) {
        this.loop = loopPattern;
    }

    public void setLocation(int pos){
       ((GenBankID) id_fasta).setLocation(pos);
    }
    
    /**
     *
     */
    public void checkPairments() {

        String seq = this.rnaHairpinSequence;
        int woobleCount = 0;
        int CG = 0, AU = 0;
        StringBuilder aux = new StringBuilder(this.viennaStructure);
        int firstIzq = 1;
        
        if (this.loop.isEmpty() || this.viennaStructure.isEmpty()) {
            return;
        }

        // Count internal loops and internalLoops
        /*mismatch = StringUtils.countMatches(this.viennaStructure, "(.(")
                + StringUtils.countMatches(this.viennaStructure, ").)");

        bulge = StringUtils.countMatches(this.viennaStructure, "..(")
                + StringUtils.countMatches(this.viennaStructure, ")..");*/

        try {
            firstIzq = this.viennaStructure.lastIndexOf('(');
            int firstDer = this.viennaStructure.indexOf(')');

            for (int i = firstIzq; i >= 0 && firstDer < viennaStructure.length(); i--) {
                if (viennaStructure.charAt(i) == '(') {
                    while (viennaStructure.charAt(firstDer) != ')') {
                        firstDer++;
                    }

                    if (viennaStructure.charAt(firstDer) == ')') {
                        if (SequenceAnalizer.isComplementaryRNAWooble(seq.charAt(i), seq.charAt(firstDer))) {
                            aux.setCharAt(i, '{');
                            aux.setCharAt(firstDer, '}');
                            woobleCount++;
                        } else {
                            if ((seq.charAt(i) == 'U' && seq.charAt(firstDer) == 'A')
                                    || (seq.charAt(i) == 'A' && seq.charAt(firstDer) == 'U')) {
                                AU++;
                            }

                            if ((seq.charAt(i) == 'C' && seq.charAt(firstDer) == 'G')
                                    || (seq.charAt(i) == 'G' && seq.charAt(firstDer) == 'C')) {
                                CG++;
                            }
                        }
                    }

                    firstDer++;
                }
            }

        } catch (Exception e) {
            System.out.println("checkPairments-ERROR: " + e.getMessage());
        }
        this.hairpinStructure = aux.toString();
        this.percent_AU = (double) AU / (double) (AU + CG + woobleCount);
        this.percent_CG = (double) CG / (double) (AU + CG + woobleCount);
        this.percent_GU = (double) woobleCount / (double) (AU + CG + woobleCount);
        //this.setMismatches(mismatch);
        //this.internalLoops = bulge;

    }
    
    /*public static void main(String []args){
       String a =  "((.((...(((.........)))...)).))";
       String b =   "((((.(((((((.((....)).)))))))))))";
       String c =   "((((((((...((((....))))...))))))))";
       String d =   "(((((((..(((((.....))))))))))))";
       String e =   "(((((((.....)))).)))";

       StemLoop sl = new StemLoop(InputSequence.BIOMART);
       sl.viennaStructure = e;
       sl.checkInternalLoops();
       
       out.println("Internal: " + sl.internalLoops);
       out.println("Bulges: " + sl.bulge);
    }*/
    
    /**
     * 
     */
    public void checkInternalLoops(){   
        this.internalLoops = 0;
        this.bulge = 0;
        String auxStruct = this.viennaStructure;
        int len, aux, auxI;
        auxStruct = auxStruct.replaceAll("\\.", "a");
        auxStruct = auxStruct.replaceAll("([a-z])\\1+", "$1");
        
        len = auxStruct.length();
        aux = auxStruct.indexOf(")");
        auxI = auxStruct.lastIndexOf("(");
        
        //out.println(auxStruct);
        //out.println(auxStruct.length());
 
        for(int i = 0; (auxI - i) >= 0 && (aux + i) < len; i++ ){
            
            if( auxStruct.charAt(auxI - i) == 'a' && 
                    auxStruct.charAt(aux + i) == 'a' ){
                this.internalLoops++;
            } else if((auxStruct.charAt(auxI - i) == 'a' && 
                    auxStruct.charAt(aux + i) != 'a') ||
                    (auxStruct.charAt(auxI - i) != 'a' && 
                    auxStruct.charAt(aux + i) == 'a')){
                /*out.println("auxI-i: " + (auxI - i) + " : " 
                        + auxStruct.charAt(auxI-i) );
                out.println("aux+i: " +(aux + i)+ ": " 
                        + auxStruct.charAt(aux + i) );*/
                this.bulge++;
            }
        }
        
    }

    /**
     *
     * @param startPosLoop
     */
    public void setNLoop(int startPosLoop) {
        char n2 = ' ', n5 = ' ', n6 = ' ', n7 = ' ', n8 = ' ';
        char precedes = ' ', precedes2 = ' ';
        int matchFirst, matchLast;

        if (reversed) {
            this.rnaHairpinSequence = reverseIt(this.rnaHairpinSequence);
            
            matchFirst = rnaHairpinSequence.indexOf(reverseIt(loop));
            matchLast = rnaHairpinSequence.lastIndexOf(reverseIt(loop));
            
            if(matchFirst == matchLast)
                startPosLoop = matchFirst;
            else
                startPosLoop = matchLast;
            
        }
        
        if (this.rnaHairpinSequence != null) {
            
            n2 = this.rnaHairpinSequence.charAt(startPosLoop + 1);
            precedes = this.rnaHairpinSequence.charAt(startPosLoop - 1);
            precedes2 = this.rnaHairpinSequence.charAt(startPosLoop - 2);
            
            if (this.loop.length() > 4) {
                n5 = this.rnaHairpinSequence.charAt(startPosLoop + 4);

                if (this.loop.length() > 5) {
                    n6 = this.rnaHairpinSequence.charAt(startPosLoop + 5);

                    if (this.loop.length() > 6) {
                        n7 = this.rnaHairpinSequence.charAt(startPosLoop + 6);

                        if (this.loop.length() > 7) {
                            n8 = this.rnaHairpinSequence.charAt(startPosLoop + 7);
                        }

                    }
                }
            }
        }

        this.n2Loop = n2;
        this.n5Loop = n5;
        this.n6Loop = n6;
        this.n7Loop = n7;
        this.n8Loop = n8;
        this.predecessorLoop = precedes;
        this.predecessor2Loop = precedes2;

        if (reversed) {
            this.rnaHairpinSequence = reverseIt(this.rnaHairpinSequence);
        }
    }

    /**
     *
     * @param percA_sequence
     */
    public void setPercA_sequence(float percA_sequence) {
        this.percA_sequence = percA_sequence;
    }

    /**
     *
     * @param percC_sequence
     */
    public void setPercC_sequence(float percC_sequence) {
        this.percC_sequence = percC_sequence;
    }

    /**
     *
     */
    public void setPercent_AG() {

        float myPercent_AG = 0;
        int count;
        count = 0;
        String aux;
        aux = rnaHairpinSequence;

        if (aux != null) {
            count = aux.length() - aux.replace("A", "").length(); //$NON-NLS-1$ //$NON-NLS-2$
            count += aux.length() - aux.replace("G", "").length(); //$NON-NLS-1$ //$NON-NLS-2$

            myPercent_AG = new Float(count) / new Float(aux.length());
        }

        this.percent_AG = myPercent_AG;
    }

    /**
     *
     * @return
     */
    public double getMfe() {
        return mfe;
    }

    /**
     *
     * @param mfe
     */
    public void setMfe(Double mfe) {
        this.mfe = mfe;
    }

    /**
     *
     * @return
     */
    public String getStructure() {
        return viennaStructure;
    }

    /**
     *
     * @param structure
     */
    public void setStructure(String structure) {
        this.viennaStructure = structure;
    }

    /**
     *
     */
    public void setPercent_AU() {
        float percentAU = 0;
        int count = 0;
        int size = 0;
        String aux = rnaHairpinSequence;
        size = (aux.length() - loop.length()) / 2;

        for (int i = 0; i < size; i++) {

            if ((aux.charAt(i) == 'A' && aux.charAt(aux.length() - 1 - i) == 'U')
                    || (aux.charAt(i) == 'U' && aux.charAt(aux.length() - 1 - i) == 'A')) {
                count++;
            }
        }

        percentAU = new Float(count) / new Float(size);

        this.percent_AU = percentAU;
    }

    /**
     *
     */
    public void setPercent_CG() {
        float percentCG;
        percentCG = 0;
        int count = 0;
        int size;
        size = 0;
        String aux = rnaHairpinSequence;
        size = (aux.length() - loop.length()) / 2;

        for (int i = 0; i < size; i++) {

            if ((aux.charAt(i) == 'C' && aux.charAt(aux.length() - 1 - i) == 'G')
                    || (aux.charAt(i) == 'G' && aux.charAt(aux.length() - 1 - i) == 'C')) {
                count++;
            }
        }

        percentCG = new Float(count) / new Float(size);

        this.percent_CG = percentCG;
    }

    /**
     *
     * @param wooble
     */
    public void setPercent_GU(int wooble) {
        float percentGU;
        int size;
        String aux = rnaHairpinSequence;
        size = (aux.length() - this.loop.length()) / 2;

        aux = aux.substring(0, size);

        percentGU = new Float(wooble) / new Float(aux.length());
        this.percent_GU = percentGU;
    }

    /**
     *
     * @param percG_sequence
     */
    public void setPercG_sequence(float percG_sequence) {
        this.percG_sequence = percG_sequence;
    }

    /**
     *
     * @param percU_sequence
     */
    public void setPercU_sequence(float percU_sequence) {
        this.percU_sequence = percU_sequence;
    }

    /**
     *
     * @param pumilioLocations
     */
    public void setAdditionalSeqLocations(List<Integer> pumilioLocations) {
        this.additionalSeqLocations = pumilioLocations;
    }

    /**
     *
     * @param rnaHairpinSequence
     */
    public void setRnaHairpinSequence(String rnaHairpinSequence) {
        this.rnaHairpinSequence = rnaHairpinSequence;
    }

    /**
     *
     * @param sequenceLength
     */
    public void setSequenceLength(int sequenceLength) {
        this.sequenceLength = sequenceLength;
    }

    /**
     *
     * @param startsAt
     */
    public void setStartsAt(int startsAt) {
        this.startsAt = startsAt;
    }

    /**
     *
     * @param id
     */
    public void setTags(String id) {
        switch (this.mode) {
            case ENSEMBL:
                ((EnsemblFastaID) id_fasta).setEnsemblTags(id);
                break;
            case FLYBASE:
                ((FlyBaseFastaID) id_fasta).setFlyBaseTags(id);
                break;
            case BIOMART:
                ((BioMartFastaID) id_fasta).setBioMartTags(id);
                break;
            case GENERIC:
                ((GenericID) id_fasta).setGenericTags(id);
                break;
        }
    }
    
    /**
     * 
     * @param gene
     * @param synonym
     * @param transcript
     * @param description
     * @param cds 
     */
    public void setTags(String gene, String synonym, String transcript, 
            String description,
            String cds) {
        id_fasta.setGeneID(gene);
        id_fasta.setTranscriptID(transcript);
        ((GenBankID) id_fasta).setDescription(description);
        ((GenBankID) id_fasta).setCDS(cds);
        ((GenBankID) id_fasta).setSynonym(synonym);
    }

    /**
     *
     * @return
     */
    public String toRowCSV() {

        return this.id_fasta.toRowCSV()
                + this.getLoopPattern() + SourceFile.ROW_DELIMITER
                + this.getLoopID() + SourceFile.ROW_DELIMITER
                + this.getTerminalPair() + SourceFile.ROW_DELIMITER
                + this.isReverse() + SourceFile.ROW_DELIMITER
                + this.predecessor2Loop + SourceFile.ROW_DELIMITER //n-2
                + this.predecessorLoop + SourceFile.ROW_DELIMITER
                + this.n2Loop + SourceFile.ROW_DELIMITER
                + this.n5Loop + SourceFile.ROW_DELIMITER
                + this.n6Loop + SourceFile.ROW_DELIMITER
                + this.n7Loop + SourceFile.ROW_DELIMITER
                + this.n8Loop + SourceFile.ROW_DELIMITER
                + this.getLoop() + SourceFile.ROW_DELIMITER
                + this.rnaHairpinSequence + SourceFile.ROW_DELIMITER
                + this.additional5Seq + SourceFile.ROW_DELIMITER
                + this.additional3Seq + SourceFile.ROW_DELIMITER
                + this.getHairpinStructure() + SourceFile.ROW_DELIMITER
                + this.getStructure() + SourceFile.ROW_DELIMITER
                + this.getPairments() + SourceFile.ROW_DELIMITER /* para que me de apareamientos */
                + this.getGUPairs() + SourceFile.ROW_DELIMITER
                + this.getMismatches() + SourceFile.ROW_DELIMITER
                + this.internalLoops + SourceFile.ROW_DELIMITER
                + this.sequenceLength + SourceFile.ROW_DELIMITER
                + this.startsAt + SourceFile.ROW_DELIMITER
                + this.endsAt + SourceFile.ROW_DELIMITER
                + getFormattedNumber(this.percA_sequence) + SourceFile.ROW_DELIMITER
                + getFormattedNumber(this.percC_sequence) + SourceFile.ROW_DELIMITER
                + getFormattedNumber(this.percG_sequence) + SourceFile.ROW_DELIMITER
                + getFormattedNumber(this.percU_sequence) + SourceFile.ROW_DELIMITER
                + getFormattedNumber(this.percent_AU) + SourceFile.ROW_DELIMITER
                + getFormattedNumber(this.percent_CG) + SourceFile.ROW_DELIMITER
                + getFormattedNumber(this.percent_GU) + SourceFile.ROW_DELIMITER
                + getFormattedNumber(this.percent_AG) + SourceFile.ROW_DELIMITER
                + getFormattedNumber(this.mfe,5) + SourceFile.ROW_DELIMITER
                + getFormattedNumber(this.relativePos) + SourceFile.ROW_DELIMITER
                + this.getAdditionalSequenceCount() + SourceFile.ROW_DELIMITER
                + this.getAdditionalSequenceLocations() + SourceFile.ROW_DELIMITER;
    }

}
