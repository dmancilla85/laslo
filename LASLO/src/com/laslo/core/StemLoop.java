package com.laslo.core;

import com.tools.io.InputSequence;
import com.tools.io.BioMartFastaID;
import com.tools.io.GenericID;
import com.tools.io.SourceFile;
import com.tools.io.FlyBaseFastaID;
import com.tools.io.EnsemblFastaID;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;
import org.apache.commons.lang.StringUtils;

/**
 * 13/12/2016
 *
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
    protected String structure;
    protected int sequenceLength;
    protected int startsAt;
    protected int endsAt;
    protected int mismatches;
    protected int bulges;
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
    protected double mfe;
    protected List<Integer> additionalSeqLocations;

    public String getAdditional5Seq() {
        return additional5Seq;
    }

    public void setAdditional5Seq(String additional5Seq) {
        this.additional5Seq = additional5Seq;
    }

    public String getAdditional3Seq() {
        return additional3Seq;
    }

    public void setAdditional3Seq(String additional3Seq) {
        this.additional3Seq = additional3Seq;
    }
    
    
    
    private String getFormattedNumber(Double number) {
        NumberFormat numberFormatter
                = NumberFormat.getNumberInstance(Locale.getDefault());
        numberFormatter.setMinimumFractionDigits(3);
        return numberFormatter.format(number);
    }

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
        this.mismatches = 0;
        this.bulges = 0;
        this.relativePos = 0.0;
        this.predecessorLoop = 0;
        this.predecessor2Loop = 0;
        this.n2Loop = this.n5Loop = this.n6Loop
                = this.n7Loop = this.n8Loop = 0;

        this.mfe = (float) 0.0;
        this.structure = "";
        this.additional5Seq = "";
        this.additional3Seq = "";
        this.additionalSeqLocations = new ArrayList<>();
    }

    public int getMismatches() {
        return mismatches;
    }

    public void setMismatches(int mismatches) {
        this.mismatches = mismatches;
    }

    public void setLoopPattern(String pattern) {
        this.loopPattern = pattern;
    }

    public double getRelativePos() {
        return this.relativePos;
    }

    public void setRelativePos(double relativePos) {
        this.relativePos = relativePos;
    }

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
                header = BioMartFastaID.getHeader();
                break;    

            case GENERIC:
                header = GenericID.getHeader();
                break;
        }

        return header
                + "LoopPattern" + SourceFile.ROW_DELIMITER
                + "LoopID" + SourceFile.ROW_DELIMITER
                + "TerminalPair" + SourceFile.ROW_DELIMITER
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
                + "Mismatches" + SourceFile.ROW_DELIMITER
                + "Bulges" + SourceFile.ROW_DELIMITER
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

    public int getEndsAt() {
        return endsAt;
    }

    public String getGeneID() {
        return this.id_fasta.getGeneID();
    }

    public String getGeneSymbol() {
        if (this.id_fasta.getGeneSymbol() != null) {
            return this.id_fasta.getGeneSymbol();
        } else {
            return "";
        }
    }

    public String getGUPairs() {
        long pairs = Math.round(this.percent_GU * this.getPairments());

        return pairs + "";
    }

    public String getHairpinStructure() {
        return hairpinStructure; //drawHairpinStructure();
    }

    public String getLoop() {
        return loop;
    }

    public String getLoopID() {
        /*return this.predecessorLoop + this.loop
                + ":" + this.getGUPairs() + "|" + this.getPairments();*/
        return this.predecessorLoop + this.loop
                + "(" + this.getTerminalPair() + ")|" + this.getPairments();
    }

    public String getLoopPattern() {

        return this.loopPattern;
    }

    public InputSequence getMode() {
        return mode;
    }

    public char getN2Loop() {
        return n2Loop;
    }

    public char getN5Loop() {
        return n5Loop;
    }

    public char getN6Loop() {
        return n6Loop;
    }

    public char getN7Loop() {
        return n7Loop;
    }

    public char getN8Loop() {
        return n8Loop;
    }

    public double getPercA_sequence() {
        return percA_sequence;
    }

    public double getPercC_sequence() {
        return percC_sequence;
    }

    public double getPercent_AG() {
        return percent_AG;
    }

    public double getPercent_AU() {
        return percent_AU;
    }

    public double getPercent_CG() {
        return percent_CG;
    }

    public double getPercent_GU() {
        return percent_GU;
    }

    public double getPercG_sequence() {
        return percG_sequence;
    }

    public double getPercU_sequence() {
        return percU_sequence;
    }

    public char getPredecessorLoop() {
        return predecessorLoop;
    }

    public String getAdditionalSequenceCount() {

        Integer count = 0;

        if (this.additionalSeqLocations != null) {
            count = this.additionalSeqLocations.size();
        }

        return count.toString();
    }

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

    public String getRnaHairpinSequence() {
        return rnaHairpinSequence;
    }

    public int getSequenceLength() {
        return sequenceLength;
    }

    public int getStartsAt() {
        return startsAt;
    }

    public int getPairments() {
        Integer l = (rnaHairpinSequence.length() - loop.length()) / 2;

        return l;
    }

    public String getTerminalPair() {
        Character a = rnaHairpinSequence.charAt(structure.lastIndexOf("("));
        Character b = rnaHairpinSequence.charAt(structure.indexOf(")"));

        return a.toString() + b.toString();
    }

    public String getTranscriptID() {
        return this.id_fasta.getTranscriptID();
    }

    public void setEndsAt(int endsAt) {
        this.endsAt = endsAt;
    }

//	public void setKozakLocations(List<Integer> kozakLocations) {
//		this.kozakLocations = kozakLocations;
//	}
    /*	public void setHairpinStructure(int loopLength, int irLength, List<Integer> mismatchs) {
		this.hairpinStructure = drawHairpinStructure(loopLength, irLength, mismatchs);
	}*/
    public void setLoop(String loopPattern) {
        this.loop = loopPattern;
    }

    public void checkPairments() {

        String seq = this.rnaHairpinSequence;
        int woobleCount = 0, mismatch = 0, bulge = 0;
        int CG = 0, AU = 0;
        StringBuilder aux = new StringBuilder(this.structure);
        int firstIzq = 1;

        if (this.loop.isEmpty() || this.structure.isEmpty()) {
            return;
        }

        // Count mismatchs and bulges
        mismatch = StringUtils.countMatches(this.structure, "(.(")
                + StringUtils.countMatches(this.structure, ").)");

        bulges = StringUtils.countMatches(this.structure, "..(")
                + StringUtils.countMatches(this.structure, ")..");

        try {
            firstIzq = this.structure.lastIndexOf('(');
            int firstDer = this.structure.indexOf(')');

            for (int i = firstIzq; i >= 0 && firstDer < structure.length(); i--) {
                if (structure.charAt(i) == '(') {
                    while (structure.charAt(firstDer) != ')') {
                        firstDer++;
                    }

                    if (structure.charAt(firstDer) == ')') {
                        if (PairmentAnalizer.isComplementaryRNAWooble(seq.charAt(i), seq.charAt(firstDer))) {
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
        this.setMismatches(mismatch);
        this.bulges = bulge;

    }

    public void setNLoop(int startPosLoop) {
        char n2 = ' ', n5 = ' ', n6 = ' ', n7 = ' ', n8 = ' ';

        if (this.rnaHairpinSequence != null) {
            n2 = this.rnaHairpinSequence.charAt(startPosLoop + 1);

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
    }

    public void setPercA_sequence(float percA_sequence) {
        this.percA_sequence = percA_sequence;
    }

    public void setPercC_sequence(float percC_sequence) {
        this.percC_sequence = percC_sequence;
    }

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

    public double getMfe() {
        return mfe;
    }

    public void setMfe(Double mfe) {
        this.mfe = mfe;
    }

    public String getStructure() {
        return structure;
    }

    public void setStructure(String structure) {
        this.structure = structure;
    }

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

    public void setPercent_GU(int wooble) {
        float percentGU;
        int size;
        String aux = rnaHairpinSequence;
        size = (aux.length() - this.loop.length()) / 2;

        aux = aux.substring(0, size);

        percentGU = new Float(wooble) / new Float(aux.length());
        this.percent_GU = percentGU;
    }

    public void setPercG_sequence(float percG_sequence) {
        this.percG_sequence = percG_sequence;
    }

    public void setPercU_sequence(float percU_sequence) {
        this.percU_sequence = percU_sequence;
    }

    /*public void setPolyAdenylationStart(List<Integer> poliALocations, String sequence) {
		
		List<Integer> guSites ,cf1Sites;
		
		// 1. Use 3' utr average
		Iterator<Integer> it = poliALocations.iterator();
		
		while(it.hasNext()){
			
			Integer element = it.next();
			
			if(element < this.sequenceLength * 0.333 - THREE_UTR_AVR){
				poliALocations.set(poliALocations.indexOf(element) , 0) ;
			}
			else{
				// 2. Check the GU - portion
				guSites = Laslo.getPatternLocations(sequence.substring(element), BiologyPatterns.GU_RICH_PATTERN );
				
				Iterator<Integer> it2 = guSites.iterator();
				
				while(it2.hasNext()){
					Integer gu_pos = it2.next();
					
					if( 25 <= gu_pos - element && gu_pos - element <= 45  )
					{
						// 3. Check for the Cleavage Factor 1 signal
						cf1Sites = Laslo.getPatternLocations(sequence.substring(gu_pos), BiologyPatterns.CLEAVAGE_FACTOR1_SITE);
						
						if(cf1Sites.isEmpty())
							guSites.set(guSites.indexOf(gu_pos) , 0) ;
					} 
					else guSites.set(guSites.indexOf(gu_pos) , 0) ;
				}
			
				while(guSites.indexOf(0) > 0 ){
					guSites.remove(guSites.indexOf(0));
				}
				
				if(guSites.isEmpty())
					poliALocations.set(poliALocations.indexOf(element) , 0);
				
			}
		}
		
		// Clean the final list
//		for(int i = 0; i < poliALocations.size(); i++ )
//			if( poliALocations.get(i) != 0)
//				this.poliAdenylationSignalLocations.add(poliALocations.get(i));
	}*/
    public void setPredecessorLoop(int irLength) {

        char precedes = ' ';

        if (this.rnaHairpinSequence != null) {
            precedes = this.rnaHairpinSequence.charAt(irLength - 1);
        }

        this.predecessorLoop = precedes;
    }

    public void setAdditionalSeqLocations(List<Integer> pumilioLocations) {
        this.additionalSeqLocations = pumilioLocations;
    }

    public void setRnaHairpinSequence(String rnaHairpinSequence) {
        this.rnaHairpinSequence = rnaHairpinSequence;
    }

    public void setSequenceLength(int sequenceLength) {
        this.sequenceLength = sequenceLength;
    }

    public void setStartsAt(int startsAt) {
        this.startsAt = startsAt;
    }

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

    public String toRowCSV() {

        return this.id_fasta.toRowCSV()
                + this.getLoopPattern() + SourceFile.ROW_DELIMITER
                + this.getLoopID() + SourceFile.ROW_DELIMITER
                + this.getTerminalPair() + SourceFile.ROW_DELIMITER
                + this.predecessorLoop + SourceFile.ROW_DELIMITER //n-2
                + this.predecessorLoop + SourceFile.ROW_DELIMITER
                + this.n2Loop + SourceFile.ROW_DELIMITER
                + this.n5Loop + SourceFile.ROW_DELIMITER
                + this.n6Loop + SourceFile.ROW_DELIMITER
                + this.n7Loop + SourceFile.ROW_DELIMITER
                + this.n8Loop + SourceFile.ROW_DELIMITER
                + this.loop + SourceFile.ROW_DELIMITER
                + this.rnaHairpinSequence + SourceFile.ROW_DELIMITER
                + this.additional5Seq + SourceFile.ROW_DELIMITER
                + this.additional3Seq + SourceFile.ROW_DELIMITER
                + this.getHairpinStructure() + SourceFile.ROW_DELIMITER
                + this.getStructure() + SourceFile.ROW_DELIMITER
                + this.getPairments() + SourceFile.ROW_DELIMITER /* para que me de apareamientos */
                + this.getGUPairs() + SourceFile.ROW_DELIMITER
                + this.getMismatches() + SourceFile.ROW_DELIMITER
                + this.bulges + SourceFile.ROW_DELIMITER
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
                + getFormattedNumber(this.mfe) + SourceFile.ROW_DELIMITER
                + getFormattedNumber(this.relativePos) + SourceFile.ROW_DELIMITER
                + this.getAdditionalSequenceCount() + SourceFile.ROW_DELIMITER
                + this.getAdditionalSequenceLocations() + SourceFile.ROW_DELIMITER;
    }

}
