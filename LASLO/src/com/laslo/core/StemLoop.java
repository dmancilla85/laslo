package com.laslo.core;

import com.tools.fasta.InputSequence;
import com.tools.fasta.BioMartFastaID;
import com.tools.fasta.GenericID;
import com.tools.fasta.FastaID;
import com.tools.fasta.FlyBaseFastaID;
import com.tools.fasta.EnsemblFastaID;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * 13/12/2016
 *
 */
public class StemLoop {

    final static String ROW_DELIMITER = ";";

    protected FastaID id_fasta;
    protected String rnaHairpinSequence;
    protected String loop;
    protected String hairpinStructure;
    protected int sequenceLength;
    protected int startsAt;
    protected InputSequence mode;
    protected int endsAt;
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
    protected String loopPattern;
    protected double percU_sequence;
    protected double relativePos;
    protected int    mismatches;
    protected double mfe;
    protected String structure;
    
    // Subsets
    protected List<Integer> pumilioLocations;
    //protected List<Integer> kozakLocations;
    //protected List<Integer> poliAdenylationSignalLocations;
    protected List<Integer> consensusLocations[];

    protected List<Integer> wooble;

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
        this.relativePos = 0.0;
        this.predecessorLoop = 0;
        this.n2Loop = this.n5Loop = this.n6Loop = 
                this.n7Loop = this.n8Loop = 0;

        this.mfe = (float) 0.0;
        this.structure = "";
        
        //this.kozakLocations = new ArrayList<Integer>();
        this.pumilioLocations = new ArrayList<>();
        this.wooble = new ArrayList<>();
        //this.poliAdenylationSignalLocations = new ArrayList<Integer>();
    }

    public int getMismatches() {
        return mismatches;
    }

    public void setMismatches(int mismatches) {
        this.mismatches = mismatches;
    }

    public void setLoopPattern(String pattern){
        this.loopPattern = pattern;
    }

    public double getRelativePos() {
        return this.relativePos;
    }

    public void setRelativePos(double relativePos) {
        this.relativePos = relativePos;
    }
    
    
    
    private String drawHairpinStructure() {
        String theStructure = ""; //$NON-NLS-1$
        String rna = rnaHairpinSequence;
        int pairs = (rna.length() - this.loop.length()) / 2;

        for (int i = 0; i < rna.length(); i++) {
            if (i < pairs) {
                if (PairmentAnalizer.isComplementaryRNA(rna.charAt(i), rna.charAt(rna.length() - 1 - i))) {
                    theStructure = theStructure.concat("(");
                } else if (PairmentAnalizer.isComplementaryRNAWooble(rna.charAt(i), rna.charAt(rna.length() - 1 - i))) {
                    theStructure = theStructure.concat("*");
                } else {
                    theStructure = theStructure.concat(".");
                }
            } else if (i > pairs + this.loop.length() - 1) {
                if (PairmentAnalizer.isComplementaryRNA(rna.charAt(i), rna.charAt(rna.length() - 1 - i))) {
                    theStructure = theStructure.concat(")");
                } else if (PairmentAnalizer.isComplementaryRNAWooble(rna.charAt(i), rna.charAt(rna.length() - 1 - i))) {
                    theStructure = theStructure.concat("*");
                } else {
                    theStructure = theStructure.concat(".");
                }
            } else {
                theStructure = theStructure.concat(".");
            }
        }

        return theStructure;
    }
    
     public String drawHairpinViennaStructure() {
        String theStructure = ""; //$NON-NLS-1$
        String rna = rnaHairpinSequence;
        int pairs = (rna.length() - this.loop.length()) / 2;

        for (int i = 0; i < rna.length(); i++) {
            if (i < pairs) {
                if (PairmentAnalizer.isMismatch(rna.charAt(i), rna.charAt(rna.length() - 1 - i))) {
                    theStructure = theStructure.concat(".");
                } else {
                    theStructure = theStructure.concat("(");
                }
            } else if (i > pairs + this.loop.length() - 1) {
                if (PairmentAnalizer.isMismatch(rna.charAt(i), rna.charAt(rna.length() - 1 - i))) {
                    theStructure = theStructure.concat(".");
                } else {
                    theStructure = theStructure.concat(")");
                }
            } else {
                theStructure = theStructure.concat(".");
            }
        }

        return theStructure;
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

            case GENERIC:
                header = GenericID.getHeader();
                break;
        }

        return header
                + "SeqLength" + FastaID.ROW_DELIMITER //$NON-NLS-1$
                //+ "# Stem-loop" + ROW_DELIMITER //$NON-NLS-1$
                + "StartsAt" + FastaID.ROW_DELIMITER //$NON-NLS-1$
                + "EndsAt" + FastaID.ROW_DELIMITER //$NON-NLS-1$
                + "StemLength" + FastaID.ROW_DELIMITER //$NON-NLS-1$
                + "LoopPattern" + FastaID.ROW_DELIMITER
                + "LoopID" + FastaID.ROW_DELIMITER
                + "GU_Pairs" + FastaID.ROW_DELIMITER
                + "Mismatches" + FastaID.ROW_DELIMITER
                + "Loop" + FastaID.ROW_DELIMITER //$NON-NLS-1$
                + "Terminal pair" + FastaID.ROW_DELIMITER
                + "Seq(RNA)" + FastaID.ROW_DELIMITER //$NON-NLS-1$
                + "Seq(Structure)" + FastaID.ROW_DELIMITER //$NON-NLS-1$
                + "N(-1)" + FastaID.ROW_DELIMITER //$NON-NLS-1$
                + "N(2)" + FastaID.ROW_DELIMITER //$NON-NLS-1$
                + "N(5)" + FastaID.ROW_DELIMITER //$NON-NLS-1$
                + "N(6)" + FastaID.ROW_DELIMITER //$NON-NLS-1$
                + "N(7)" + FastaID.ROW_DELIMITER //$NON-NLS-1$
                + "N(8)" + FastaID.ROW_DELIMITER //$NON-NLS-1$
                + "A_SeqPercent" + FastaID.ROW_DELIMITER //$NON-NLS-1$
                + "G_SeqPercent" + FastaID.ROW_DELIMITER //$NON-NLS-1$
                + "C_SeqPercent" + FastaID.ROW_DELIMITER //$NON-NLS-1$
                + "U_SeqPercent" + FastaID.ROW_DELIMITER //$NON-NLS-1$
                + "AU_percent" + FastaID.ROW_DELIMITER //$NON-NLS-1$
                + "CG_percent" + FastaID.ROW_DELIMITER //$NON-NLS-1$
                + "GU_percent" + FastaID.ROW_DELIMITER //$NON-NLS-1$
                + "Purine_percent" + FastaID.ROW_DELIMITER //$NON-NLS-1$
                + "MFE_InSilico" + FastaID.ROW_DELIMITER
                + "RnaFold_Structure" + FastaID.ROW_DELIMITER
                + "PUMILIOs" + FastaID.ROW_DELIMITER //$NON-NLS-1$
                + "PUMILIO_positions" + FastaID.ROW_DELIMITER //$NON-NLS-1$
                + "Position from start" + FastaID.ROW_DELIMITER;
        /*+ "PolyASignalCount" + ROW_DELIMITER  //$NON-NLS-1$
				+ "PolyASignal_positions" + ROW_DELIMITER //$NON-NLS-1$;*/
        //+ "+3'UTR value" + ROW_DELIMITER
        /*+ "KOZAK" + ROW_DELIMITER
				+ "KOZAK_positions" + ROW_DELIMITER*/
        //+ "+5'UTR value";
    }

    
    
    public String getConsensusCount(int index) {
        Integer count = 0;

        if (this.consensusLocations[index] != null) {
            count = this.consensusLocations[index].size();
        }

        return count.toString();
    }

    public String getConsensusLocations(int index) {

        String locations = ""; //$NON-NLS-1$

        if (this.consensusLocations[index] == null) {
            return "";
        }

        Iterator<Integer> itr = this.consensusLocations[index].iterator();

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
        long pairs = Math.round(this.percent_GU * this.getStemLength());

        return pairs + "";
    }

    public String getHairpinStructure() {
        return drawHairpinStructure();
    }
    

    public String getKozakCount() {

        Integer count = 0;

        /*if(this.kozakLocations != null)
			count = this.kozakLocations.size();*/
        return count.toString();
    }

    public String getKozakLocations() {

        String locations = ""; //$NON-NLS-1$

//		if(this.kozakLocations == null)
//			return "";
//		
//		Iterator<Integer> itr = this.kozakLocations.iterator();
//
//		while (itr.hasNext()) {
//
//			Integer element = itr.next();
//
//			if(itr.hasNext())
//				locations = locations.concat(element + ","); //$NON-NLS-1$
//			else
//				locations = locations.concat(Integer.toString(element));
//		}
        return locations;
    }

    public String getLoop() {
        return loop;
    }

    public String getLoopID() {
        /*return this.predecessorLoop + this.loop
                + ":" + this.getGUPairs() + "|" + this.getStemLength();*/
        return this.predecessorLoop + this.loop
                + "(" + this.getTerminalPair() + ")|" + this.getStemLength();
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

    public String getPoliAdenylationSignalCount() {
        Integer count = 0;

        /*if(this.poliAdenylationSignalLocations != null)
			count = this.poliAdenylationSignalLocations.size();*/
        return count.toString();
    }

    public String getPoliAdenylationSignalLocations() {

        String locations = ""; //$NON-NLS-1$

//		if(this.poliAdenylationSignalLocations == null)
//			return "";
//		
//		Iterator<Integer> itr = this.poliAdenylationSignalLocations.iterator();
//
//		while (itr.hasNext()) {
//
//			Integer element = itr.next();
//
//			if(itr.hasNext())
//				locations = locations.concat(element + ","); //$NON-NLS-1$
//			else
//				locations = locations.concat(Integer.toString(element));
//		}
        return locations;
    }

    public char getPredecessorLoop() {
        return predecessorLoop;
    }

    public String getPumilioCount() {

        Integer count = 0;

        if (this.pumilioLocations != null) {
            count = this.pumilioLocations.size();
        }

        return count.toString();
    }

    public String getPumilioLocations() {

        String locations = ""; //$NON-NLS-1$

        if (this.pumilioLocations == null) {
            return "";
        }

        Iterator<Integer> itr = this.pumilioLocations.iterator();

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

    public int getStemLength() {
        Integer l = (rnaHairpinSequence.length() - loop.length()) / 2;

        return l;
    }
    
    public String getTerminalPair() {
        Character a = rnaHairpinSequence.charAt((rnaHairpinSequence.length() - loop.length()) / 2 - 1);
        Character b = rnaHairpinSequence.charAt((rnaHairpinSequence.length() - loop.length()) / 2 + loop.length());

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

    public void checkPairments(){
        
        int lengthIR = this.rnaHairpinSequence.length() - this.loop.length();
        String seq = rnaHairpinSequence;
        int wooble = 0, mismatch = 0;
        
        if(this.loop.isEmpty() || this.rnaHairpinSequence.isEmpty())
            return;
        
        for(int i = 0; i <(lengthIR)/2; i++){
            if(PairmentAnalizer.isComplementaryRNAWooble(seq.charAt(i), seq.charAt(seq.length() - 1- i)))
                wooble++;
            else
                if(!PairmentAnalizer.isComplementaryRNA(seq.charAt(i), seq.charAt(seq.length() - 1- i)))
                    mismatch++;
        }
        
        this.setMismatches(mismatch);
        this.setPercent_GU(wooble);

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

    public void setPumilioLocations(List<Integer> pumilioLocations) {
        this.pumilioLocations = pumilioLocations;
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

    public void setWoobleLocations(List<Integer> woobleLocations) {
        this.wooble = woobleLocations;
    }

    public String toRowCSV() {

        //long minKozac = 0, maxPoliA = 0;
        /*	for(int i = 0; i < this.kozakLocations.size(); i++){
			if(i == 0)
				minKozac = this.kozakLocations.get(i);
			else
				if(minKozac > this.kozakLocations.get(i) && this.kozakLocations.get(i) !=0 )
					minKozac = this.kozakLocations.get(i);
		}*/
 /*for(int i = 0; i < this.poliAdenylationSignalLocations.size(); i++){
			if(i == 0)
				maxPoliA = this.poliAdenylationSignalLocations.get(i);
			else
				if(maxPoliA < this.poliAdenylationSignalLocations.get(i) 
						&& this.poliAdenylationSignalLocations.get(i) !=0 )
					maxPoliA = this.poliAdenylationSignalLocations.get(i);
		}*/
        return this.id_fasta.toRowCSV() + sequenceLength + FastaID.ROW_DELIMITER
                + startsAt + FastaID.ROW_DELIMITER + endsAt + FastaID.ROW_DELIMITER
                + this.getStemLength() + FastaID.ROW_DELIMITER /* para que me de apareamientos */
                + this.getLoopPattern() + FastaID.ROW_DELIMITER
                + this.getLoopID() + FastaID.ROW_DELIMITER
                + this.getGUPairs() + FastaID.ROW_DELIMITER
                + this.getMismatches() + FastaID.ROW_DELIMITER
                + this.loop + FastaID.ROW_DELIMITER
                + this.getTerminalPair() + FastaID.ROW_DELIMITER
                + this.rnaHairpinSequence + FastaID.ROW_DELIMITER
                + this.getHairpinStructure() + FastaID.ROW_DELIMITER
                + this.predecessorLoop + FastaID.ROW_DELIMITER
                + this.n2Loop + FastaID.ROW_DELIMITER
                + this.n5Loop + FastaID.ROW_DELIMITER
                + this.n6Loop + FastaID.ROW_DELIMITER
                + this.n7Loop + FastaID.ROW_DELIMITER
                + this.n8Loop + FastaID.ROW_DELIMITER
                + Double.toString(this.percA_sequence).replace('.', ',')
                + FastaID.ROW_DELIMITER
                + Double.toString(this.percG_sequence).replace('.', ',')
                + FastaID.ROW_DELIMITER
                + Double.toString(this.percC_sequence).replace('.', ',')
                + FastaID.ROW_DELIMITER
                + Double.toString(this.percU_sequence).replace('.', ',')
                + FastaID.ROW_DELIMITER
                + Double.toString(this.percent_AU).replace('.', ',')
                + FastaID.ROW_DELIMITER
                + Double.toString(this.percent_CG).replace('.', ',')
                + FastaID.ROW_DELIMITER
                + Double.toString(this.percent_GU).replace('.', ',')
                + FastaID.ROW_DELIMITER
                + Double.toString(this.percent_AG).replace('.', ',')
                + FastaID.ROW_DELIMITER
                + Double.toString(this.getMfe()).replace('.', ',')
                + FastaID.ROW_DELIMITER
                + this.getStructure()
                + FastaID.ROW_DELIMITER
                + this.getPumilioCount() + FastaID.ROW_DELIMITER
                + this.getPumilioLocations() + FastaID.ROW_DELIMITER
                + Double.toString(this.getRelativePos()).replace('.', ',') + FastaID.ROW_DELIMITER;
        /*+ this.getPoliAdenylationSignalCount() + FastaID.ROW_DELIMITER
				+ this.getPoliAdenylationSignalLocations() +FastaID.ROW_DELIMITER
				+ maxPoliA + FastaID.ROW_DELIMITER
				+ this.getKozakCount() + FastaID.ROW_DELIMITER
				+ this.getKozakLocations() +FastaID.ROW_DELIMITER
				+ minKozac;*/
    }

}
