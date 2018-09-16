package com.laslo.turner1999;

import com.tools.BaseRelation;
import java.util.ArrayList;

class StackingEnergies{
// STACKING ENERGIES: TERMINAL MISMATCHES AND BASE-PAIRS
static final double 	
		ua_aa = -0.8,
		ua_ac = -0.6,
		ua_ag = -0.8,
		ua_au = -0.6,
		uc_aa = -1.0,
		uc_ac = -0.7,
		uc_ag = -1.0,
		uc_au = -0.8,
		ug_aa = -0.8,
		ug_ac = -0.6,
		ug_ag = -0.8,
		ug_au = -0.6,
		uu_aa = -1.0,
		uu_ac = -0.7,
		uu_ag = -1.0,
		uu_au = -0.8,
		ga_ca = -1.5,
		ga_cc = -1.0,
		ga_cg = -1.4,
		ga_cu = -1.0,
		gc_ca = -1.5,
		gc_cc = -1.1,
		gc_cg = -1.5,
		gc_cu = -1.4,
		gg_ca = -1.4,
		gg_cc = -1.0,
		gg_cg = -1.6,
		gg_cu = -1.0,
		gu_ca = -1.5,
		gu_cc = -0.8,
		gu_cg = -1.5,
		gu_cu = -1.2,
		ca_ga = -1.1,
		ca_gc = -1.1,
		ca_gg = -1.6,
		ca_gu = -1.1,
		cc_ga = -1.5,
		cc_gc = -0.7,
		cc_gg = -1.5,
		cc_gu = -1.0,
		cg_ga = -1.3,
		cg_gc = -1.1,
		cg_gg = -1.4,
		cg_gu = -1.1,
		cu_ga = -1.5,
		cu_gc = -0.5,
		cu_gg = -1.5,
		cu_gu = -0.7,
		ua_ga = -0.3,
		ua_gc = -0.6,
		ua_gg = -0.6,
		ua_gu = -0.6,
		uc_ga = -1.0,
		uc_gc = -0.7,
		uc_gg = -1.0,
		uc_gu = -0.8,
		ug_ga = -0.8,
		ug_gc = -0.6,
		ug_gg = -0.8,
		ug_gu = -0.6,
		uu_ga = -1.0,
		uu_gc = -0.7,
		uu_gg = -1.0,
		uu_gu = -0.6,
		aa_ua = -1.0,
		aa_uc = -0.7,
		aa_ug = -1.1,
		aa_uu = -0.7,
		ac_ua = -0.8,
		ac_uc = -0.6,
		ac_ug = -0.8,
		ac_uu = -0.6,
		ag_ua = -1.1,
		ag_uc = -0.7,
		ag_ug = -1.2,
		ag_uu = -0.7,
		au_ua = -0.8,
		au_uc = -0.5,
		au_ug = -0.8,
		au_uu = -0.5,
		ga_ua = -1.0,
		ga_uc = -0.7,
		ga_ug = -0.5,
		ga_uu = -0.7,
		gc_ua = -0.8,
		gc_uc = -0.6,
		gc_ug = -0.8,
		gc_uu = -0.6,
		gg_ua = -1.1,
		gg_uc = -0.7,
		gg_ug = -0.8,
		gg_uu = -0.7,
		gu_ua = -0.8,
		gu_uc = -0.5,
		gu_ug = -0.8,
		gu_uu = -0.5;
}

class StackingEnthalpies{
	static final double 	ua_aa = -3.9,
			ua_ac = -2.3,
			ua_ag = -3.1,
			ua_au = -2.3,
			uc_aa = 2.0,
			uc_ac = 6.0,
			uc_ag = 2.0,
			uc_au = 4.6,
			ug_aa = -3.5,
			ug_ac = -2.3,
			ug_ag = -3.5,
			ug_au = -2.3,
			uu_aa = 2.0,
			uu_ac = -0.3,
			uu_ag = 2.0,
			uu_au = -1.7,
			ca_ga = -5.2,
			ca_gc = -7.2,
			ca_gg = -7.1,
			ca_gu = -7.2,
			cc_ga = -4.0,
			cc_gc = 0.5,
			cc_gg = -4.0,
			cc_gu = -0.3,
			cg_ga = -5.6,
			cg_gc = -7.2,
			cg_gg = -6.2,
			cg_gu = -7.2,
			cu_ga = -4.0,
			cu_gc = -4.2,
			cu_gg = -4.0,
			cu_gu = -5.0,
			ua_ga = -3.4,
			ua_gc = -2.3,
			ua_gg = -0.6,
			ua_gu = -2.3,
			uc_ga = 2.0,
			uc_gc = 6.0,
			uc_gg = 2.0,
			uc_gu = 4.6,
			ug_ga = -3.5,
			ug_gc = -2.3,
			ug_gg = -3.5,
			ug_gu = -2.3,
			uu_ga = 2.0,
			uu_gc = -0.3,
			uu_gg = 2.0,
			uu_gu = 1.6,
			aa_ua = -4.0,
			aa_uc = -4.3,
			aa_ug = -3.8,
			aa_uu = -4.3,
			ac_ua = -6.3,
			ac_uc = -5.1,
			ac_ug = -6.3,
			ac_uu = -1.4,
			ag_ua = -8.9,
			ag_uc = -4.3,
			ag_ug = -8.9,
			ag_uu = -4.3,
			au_ua = -6.3,
			au_uc = -1.8,
			au_ug = -6.3,
			au_uu = 1.4,
			ga_ua = -4.8,
			ga_uc = -4.3,
			ga_ug = 3.1,
			ga_uu = -4.3,
			gc_ua = -6.3,
			gc_uc = -5.1,
			gc_ug = -6.3,
			gc_uu = -1.4,
			gg_ua = -8.9,
			gg_uc = -4.3,
			gg_ug = -1.5,
			gg_uu = -4.3,
			gu_ua = -6.3,
			gu_uc = -1.8,
			gu_ug = -6.3,
			gu_uu = 1.4;
			
}

class Pair{
	String pair;
	BaseRelation kind;
	boolean isTerminalMismatch;
	int id;
	
	Pair(){
		pair = "";
		kind = BaseRelation.WATSON_CRICK;
		isTerminalMismatch = false;
		id = 1;
	}
	
	Pair(String pair, BaseRelation kind, boolean isTerminal, int id){
		this.pair = pair;
		this.kind = kind;
		this.isTerminalMismatch = isTerminal;
		this.id = id;
	}
	
        @Override
	public String toString(){
		String kind = "";
		String end = ""; 
		
		if(isTerminalMismatch)
			end = " - Terminal";
		
		
		if(this.kind == BaseRelation.WATSON_CRICK)
			kind = "Watson-Crick";
		if(this.kind == BaseRelation.MISMATCH)
			kind = "Mismatch";
		if(this.kind == BaseRelation.WOOBLE)
			kind = "Wooble";
		
		return this.id + " - " + this.pair + "(" + kind + ")" + " " + end;
	}
}

/**
 *
 * @author David
 */
public class FoldingFreeEnergyChange {

    /**
     *
     */
    protected char strand;

    /**
     *
     */
    protected String sequence;

    /**
     *
     */
    protected int	loopStart;

    /**
     *
     */
    protected int 	loopEnd;

    /**
     *
     */
    public ArrayList<Pair> pairs = null;
	
    /**
     *
     * @param strand
     * @param sequence
     * @param loopStart
     * @param loopEnd
     */
    public FoldingFreeEnergyChange(char strand, String sequence, int loopStart, int loopEnd) {
		super();
		this.strand = strand;
		this.sequence = sequence;
		this.loopStart = loopStart;
		this.loopEnd = loopEnd;
		
		makePairs();
	}

    /**
     *
     * @return
     */
    public char getStrand() {
		return strand;
	}

    /**
     *
     * @param strand
     */
    public void setStrand(char strand) {
		
		if(strand != '+' && strand != '-')
			this.strand = '*';
		else
			this.strand = strand;
	}

    /**
     *
     * @return
     */
    public String getSequence() {
		return sequence;
	}

    /**
     *
     * @param sequence
     */
    public void setSequence(String sequence) {
		this.sequence = sequence;
	}

    /**
     *
     * @return
     */
    public int getLoopStart() {
		return loopStart;
	}

    /**
     *
     * @param loopStart
     */
    public void setLoopStart(int loopStart) {
		this.loopStart = loopStart;
	}

    /**
     *
     * @return
     */
    public int getLoopEnd() {
		return loopEnd;
	}

    /**
     *
     * @param loopEnd
     */
    public void setLoopEnd(int loopEnd) {
		this.loopEnd = loopEnd;
	}
	
    /**
     *
     * @return
     */
    public double getTerminalMismatchesEnergy(){
		double value = 0.0;
		
		
		
		return value;
	}
	
    /**
     *
     * @param sequence
     * @return
     */
    public double formula(String sequence){
		/**
		 *  ΔG°37 hairpin (>3 nucleotides in loop) = ΔG°37 initiation (n) + ΔG°37 (terminal mismatch) 
		 *  + ΔG°37 (UU or GA first mismatch) + ΔG°37 (GG first mismatch) + ΔG°37 (special GU closure) 
		 *  + ΔG°37 penalty (all C loops)
		 */
		
		// DESTABILIZING ENERGIES BY SIZE OF LOOP (INTERPOLATE WHERE NEEDED) Size N-1
		double initDestabilizingEnergiesHairpin[] = {0, 0, 5.4, 5.6, 5.7, 5.4, 6.0, 5.5, 6.4, 6.5,
				6.6, 6.7, 6.8, 6.9, 6.9, 7.0, 7.1, 7.1, 7.2, 7.2, 7.3, 7.3, 7.4,
				7.4, 7.5, 7.5, 7.5, 7.6, 7.6, 7.7}; 
		// DESTABILIZING ENTHALPIES BY SIZE OF LOOP (INTERPOLATE WHERE NEEDED)
		@SuppressWarnings("unused")
		double initDestabilizingEnthalpiesHairpin[] =  { 0,0, 1.3 , 4.8, 3.6, 2.9, 1.3, 
				2.9, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 
				5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0};
		
		double deltaG37 = 0.0;
		int loopLength = 0;
		
		if(loopLength > 3){
			deltaG37 += initDestabilizingEnergiesHairpin[loopLength - 1];
		}
		
		
		return deltaG37;
	}
	
    /**
     *
     */
    public void makePairs(){
		
		ArrayList<Pair> pairs = new ArrayList<Pair>();
		int pairsLength = sequence.length();
		String aux = "";
		char base1, base2;
		int end = 0;
		int id = 1;
		boolean isTerminal = false;
		
		if(sequence == null || sequence.length() == 0)
			return;
		
		for(int i = 0; i < loopStart; i++){
			end = pairsLength - i - 1;
			
			base1 = sequence.charAt(i);
			base2 = sequence.charAt(end);
			aux = sequence.substring(i, i + 1) + sequence.substring(end, end + 1);
			
			if(i + 2 == loopStart){
				isTerminal = true;
			}
			
			if( (base1 == 'A' && base2 == 'U') ||
				(base1 == 'U' && base2 == 'A') ||
				(base1 == 'C' && base2 == 'G') ||
				(base1 == 'G' && base2 == 'C'))
				pairs.add(new Pair(aux, BaseRelation.WATSON_CRICK, isTerminal, id++));
			else {
				if((base1 == 'G' && base2 == 'U') ||
				   (base1 == 'U' && base2 == 'G'))
					pairs.add(new Pair(aux, BaseRelation.WOOBLE, isTerminal, id++));
				else
					pairs.add(new Pair(aux, BaseRelation.MISMATCH, isTerminal, id++));
			}
		}
		
		this.pairs = pairs;
	}
	
    /**
     *
     * @param args
     */
    public static void main(String[] args) {
		
		String stemLoop = "CGCIGCUCUGGCGGCAGCG";
		int start = 7;
		int end = 11;
		
		//FoldingFreeEnergyChange test = new FoldingFreeEnergyChange('+', stemLoop, start, end);
		
		/*for(int i = 0; i < test.pairs.size(); i++)
			System.out.println(test.pairs.get(i));*/
                
                // The sequence to be predicted
    byte[] seq = new String("GUGACGUGUGCAAAUGUGACGUGUGCAAAUGUGACGUGUGCAAAU").getBytes();
		
    // Create a new instance of RNAFold4J
    //RNAFoldAPI rfa = new RNAFoldAPI();
    
    // Predict the MFE and corresponding structure
    //MFEData mfe = rfa.getMFE(seq);
    //System.out.println(String.format("Sequence:  %s\nStructure: %s \nMFE:       %f", new String(seq), new String(mfe.structure), mfe.mfe));

    // Predict the base pair probability matrix (equivalent to using the -p option).
    // Returns linearized upper trianguar matrix
    //double[] bppm = rfa.getBppm(seq);
	}
}
