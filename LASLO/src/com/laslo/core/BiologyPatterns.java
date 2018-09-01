/**
 * 
 */
package com.laslo.core;

/**
 * @author David
 *
 */
public class BiologyPatterns {
	// PUMILIO
	public static final String PUM1 = "UGUANAUA"; //$NON-NLS-1$
        public static final String PUM2 = "UGUANAUARNNNNBBBBSCCS";
	// Methyl-6-Adenosine Degenerated Consensus
	public static final String MET6_A_DEG_CONSENSUS = "RRACH";
	// Polyadenilation signal 1
	public final static String POLY_A_SGN_1 = "AUUAAA"; // CPSF Binding site
	// Polyadenilation signal 2
	public final static String POLY_A_SGN_2 = "AUAAAA"; // CPSF Binding site
	
	/** Kozak-like sequences in various eukaryotes
	 * 	Vertebrate						gccRccATGG[2]
		Fruit fly			Arthropoda	cAAacATG[7]
		Yeast 				Ascomycota 	aAaAaAATGTCt[8]
		Slime mold 			Amoebozoa 	aaaAAAATGRna[9]
		Ciliate				Ciliophora	nTaAAAATGRct[9]
		Plasmodium spp.		Apicomplexa	taaAAAATGAan[9]
		Toxoplasma 			Apicomplexa	gncAaaATGg[10]
		Trypanosomatidae	Euglenozoa	nnnAnnATGnC[9]
		Terrestrial plants				AACAATGGC[11]
	*
	*/
	public final static String KOZAK_LIKE_STRONG = "SVVSVVRVVAUGG";
	public final static String KOZAK_LIKE_WEAK = "SVNVNNVVVAUGR"; //
			// "SVNSNNRVVAUGR"; // To match with NM_079620.3 (NCBI)
			// VNNSNVRVVAUGN to match with NM_001297949.1?
	public final static String GU_RICH_PATTERN = "KKKKKK"; // CSTF Binding site
	public final static String CLEAVAGE_FACTOR1_SITE = "UGUAA";
}
