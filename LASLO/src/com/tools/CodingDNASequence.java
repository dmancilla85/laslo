/**
 * 
 */
package com.tools;

import java.util.ArrayList;

/**
 * @author David
 *
 */
public class CodingDNASequence{
	public final static String START_CODON_MET = "AUG"; //$NON-NLS-1$
	public final static String START_CODON_LEU = "CUG"; //$NON-NLS-1$
	public final static String STOP_CODON_AMBER = "UAG"; //$NON-NLS-1$
	public final static String STOP_CODON_OCHRE = "UAA"; //$NON-NLS-1$
	public final static String STOP_CODON_OPAL = "UGA"; //$NON-NLS-1$
	
	protected int start;
	protected int end;
	protected String geneID;
	
	CodingDNASequence(){
		this.start = 0;
		this.end = 0;
		this.geneID = ""; //$NON-NLS-1$
	}
	
	CodingDNASequence(int start, int end, String id){
		this.start = start;
		this.end = end;
		this.geneID = id;
	}

	public int getStart() {
		return start;
	}

	public void setStart(int start) {
		this.start = start;
	}

	public int getEnd() {
		return end;
	}

	public void setEnd(int end) {
		this.end = end;
	}
	
	public int length(){
		if ((end - start) > 0)
			return end - start;
		else
			return 0;
	}

	public String getGeneID() {
		return geneID;
	}

	public void setGeneID(String geneID) {
		this.geneID = geneID;
	}

	@Override
	public String toString() {
		return "CodingDNASequence [SIZE= " + length() + //$NON-NLS-1$
				", start=" + start + ", end=" + end + ", geneID=" + geneID + "]"; //$NON-NLS-1$ //$NON-NLS-2$ //$NON-NLS-3$ //$NON-NLS-4$
	}
	
	public static int indexOf(String sequence, String codon, int initialPos){
		int index = -1;
		int i = initialPos;
		boolean found = false;
				
		while(i + 3 < sequence.length() && !found){
			
			if( sequence.substring(i, i + 3).equals(codon) ){
				found = true;
				index = i;
			}
			
			i += 3;
		}
		
		return index;
	}
	
	public static int lastIndexOf(String sequence, String codon, int initialPos){
		int index = -1;
		int i = initialPos;
				
		while(i + 3 < sequence.length()){
			
			if( sequence.substring(i, i + 3).equals(codon) ){
				index = i;
			}
			
			i += 3;
		}
		
		return index;
	}
	
	public static ArrayList<CodingDNASequence> getCDS(String sequence, String id){
		
		ArrayList<CodingDNASequence> cds = new ArrayList<CodingDNASequence>();
		String aux = "", aux2 = ""; //$NON-NLS-1$ //$NON-NLS-2$
		boolean canContinue = sequence.length() > 0;
		int startPos = 0, endPos = 0, startAux = 0, startPosOld;
		
		aux = sequence.toString();
		
		while(canContinue){
			// Find possible start codon for Met
			startPos = aux.indexOf(START_CODON_MET, endPos + 3);
		
			if(startPos > 0){
				// Find possible stop codon - variant 1	
				startPosOld = startPos + 3;
				
				do {
					startAux = aux.indexOf(START_CODON_MET, startPosOld + 3);
					
					if(startAux > 0){
						aux2 = aux.substring(startPos, startAux);
						startPosOld = startAux + 3;
					}
				} while( indexOf(aux2, STOP_CODON_AMBER, 0) < 0 && startAux > 0);
				
				if(startAux > 0){
					aux2 = aux.substring(startPos, startAux);
					
					endPos = startPos + lastIndexOf(aux2, STOP_CODON_AMBER, 0);
					
					// Find possible stop codon - variant 2
					if(endPos <= startPos)
						endPos = startPos + lastIndexOf(aux2, STOP_CODON_OCHRE, 0);
					
					// Find possible stop codon - variant 3
					if(endPos <= startPos)
						endPos = startPos + lastIndexOf(aux2, STOP_CODON_OPAL, 0);
					
				} else {
					endPos = lastIndexOf(aux, STOP_CODON_AMBER, 0);
					
					// Find possible stop codon - variant 2
					if(endPos <= 0)
						endPos = lastIndexOf(aux, STOP_CODON_OCHRE, 0);
					
					// Find possible stop codon - variant 3
					if(endPos <= 0)
						endPos = lastIndexOf(aux, STOP_CODON_OPAL, 0);
				}
				
				if(endPos - startPos> 0 && endPos > 0){
					
					cds.add(new CodingDNASequence(startPos + 1, endPos + 1, id) );
					
				} else canContinue = false;
				
			} else canContinue = false;
			
		}
		
		return cds;
	}
	
}