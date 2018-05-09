package com.tools.joewandy;

public class Dna implements Sequence {

	private String sequenceString;
	private String id;
	/**
	 *  A = 00
	 *  C = 01
	 *  G = 10
	 *  U = 11
	 *  T = 11
	 * */
	
	//byte sequence[][];

	public Dna(String id, String dnaString) {
		this(dnaString);
		
		this.id = id;
	}

	public Dna(String sequenceString) {
				
		//this.sequence = new byte[sequenceString.length()][2];
		this.setSequenceString(sequenceString);
		
	}

	public static boolean isComplementaryRNA(byte base1, byte base2){
		
		boolean ret = false;
		
		ret = (base1 == 0x0 && base2 == 0x3) || (base1 == 0x3 && base2 == 0x0) ||
				(base1 == 0x1 && base2 == 0x2) || (base1 == 0x2 && base2 == 0x1);
		
		return ret;
	}
	
	public static boolean isComplementaryRNAWooble(byte base1, byte base2){
		
		boolean ret = false;
		
		ret = (base1 == 0x0 && base2 == 0x3) || (base1 == 0x3 && base2 == 0x0) ||
				(base1 == 0x1 && base2 == 0x2) || (base1 == 0x2 && base2 == 0x1) ||
				(base1 == 0x3 && base2 == 0x2) || (base1 == 0x2 && base2 == 0x3);
		
		return ret;
	}
	
	/*public String getSequenceBits(){
		
		String cadena = "";
		
		System.out.println("Size: " + sequence.length);
		
		for(int i = 0; i < sequence.length; i++){
			for (int j = 0; i < 2; i++) {
		         cadena = cadena.concat(sequence[i][j] == 1 ? "1" : "0");
		      }
		}
		return cadena;
	}*/
	
	@Override
	public String getSequenceString() {
		//return this.sequenceString;
		//String sequenceStr = "";
		
//		for(int i = 0; i < sequence.length; i ++){
//			
//			if(sequence[i][0] == 1){
//				if(sequence[i][1] == 1){
//					sequenceStr = sequenceStr.concat("U");
//				} else
//					sequenceStr = sequenceStr.concat("G");
//			} else {
//				if(sequence[i][1] == 1){
//					sequenceStr = sequenceStr.concat("C");
//				} else
//					sequenceStr = sequenceStr.concat("A");
//			}
//				
//		}
		
		return this.sequenceString;
	}

	@Override
	public void setSequenceString(String sequenceString) {
		
		this.sequenceString = 
				sequenceString.replaceAll("\\s", "").toUpperCase(); //$NON-NLS-1$ //$NON-NLS-2$
		
		/*for(int i = 0; i < sequenceString.length(); i++){
			
			switch(sequenceString.charAt(i)){
				case 'A':
					sequence[i][0] = 0x0;
					sequence[i][1] = 0x0;
					break;
				case 'G':
					sequence[i][0] = 0x1;
					sequence[i][1] = 0x0;
					break;
				case 'C':
					sequence[i][0] = 0x0;
					sequence[i][1] = 0x1;
					break;
				case 'U':
					sequence[i][0] = 0x1;
					sequence[i][1] = 0x1;
					break;
				case 'T':
					sequence[i][0] = 0x1;
					sequence[i][1] = 0x1;
					break;
			}
		}*/
	}

	@Override
	public String getId() {
		return this.id;
	}

	@Override
	public void setId(String id) {
		this.id = id;
	}

	@Override
	public String toString() {
		return "Dna [id=" + this.id + ", sequenceString=" 
	+ this.getSequenceString() + "]"; //$NON-NLS-1$ 
	}

}
