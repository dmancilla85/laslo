/**
 *
 */
package com.tools.fasta;

/**
 * @author David
 *
 */
public class FastaID {

    final static String FIELD_DELIMITER = "|"; //$NON-NLS-1$
    public final static String ROW_DELIMITER = ";"; //$NON-NLS-1$

    protected String transcriptID;
    protected String geneID;
    protected String geneSymbol;
    protected String fieldDelimiter;
    protected String rowDelimiter;
    public static final String LOG_EXT = ".log";
    public static final String CSV_EXT = ".csv";
    public static final String FASTA_EXT = ".fasta";
    public static final String FASTA_EXT_2 = ".fa";

	public static InputSequence detectHeader(String header){
		
		InputSequence out = InputSequence.GENERIC;
		
		if(header.contains("gene:") &&
			header.contains("gene_biotype:") &&
			header.contains("transcript_biotype:") &&
			header.contains("gene_symbol:") &&
			header.contains("description:"))
			out = InputSequence.ENSEMBL;
		else {
			if(header.contains("type=") &&
				header.contains("loc=") &&
				header.contains("name=") &&
				header.contains("dbxref=") &&
				header.contains("score=") &&
				header.contains("MD5=") &&
				header.contains("release="))
				out = InputSequence.FLYBASE;
				
			else {
				if(header.contains("|"))
					out = InputSequence.BIOMART;
			}
		}
			
			
		return out;
	}
	
    public String getTranscriptID() {
        return transcriptID;
    }

    public void setTranscriptID(String transcriptID) {
        this.transcriptID = transcriptID;
    }

    public String getGeneID() {
        return geneID;
    }

    public void setGeneID(String geneID) {
        this.geneID = geneID;
    }

    public String getGeneSymbol() {
        return geneSymbol;
    }

    public void setGeneSymbol(String geneSymbol) {
        this.geneSymbol = geneSymbol;
    } 
    
    public String getRowDelimiter() {
        return rowDelimiter;
    }

    public void setRowDelimiter(String rowDelimiter) {
        this.rowDelimiter = rowDelimiter;
    }

    public String getFieldDelimiter() {
        return fieldDelimiter;
    }

    public void setFieldDelimiter(String fieldDelimiter) {
        this.fieldDelimiter = fieldDelimiter;
    }
    
    public FastaID() {
        this.fieldDelimiter = FIELD_DELIMITER;
        this.rowDelimiter = ROW_DELIMITER;
    }

    public static String getHeader() {
        return "NO HEADER";
    }

    @Override
    public String toString() {
        return "geneID = " + geneID + " - "
                + "transcriptID = " + transcriptID;
    }

    public String toRowCSV() {
        return geneID + ROW_DELIMITER
                + transcriptID + ROW_DELIMITER;
    }
}
