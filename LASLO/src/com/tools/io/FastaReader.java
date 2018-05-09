package com.tools.io;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.tools.joewandy.Dna;
import com.tools.joewandy.Protein;
import com.tools.joewandy.Sequence;

/**
 * A quick and dirty FASTA reader
 * 
 * @author joewandy
 * 
 */
public class FastaReader {

	// TODO: validate acceptable symbols in file parsed based on content type
	public enum FastaFileContent {
		DNA_SEQUENCE, PROTEIN_SEQUENCE
	};

	private FastaFileContent contentType;
	private String inputFile;
	private List<Sequence> sequenceList;

	public FastaReader(FastaFileContent contentType, String inputFile) {
		this.contentType = contentType;
		this.inputFile = inputFile;
		this.sequenceList = new ArrayList<Sequence>();
	}

	public Sequence getFirstSequence() {
		if (sequenceList.size() < 1) {
			return null;
		}
		return sequenceList.get(0);
	}

	public List<Sequence> processFile() {

		BufferedReader fileReader = null;
		String line = null, aux = "";
		String id = ""; //$NON-NLS-1$
		StringBuffer contentBuffer = new StringBuffer();

		try {

			fileReader = new BufferedReader(new FileReader(inputFile));
			do {

				line = fileReader.readLine();
				if (line != null) {

					line = line.trim();
					if (line.isEmpty()) {
						continue;
					}
					char firstChar = line.charAt(0);
					if (firstChar == '>') {

						// save the previous sequence read
						addToSequenceList(id, contentBuffer);

						// now can get the new id > ..
						id = line.substring(1).trim();

						// start a new content buffer
						contentBuffer = null;
						contentBuffer = new StringBuffer();

					} else if (firstChar == ';') {

						// comment line, skip it

					} else {

						// carry on reading sequence content
						aux = line.trim();
						contentBuffer.append(aux);

					}

				} else {

					// save the final sequence content
					addToSequenceList(id, contentBuffer);

				}

			} while (line != null);

		} catch (FileNotFoundException e) {
			System.out.println("File " + inputFile + " is not found"); //$NON-NLS-1$ //$NON-NLS-2$
			System.exit(1);
		} catch (IOException e) {
			System.out.println("An IO error has occured: " + e.getMessage()); //$NON-NLS-1$
			System.exit(1);
		} finally {
			if (fileReader != null) {
				try {
					fileReader.close();
				} catch (IOException e) {
					// do nothing
				}
			}
		}

		contentBuffer = null;
		return sequenceList;

	}

	private void addToSequenceList(String id, StringBuffer sb) {
		if (sb.length() != 0) {
			Sequence s = null;
			
			try{
			//String content = sb.toString();
			
			if (this.contentType == FastaFileContent.DNA_SEQUENCE) {
				s = new Dna(id, sb.toString());
			} else if (this.contentType == FastaFileContent.PROTEIN_SEQUENCE) {
				s = new Protein(id, sb.toString());
			}
			} catch(Exception e){
				System.err.println("Error en secuencia " + id
						+ ". Longitud: " + sb.length());
			}
			this.sequenceList.add(s);
		}
	}

}
