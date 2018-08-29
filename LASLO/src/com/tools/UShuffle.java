/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.tools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.biojava.nbio.core.sequence.DNASequence;
import static org.biojava.nbio.core.sequence.io.FastaReaderHelper.readFastaDNASequence;

/**
 * fastaSeq.getOriginalHeader()
 * @author David A. Mancilla
 */
public class UShuffle {

    private static final String COMMAND_SHUFFLE = "./ext/ushuffle.exe -s ";
    private static final String RANDOM_PATH = "\\shuffled";
    public static void main(String[] args) {
  
            String path = "C:\\Users\\David\\Desktop\\profile\\";
            String ext = ".fasta";
            String fileName = "repressed_non_bound";
            
			// Generation of the iterator of {id,sequence}
			LinkedHashMap<String, DNASequence> fasta = null;
        try {
            fasta = readFastaDNASequence(new File(fileName), false);
        } catch (IOException ex) {
            Logger.getLogger(UShuffle.class.getName()).log(Level.SEVERE, null, ex);
        }
            makeShuffleSequences(path, fasta, 1, 2);

    }

    /**
     * 
     * @return 
     */
    public static String getRandomDir(){
        return RANDOM_PATH;
    }
    
    /**
     * 
     * @param path
     * @param fasta
     * @param k
     * @param nRandoms 
     */
    public static void makeShuffleSequences(String path, LinkedHashMap<String, 
            DNASequence> fasta, int nRandoms, int k) {
        Runtime rt;
        String aux="";
        rt = Runtime.getRuntime();
        String fileName = ""; // fileName here
        String destiny = "", sequence = "", header = "";
        boolean mkdirs, delete;
        mkdirs = new File(path + RANDOM_PATH).mkdirs();
		
        for (int i = 1; i <= nRandoms; i++) {
            try {

                destiny = path + RANDOM_PATH + '\\' + fileName + "_rnd" + i + ".fa";
				File file = new  File(destiny);
				
                if (file.exists()) 
                    delete = (new File(destiny)).delete();
                
				file.createNewFile();
				
				FileWriter fw = new FileWriter(file.getAbsoluteFile(), true);
				BufferedWriter bw = new BufferedWriter(fw);
				
				for (Map.Entry<String, DNASequence> entry : fasta.entrySet()){
					DNASequence element = entry.getValue();
					header = element.getOriginalHeader();
                                        sequence = element.getSequenceAsString();
					Process pr = rt.exec(COMMAND_SHUFFLE + sequence
						+ "-n " + nRandoms + "-k" + k);

					try (InputStream in = pr.getInputStream()) {
						int c;

						while ((c = in.read()) != -1) {

							if (c != '\n') {
								aux += (char) c;
							} else {
								i++;
								aux = "";
								bw.write(header);
								bw.newLine();
								bw.write(aux);
								bw.newLine();
							}

						}
					}

					 int exitVal = pr.waitFor();
                
					System.out.println("Generating randomized sequences in " + destiny
					+ " - ["+ exitVal + "]");
					
				}
				
				
            } catch (IOException | InterruptedException ex) {
                Logger.getLogger(ShuffleSeq.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

    }
}
