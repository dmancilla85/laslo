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
  
            String path = "C:\\Users\\David\\Documents\\NetBeansProjects\\lasloProject\\LASLO\\ext";
            String ext = ".fasta";
            String fileName = "fruitfly_chen";
            
            System.out.println("Iniciando...");
			// Generation of the iterator of {id,sequence}
			LinkedHashMap<String, DNASequence> fasta = null;
        try {
            fasta = readFastaDNASequence(new File("./ext/" +fileName + ext), false);
        } catch (IOException ex) {
            Logger.getLogger(UShuffle.class.getName()).log(Level.SEVERE, null, ex);
        }
            makeShuffleSequences(path, fileName, fasta, 5, 2);

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
 * @param filename
 * @param fasta
 * @param nRandoms
 * @param k 
 */
    public static void makeShuffleSequences(String path, String filename,
            LinkedHashMap<String, DNASequence> fasta, int nRandoms, int k) {
        Runtime rt;
        String aux="";
        rt = Runtime.getRuntime();
        String fileNameWithOutExt = filename.replaceFirst("[.][^.]+$", "");
        String destiny = "", sequence = "", header = "";
        boolean mkdirs, delete;
        mkdirs = new File(path + RANDOM_PATH).mkdirs();
		
        for (int i = 1; i <= nRandoms; i++) {
            try {

                destiny = path + RANDOM_PATH + '\\' + fileNameWithOutExt 
                        + "_rnd" + i + ".fa";
		File file = new  File(destiny);
				
                if (file.exists()) {
                    delete = (new File(destiny)).delete();
                }
                
		file.createNewFile();
				
		FileWriter fw = new FileWriter(file.getAbsoluteFile(), true);
		BufferedWriter bw = new BufferedWriter(fw);
		int j = 1;		
                for (Map.Entry<String, DNASequence> entry : fasta.entrySet()){
                        DNASequence element = entry.getValue();
                        header = element.getOriginalHeader();
                        sequence = element.getSequenceAsString();
                        String cmd = COMMAND_SHUFFLE + sequence
                                + " -n 1 -k " + k;
                        Process pr = rt.exec(cmd);

                        try (InputStream in = pr.getInputStream()) {
                                int c;

                                while ((c = in.read()) != -1) {

                                        if (c != '\n') {
                                                aux += (char) c;
                                        } else {
                                                bw.write(">" + header);
                                                bw.newLine();
                                                bw.write(aux);
                                                bw.newLine();
                                                bw.newLine();
                                                aux = "";
                                        }

                                }
                        }

                         int exitVal = pr.waitFor();

                        /*System.out.println(destiny + ": Sequence " + j++
                        + " - ["+ exitVal + "]");*/

                }

                bw.close();
                bw = null;
                fw = null;
				
            } catch (IOException | InterruptedException ex) {
                System.out.println("Error: " + ex.getMessage());
            }
        }

    }
}
