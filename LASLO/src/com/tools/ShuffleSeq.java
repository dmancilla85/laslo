/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.tools;

import java.io.File;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author David A. Mancilla
 */
public class ShuffleSeq {

    private static final String COMMAND_SHUFFLESEQ = "./ext/shuffleseq.exe";
    private static final String RANDOM_PATH = "\\shuffled";

    public static void main(String[] args) {
  
            String path = "C:\\Users\\David\\Desktop\\profile\\";
            String ext = ".fasta";
            String fileName = "repressed_non_bound";
            
            makeShuffleSequences(path, fileName, ext, 1000);

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
     * @param fileName
     * @param extension
     * @param nRandoms 
     */
    public static void makeShuffleSequences(String path, String fileName,
            String extension, int nRandoms) {
        Runtime rt;
        rt = Runtime.getRuntime();
        String origin = path + "\\" + fileName + extension;
        String destiny = "";
        boolean mkdirs;
        mkdirs = new File(path + RANDOM_PATH).mkdirs();
        
        for (int i = 1; i <= nRandoms; i++) {
            try {

                destiny = path + RANDOM_PATH + '\\' + fileName + "_rnd" + i + extension;

                if (new File(destiny).exists()) {
                    boolean delete = (new File(destiny)).delete();
                }

                Process pr = rt.exec(COMMAND_SHUFFLESEQ + " -sequence "
                        + origin + " -outseq " + destiny);
                int exitVal = pr.waitFor();
                System.out.println("Generating randomized sequences in " + destiny
                + " - ["+ exitVal + "]");
            } catch (IOException | InterruptedException ex) {
                Logger.getLogger(ShuffleSeq.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

    }
}
