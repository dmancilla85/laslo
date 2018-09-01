/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.tools;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author David A. Mancilla
 */
public class RNAfold {

    private final static String COMMAND_RNAFOLD = "./ext/RNAfold.exe";
    private final static String RNAFOLD_ARGS = " -d2 --noLP --noPS";
    private String structure;
    private double mfe;
    
    /**
     * 
     */
    public RNAfold(){
        this.mfe = 0.0;
        this.structure = "";
    }
    
    /**
     * 
     * @param sequence 
     */
    public RNAfold(String sequence){
        String command = COMMAND_RNAFOLD + RNAFOLD_ARGS;
        String aux = "";
        String output = "";
        int i = 1;

        try {
            Process child = Runtime.getRuntime().exec(command);
            
            OutputStream out = child.getOutputStream();
            out.write(sequence.getBytes());
            out.write(13);
            out.write('@');
            out.write(13);
            out.close();
            
            try (InputStream in = child.getInputStream()) {
                int c;

                while ((c = in.read()) != -1) {

                    if (c != '\n') {
                        aux += (char) c;
                    } else {
                        i++;
                        aux = "";
                    }
                    
                    if(i == 2){
                        output = aux;
                    }
                }
                
                //String result[] = output.split(" ");
                
                //if(result[0].length() > 0)
                
                this.structure = output.trim()
                        .substring(0, sequence.length()) ;
                
                /*if(result.length > 2){
                    for(int j = 2; j < result.length; j++)
                        result[1] += result[j];
                }*/
                
                Pattern regex = Pattern.compile("(\\d+(?:\\.\\d+)?)");
                Matcher matcher = regex.matcher(output);
                if(matcher.find()){
                    this.mfe = new Double(matcher.group(1));
                    this.mfe *= (-1);
                }
                
                /*this.mfe = new Double(output.substring(output.) .
                        replace(")", "").
                        replace("(", "").
                        replace('\r', ' '));*/
            }
            
            out.close();
        } catch (IOException | NumberFormatException ex) {
            Logger.getLogger(RNAfold.class.getName()).log(Level.SEVERE, null, ex);
            System.out.println(Arrays.toString(ex.getStackTrace()));
        }
    }

    /**
     * 
     * @return 
     */
    public String getStructure() {
        return structure;
    }

    /**
     * 
     * @param structure 
     */
    public void setStructure(String structure) {
        this.structure = structure;
    }

    @Override
    /**
     * 
     */
    public String toString() {
        return "RnaFoldInterface{" + "structure=" + structure + ", mfe=" + mfe + '}';
    }

    /**
     * 
     * @return 
     */
    public double getMfe() {
        return mfe;
    }

    /**
     * 
     * @param mfe 
     */
    public void setMfe(float mfe) {
        this.mfe = mfe;
    }
  
    public static void main(String[] args) {

       String sequence = "UAGAGAUCUCUAUGUAUUUCCC";
       RNAfold test;
        test = new RNAfold(sequence);
       System.out.println(test);
    }

}
