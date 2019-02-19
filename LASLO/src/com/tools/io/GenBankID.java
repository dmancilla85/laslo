/*
 * Copyright (C) 2018 David A. Mancilla
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package com.tools.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.LinkedHashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.io.GenbankReaderHelper;
import org.biojava.nbio.core.sequence.io.GenbankWriterHelper;

/**
 *
 * @author David A. Mancilla
 */
public class GenBankID extends SourceFile {

    protected String description;
    protected int cdsStart;
    protected int cdsEnd;
    protected String location;
    protected String synonym;
    protected String proxyConf;
    private static final String PROXY_FILE = "proxy";

    private final static String HEADER
            = "Gen" + ROW_DELIMITER
            + "GeneSynonym" + ROW_DELIMITER
            + "Note" + ROW_DELIMITER
            + "AccessionID" + ROW_DELIMITER
            + "CDS_Start" + ROW_DELIMITER
            + "CDS_End" + ROW_DELIMITER
            + "Location" + ROW_DELIMITER;

    public GenBankID() {
        this.description = "";
        this.cdsStart = 0;
        this.cdsEnd = 0;
        this.proxyConf = "";
    }

    public GenBankID(String synonym, String description, int cdsStart, int cdsEnd) {
        this.description = description;
        this.cdsStart = cdsStart;
        this.cdsEnd = cdsEnd;
        this.synonym = synonym;
        this.proxyConf = "";
    }

    public static String makeFile(String path,
            LinkedHashMap<String, DNASequence> dnaFile) {
        File file = new File(path + "\\sequence.gb");

        try {
            GenbankWriterHelper.writeNucleotideSequence(file,
                    dnaFile.values());
        } catch (Exception ex) {
            Logger.getLogger(GenBankID.class.getName()).log(Level.SEVERE, null, ex);
        }

        return path + "\\sequence.gb";
    }

    public static String getProxyConfiguration() {

        String proxy = "";

        if(!(new File(PROXY_FILE).exists()))
            return "";
        
        try (BufferedReader br = new BufferedReader(new FileReader(PROXY_FILE))) {

            String sCurrentLine;

            while ((sCurrentLine = br.readLine()) != null) {
                proxy = sCurrentLine;
            }

        } catch (IOException e) {
            //System.out.println(e.getLocalizedMessage());
        }

        return proxy;
    }

    public static LinkedHashMap<String, DNASequence>
            downLoadSequenceForId(String genBankId) throws Exception {
        LinkedHashMap<String, DNASequence> dnaFile = null;
        String request = String
                .format("db=nuccore&id=%s&rettype=gb&retmode=text", genBankId);
        String idFormatted = genBankId;
        URL ncbiGenbank = null;

        if (idFormatted.contains(".")) {
            idFormatted = idFormatted
                    .substring(0, idFormatted.lastIndexOf("."));
        }

        try {
            String proxyConn = getProxyConfiguration();

            if (proxyConn.length() > 0) {
                String proxyParm[] = proxyConn.split(",");

                // defined a proxy connection
                System.setProperty("http.proxyHost", proxyParm[0].trim());
                System.setProperty("http.proxyPort", proxyParm[1].trim());

                // If proxy requires authentication, 
                if (proxyParm.length == 4) {
                    System.setProperty("http.proxyUser", proxyParm[2].trim());
                    System.setProperty("http.proxyPassword", proxyParm[3].trim());
                }
            }

            ncbiGenbank = new URL(
                    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
                    + request);

            dnaFile = GenbankReaderHelper
                    .readGenbankDNASequence(ncbiGenbank.openStream());
        } catch (MalformedURLException ex) {
            System.out.println("ERROR: Malformed URL Exception. Cause: "
                    + ex.getCause().getLocalizedMessage());
            dnaFile = null;
        } catch (IOException ex) {

            if (ex.getLocalizedMessage().contains("400")) {
                System.out.println("ERROR: " + genBankId + " code not found.");
            } else {
                System.out.println("ERROR: IO Exception. Cause: "
                        + ex.getLocalizedMessage());
            }
            dnaFile = null;
        }

        return dnaFile;
    }

    /**
     *
     * @param args
     */
    public static void main(String[] args) {

        try {
            downLoadSequenceForId("NM_001275794.1,NM_005690");
        } catch (Exception ex) {
            Logger.getLogger(GenBankID.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public String getSynonym() {
        return synonym;
    }

    public void setSynonym(String synonym) {
        this.synonym = synonym.replace(';', ',');
    }

    public String getDescription() {
        return description;
    }

    public void setDescription(String description) {
        this.description = description.replace(';', ',');
    }

    public String getLocation() {
        return location;
    }

    public void setLocation(int pos) {
        if (this.cdsEnd != 0) {
            if (pos < cdsStart) {
                location = "5'UTR";
            } else if (pos > cdsEnd) {
                location = "3'UTR";
            } else {
                location = "CDS";
            }
        }
    }

    public int getCdsStart() {
        return cdsStart;
    }

    public int getCdsEnd() {
        return cdsEnd;
    }

    public void setCDS(String cds) {

        if (cds == null) {
            return;
        }

        String[] parts = cds.split("\\.\\.");

        if (parts.length > 0) {
            this.cdsStart = new Integer(parts[0]);
            this.cdsEnd = new Integer(parts[1]);
        }
    }

    public static String getHeader() {
        return GenBankID.HEADER;
    }

    @Override
    public String toRowCSV() {
        return geneID.replace(';', ',') + ROW_DELIMITER
                + synonym + ROW_DELIMITER
                + transcriptID.replace(';', ',') + ROW_DELIMITER
                + description + ROW_DELIMITER
                + cdsStart + ROW_DELIMITER
                + cdsEnd + ROW_DELIMITER
                + location + ROW_DELIMITER;
    }
}
