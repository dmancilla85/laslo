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

/**
 *
 * @author David A. Mancilla
 */
public class GenBankID extends SourceFile{
    protected String description;
    protected int cdsStart;
    protected int cdsEnd;
    protected String location;
    
    private final static String HEADER = 
            "AccessionID" + ROW_DELIMITER +
            "Version" + ROW_DELIMITER +
            "Description" + ROW_DELIMITER +
            "CDS_Start" + ROW_DELIMITER +
            "CDS_End" + ROW_DELIMITER + 
            "Location" + ROW_DELIMITER;

    public GenBankID() {
        this.description = "";
        this.cdsStart = 0;
        this.cdsEnd = 0;
    }
    
    public GenBankID(String description, int cdsStart, int cdsEnd) {
        this.description = description;
        this.cdsStart = cdsStart;
        this.cdsEnd = cdsEnd;
    }

    public String getDescription() {
        return description;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public String getLocation() {
        return location;
    }

    public void setLocation(int pos) {
        if(this.cdsEnd != 0){
            if(pos < cdsStart){
                location = "5UTR";
            } else if(pos > cdsEnd){
                location = "3UTR";
            } else location = "CDS";
        }
    }
    
    public int getCdsStart() {
        return cdsStart;
    }

    public int getCdsEnd() {
        return cdsEnd;
    }
     
    public void setCDS(String cds) {
        
        if(cds == null){
            return;
        }
        
        String[] parts = cds.split("\\.\\.");
        
        if(parts.length > 0){
            this.cdsStart = new Integer(parts[0]);
            this.cdsEnd = new Integer(parts[1]);
        }
    }
    
    public static String getHeader() {
        return GenBankID.HEADER;
    }
    
    @Override
    public String toRowCSV() {
        return geneID + ROW_DELIMITER
                + transcriptID + ROW_DELIMITER
                + description + ROW_DELIMITER
                + cdsStart + ROW_DELIMITER
                + cdsEnd + ROW_DELIMITER
                + location + ROW_DELIMITER;
    }
}
