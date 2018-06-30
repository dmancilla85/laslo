package com.tools.fasta;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

// 0. Una vez que funcione esto
// 1. Aplicar file separator en todas las clases -- evitar uso de \\
// 2. Empezar a usar internationalize strings
// 3. Soporte para GenBank
public class FASTACorrector {

    protected final static String HEADER_START = ">";
    protected final static String BLANK = " ";
    protected final static String NEW_LINE = System.getProperty("line.separator");
    protected final static String HEADER_DESC = "Sequence ID";

    public static boolean formatFile(String fileName) {
        try {
            Path filePath = Paths.get(fileName);
            int count = 1;
            List<String> fileContent = new ArrayList<>(Files.readAllLines(filePath, StandardCharsets.UTF_8));
            String line;

            for (int i = 0; i < fileContent.size(); i++) {
                line = fileContent.get(i);

                // put a header if the sequence doesnt have
                if (i == 0 && !line.contains(HEADER_START)) {
                    line = HEADER_START + HEADER_DESC + " " + (count++) + NEW_LINE + line;
                } else if (!line.contains(HEADER_START)) {
                    line = line.replaceAll(BLANK, "");
                    line = line.replaceAll(NEW_LINE, "");
                } else {
                    line = NEW_LINE + line;
                }
                fileContent.set(i, line);
            }

            Files.write(filePath, fileContent, StandardCharsets.UTF_8);
            return true;
        } catch (IOException ex) {
            Logger.getLogger(FASTACorrector.class.getName()).log(Level.SEVERE, null, ex);
            return false;
        }
    }

    public static void main(String[] args) {
        formatFile("C:\\Users\\David\\Downloads\\mitoregulin.fa");
    }
}
