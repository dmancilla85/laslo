package com.tools;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

public class OSValidator {

    private static final String OS = System.getProperty("os.name").toLowerCase();

    public static void main(String[] args) {

        System.out.println(OS);

        if (isWindows()) {
            System.out.println("This is Windows");
        } else if (isMac()) {
            System.out.println("This is Mac");
        } else if (isUnix()) {
            System.out.println("This is Unix or Linux");
        } else if (isSolaris()) {
            System.out.println("This is Solaris");
        } else {
            System.out.println("Your OS is not support!!");
        }
        
        System.out.println("Number of cores: " + getNumberOfCPUCores());
    }

    public static boolean isWindows() {
        return (OS.contains("win"));
    }

    public static boolean isMac() {
        return (OS.contains("mac"));
    }

    public static boolean isUnix() {
        return (OS.contains("nix") || OS.contains("nux") || OS.contains("aix"));
    }

    public static boolean isSolaris() {
        return (OS.contains("sunos"));
    }

    public static String getOS() {
        if (isWindows()) {
            return "win";
        } else if (isMac()) {
            return "osx";
        } else if (isUnix()) {
            return "uni";
        } else if (isSolaris()) {
            return "sol";
        }

        return "err";
    }

    /**
     *
     * @return
     */
    public static int getNumberOfCPUCores() {
        //OSValidator osValidator; 
        //osValidator = new OSValidator();
        String command = "";
        if (OSValidator.isMac()) {
        } else if (OSValidator.isUnix()) {
            command = "lscpu";
        } else if (OSValidator.isWindows()) {
            command = "cmd /C WMIC CPU Get /Format:List";
        } else {
            command = "sysctl -n machdep.cpu.core_count";
        }
        Process process = null;
        int numberOfCores = 0;
        int sockets = 0;
        try {
            if (OSValidator.isMac()) {
                String[] cmd = {"/bin/sh", "-c", command};
                process = Runtime.getRuntime().exec(cmd);
            } else {
                process = Runtime.getRuntime().exec(command);
            }
        } catch (IOException e) {
        }

        BufferedReader reader;
        reader = new BufferedReader(
                new InputStreamReader(process.getInputStream()));
        String line;

        try {
            while ((line = reader.readLine()) != null) {
                if (OSValidator.isMac()) {
                    numberOfCores = line.length() > 0 ? Integer.parseInt(line) : 0;
                } else if (OSValidator.isUnix()) {
                    if (line.contains("Core(s) per socket:")) {
                        numberOfCores = Integer.parseInt(line.split("\\s+")[line.split("\\s+").length - 1]);
                    }
                    if (line.contains("Socket(s):")) {
                        sockets = Integer.parseInt(line.split("\\s+")[line.split("\\s+").length - 1]);
                    }
                } else if (OSValidator.isWindows()) {
                    if (line.contains("NumberOfCores")) {
                        numberOfCores = Integer.parseInt(line.split("=")[1]);
                    }
                }
            }
        } catch (IOException e) {
        }
        if (OSValidator.isUnix()) {
            return numberOfCores * sockets;
        }
        return numberOfCores;
    }

}
