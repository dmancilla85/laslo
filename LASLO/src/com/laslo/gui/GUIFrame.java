/*
 * Copyright (C) 2018 David
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
package com.laslo.gui;

import com.tools.fasta.InputSequence;
import com.laslo.core.LoopCatcher;
import com.tools.ReturnValue;
import java.awt.event.WindowEvent;
import java.io.File;
import java.util.ArrayList;
import java.io.PrintStream;
import static java.lang.System.out;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.ImageIcon;
import javax.swing.JFileChooser;
import javax.swing.JFrame;

/**
 *
 * @author David A. Mancilla
 */
public class GUIFrame extends javax.swing.JFrame {

    /**
     * Creates new form NewJFrame
     */
    public GUIFrame() {
        loopCatcher = new LoopCatcher();
        isRunning = false;
        
        initComponents();
        
        this.jLblError.setText("");
        this.jFTpercMismatch.setValue(25);
        this.jFTpercWooble.setValue(50);
        this.jFTfieldSep.setVisible(false);
        this.jFTnumCols.setVisible(false);
        this.jFTpercMismatch.setVisible(false);
        this.jFTpercWooble.setVisible(false);
        
        PrintStream printStream = new PrintStream(new CustomOutputStream(
                this.jTAConsole));
        System.setOut(printStream);
        System.setErr(printStream);  
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jLabel1 = new javax.swing.JLabel();
        jLabel2 = new javax.swing.JLabel();
        jLabel3 = new javax.swing.JLabel();
        jLabel4 = new javax.swing.JLabel();
        jScrollPane1 = new javax.swing.JScrollPane();
        jTALoopPatterns = new javax.swing.JTextArea();
        jLabel5 = new javax.swing.JLabel();
        jLabel6 = new javax.swing.JLabel();
        jTFPathIn = new javax.swing.JTextField();
        jTFPathOut = new javax.swing.JTextField();
        jButtonIn = new javax.swing.JButton();
        jButtonOut = new javax.swing.JButton();
        jBStart = new javax.swing.JButton();
        jSpinMinLength = new javax.swing.JSpinner();
        jSpinMaxLength = new javax.swing.JSpinner();
        jSpinWooble = new javax.swing.JSpinner();
        jSpinMismatch = new javax.swing.JSpinner();
        jLabel7 = new javax.swing.JLabel();
        jCBOrigin = new javax.swing.JComboBox<>();
        jLblError = new javax.swing.JLabel();
        jScrollPane2 = new javax.swing.JScrollPane();
        jTAConsole = new javax.swing.JTextArea();
        jFTnumCols = new javax.swing.JFormattedTextField();
        jFTfieldSep = new javax.swing.JFormattedTextField();
        jFTpercMismatch = new javax.swing.JFormattedTextField();
        jFTpercWooble = new javax.swing.JFormattedTextField();
        jcbExtended = new javax.swing.JCheckBox();
        jMenuBar1 = new javax.swing.JMenuBar();
        jMenuFile = new javax.swing.JMenu();
        jMIExit = new javax.swing.JMenuItem();
        jMenuHelp = new javax.swing.JMenu();
        jMIAbout = new javax.swing.JMenuItem();

        setDefaultCloseOperation(javax.swing.WindowConstants.DO_NOTHING_ON_CLOSE);
        setTitle("LASLO v0.9");
        setFont(new java.awt.Font("Calibri", 0, 12)); // NOI18N
        setIconImage(new ImageIcon(getClass().getResource("/resources/noun_655767_cc.png")).getImage());
        setLocationByPlatform(true);
        setMinimumSize(new java.awt.Dimension(533, 437));
        setResizable(false);
        addWindowListener(new java.awt.event.WindowAdapter() {
            public void windowClosing(java.awt.event.WindowEvent evt) {
                formWindowClosing(evt);
            }
        });

        jLabel1.setFont(new java.awt.Font("Calibri", 0, 12)); // NOI18N
        jLabel1.setText("Longitud de los stem entre: ");

        jLabel2.setFont(new java.awt.Font("Calibri", 0, 12)); // NOI18N
        jLabel2.setText("Pares wooble admitidos: ");

        jLabel3.setFont(new java.awt.Font("Calibri", 0, 12)); // NOI18N
        jLabel3.setText("LOOPS");

        jLabel4.setFont(new java.awt.Font("Calibri", 0, 12)); // NOI18N
        jLabel4.setText("UBICACION");

        jTALoopPatterns.setColumns(20);
        jTALoopPatterns.setFont(new java.awt.Font("Calibri Light", 0, 12)); // NOI18N
        jTALoopPatterns.setRows(3);
        jTALoopPatterns.setTabSize(2);
        jTALoopPatterns.setText("CNGG, CNGGN, CNGGNN, CNGGNNN, CNGGNNNN");
        jTALoopPatterns.setToolTipText("Ingrese los patrones de loop a buscar separados por coma");
        jTALoopPatterns.setBorder(javax.swing.BorderFactory.createEmptyBorder(1, 1, 1, 1));
        jScrollPane1.setViewportView(jTALoopPatterns);

        jLabel5.setFont(new java.awt.Font("Calibri", 0, 12)); // NOI18N
        jLabel5.setText("Mismatchs admitidos: ");

        jLabel6.setFont(new java.awt.Font("Calibri", 0, 12)); // NOI18N
        jLabel6.setText("DESTINO");

        jTFPathIn.setFont(new java.awt.Font("Calibri", 0, 12)); // NOI18N
        jTFPathIn.setToolTipText("Ubicaciòn de los archivos FASTA a revisar");
        jTFPathIn.setBorder(javax.swing.BorderFactory.createEtchedBorder());
        jTFPathIn.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jTFPathInActionPerformed(evt);
            }
        });

        jTFPathOut.setFont(new java.awt.Font("Calibri", 0, 12)); // NOI18N
        jTFPathOut.setToolTipText("Donde se guardarán los resultados del proceso");
        jTFPathOut.setBorder(javax.swing.BorderFactory.createEtchedBorder());
        jTFPathOut.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jTFPathOutActionPerformed(evt);
            }
        });

        jButtonIn.setIcon(new javax.swing.ImageIcon(getClass().getResource("/resources/noun_53223_cc.png"))); // NOI18N
        jButtonIn.setBorder(new javax.swing.border.SoftBevelBorder(javax.swing.border.BevelBorder.RAISED));
        jButtonIn.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jButtonInActionPerformed(evt);
            }
        });

        jButtonOut.setIcon(new javax.swing.ImageIcon(getClass().getResource("/resources/noun_53223_cc.png"))); // NOI18N
        jButtonOut.setBorder(new javax.swing.border.SoftBevelBorder(javax.swing.border.BevelBorder.RAISED));
        jButtonOut.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jButtonOutActionPerformed(evt);
            }
        });

        jBStart.setText("Comenzar");
        jBStart.setBorder(new javax.swing.border.SoftBevelBorder(javax.swing.border.BevelBorder.RAISED));
        jBStart.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jBStartActionPerformed(evt);
            }
        });

        jSpinMinLength.setPreferredSize(new java.awt.Dimension(35, 22));
        jSpinMinLength.setValue(4);
        jSpinMinLength.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                jSpinMinLengthStateChanged(evt);
            }
        });

        jSpinMaxLength.setMinimumSize(new java.awt.Dimension(35, 30));
        jSpinMaxLength.setPreferredSize(new java.awt.Dimension(45, 22));
        jSpinMaxLength.setValue(12);
        jSpinMaxLength.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                jSpinMaxLengthStateChanged(evt);
            }
        });

        jSpinWooble.setPreferredSize(new java.awt.Dimension(35, 22));
        jSpinWooble.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                jSpinWoobleStateChanged(evt);
            }
        });

        jSpinMismatch.setPreferredSize(new java.awt.Dimension(35, 22));
        jSpinMismatch.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                jSpinMismatchStateChanged(evt);
            }
        });

        jLabel7.setFont(new java.awt.Font("Calibri", 0, 12)); // NOI18N
        jLabel7.setText("ORIGEN");

        jCBOrigin.setFont(new java.awt.Font("Calibri", 0, 12)); // NOI18N
        jCBOrigin.setModel(new javax.swing.DefaultComboBoxModel<>(new String[] { "Ensembl", "FlyBase", "BioMart", "Otro" }));
        jCBOrigin.setToolTipText("Origen del archivo FASTA");
        jCBOrigin.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                jCBOriginItemStateChanged(evt);
            }
        });
        jCBOrigin.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jCBOriginActionPerformed(evt);
            }
        });

        jLblError.setFont(new java.awt.Font("Calibri", 0, 12)); // NOI18N
        jLblError.setForeground(new java.awt.Color(255, 0, 0));
        jLblError.setText("ERROR");

        jTAConsole.setEditable(false);
        jTAConsole.setBackground(new java.awt.Color(204, 204, 204));
        jTAConsole.setColumns(20);
        jTAConsole.setFont(new java.awt.Font("Consolas", 0, 13)); // NOI18N
        jTAConsole.setLineWrap(true);
        jTAConsole.setRows(5);
        jTAConsole.setBorder(javax.swing.BorderFactory.createTitledBorder("Consola"));
        jScrollPane2.setViewportView(jTAConsole);

        try {
            jFTnumCols.setFormatterFactory(new javax.swing.text.DefaultFormatterFactory(new javax.swing.text.MaskFormatter("#")));
        } catch (java.text.ParseException ex) {
            ex.printStackTrace();
        }
        jFTnumCols.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        jFTnumCols.setToolTipText("Separador de campos");

        try {
            jFTfieldSep.setFormatterFactory(new javax.swing.text.DefaultFormatterFactory(new javax.swing.text.MaskFormatter("*")));
        } catch (java.text.ParseException ex) {
            ex.printStackTrace();
        }
        jFTfieldSep.setHorizontalAlignment(javax.swing.JTextField.CENTER);
        jFTfieldSep.setToolTipText("Separador de campos");

        try {
            jFTpercMismatch.setFormatterFactory(new javax.swing.text.DefaultFormatterFactory(new javax.swing.text.MaskFormatter("##%")));
        } catch (java.text.ParseException ex) {
            ex.printStackTrace();
        }
        jFTpercMismatch.setHorizontalAlignment(javax.swing.JTextField.RIGHT);
        jFTpercMismatch.setToolTipText("Separador de campos");

        try {
            jFTpercWooble.setFormatterFactory(new javax.swing.text.DefaultFormatterFactory(new javax.swing.text.MaskFormatter("##%")));
        } catch (java.text.ParseException ex) {
            ex.printStackTrace();
        }
        jFTpercWooble.setHorizontalAlignment(javax.swing.JTextField.RIGHT);
        jFTpercWooble.setText("%");
        jFTpercWooble.setToolTipText("Separador de campos");

        jcbExtended.setText("Modo extendido");
        jcbExtended.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jcbExtendedActionPerformed(evt);
            }
        });

        jMenuFile.setText("Archivo");

        jMIExit.setText("Salir");
        jMIExit.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMIExitActionPerformed(evt);
            }
        });
        jMenuFile.add(jMIExit);

        jMenuBar1.add(jMenuFile);

        jMenuHelp.setText("Ayuda");

        jMIAbout.setText("Acerca de...");
        jMenuHelp.add(jMIAbout);

        jMenuBar1.add(jMenuHelp);

        setJMenuBar(jMenuBar1);

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jScrollPane1)
                    .addComponent(jScrollPane2)
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jLabel4)
                            .addComponent(jLabel6)
                            .addComponent(jLabel7))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(jCBOrigin, javax.swing.GroupLayout.PREFERRED_SIZE, 149, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(jFTfieldSep, javax.swing.GroupLayout.PREFERRED_SIZE, 25, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(jFTnumCols, javax.swing.GroupLayout.PREFERRED_SIZE, 25, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addGap(0, 203, Short.MAX_VALUE))
                            .addComponent(jTFPathIn)
                            .addComponent(jTFPathOut))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jButtonIn, javax.swing.GroupLayout.Alignment.TRAILING)
                            .addComponent(jButtonOut, javax.swing.GroupLayout.Alignment.TRAILING)))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(jBStart, javax.swing.GroupLayout.PREFERRED_SIZE, 81, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jLblError, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jcbExtended)
                            .addGroup(layout.createSequentialGroup()
                                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                                    .addComponent(jLabel1, javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(jLabel2, javax.swing.GroupLayout.PREFERRED_SIZE, 152, javax.swing.GroupLayout.PREFERRED_SIZE)
                                    .addComponent(jLabel5, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.PREFERRED_SIZE, 152, javax.swing.GroupLayout.PREFERRED_SIZE))
                                .addGap(10, 10, 10)
                                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                                    .addComponent(jSpinMismatch, javax.swing.GroupLayout.PREFERRED_SIZE, 45, javax.swing.GroupLayout.PREFERRED_SIZE)
                                    .addComponent(jSpinWooble, javax.swing.GroupLayout.PREFERRED_SIZE, 45, javax.swing.GroupLayout.PREFERRED_SIZE)
                                    .addComponent(jSpinMinLength, javax.swing.GroupLayout.PREFERRED_SIZE, 45, javax.swing.GroupLayout.PREFERRED_SIZE))
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(jSpinMaxLength, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                    .addComponent(jFTpercWooble, javax.swing.GroupLayout.PREFERRED_SIZE, 33, javax.swing.GroupLayout.PREFERRED_SIZE)
                                    .addComponent(jFTpercMismatch, javax.swing.GroupLayout.PREFERRED_SIZE, 33, javax.swing.GroupLayout.PREFERRED_SIZE)))
                            .addComponent(jLabel3))
                        .addGap(0, 0, Short.MAX_VALUE)))
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                        .addComponent(jLabel4)
                        .addComponent(jTFPathIn, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addComponent(jButtonIn, javax.swing.GroupLayout.PREFERRED_SIZE, 20, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                        .addComponent(jLabel6)
                        .addComponent(jTFPathOut))
                    .addComponent(jButtonOut, javax.swing.GroupLayout.PREFERRED_SIZE, 20, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel7)
                    .addComponent(jCBOrigin, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jFTfieldSep, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jFTnumCols, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(18, 18, 18)
                .addComponent(jcbExtended)
                .addGap(57, 57, 57)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel1)
                    .addComponent(jSpinMinLength, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jSpinMaxLength, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jSpinWooble, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel2)
                    .addComponent(jFTpercWooble, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jSpinMismatch, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel5)
                    .addComponent(jFTpercMismatch, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jLabel3)
                .addGap(18, 18, 18)
                .addComponent(jScrollPane1, javax.swing.GroupLayout.PREFERRED_SIZE, 50, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jBStart)
                    .addComponent(jLblError))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jScrollPane2, javax.swing.GroupLayout.PREFERRED_SIZE, 153, javax.swing.GroupLayout.PREFERRED_SIZE))
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void jTFPathInActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jTFPathInActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_jTFPathInActionPerformed

    private void jTFPathOutActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jTFPathOutActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_jTFPathOutActionPerformed

    private void jButtonOutActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jButtonOutActionPerformed
        // TODO add your handling code here:
        final JFileChooser fc = new JFileChooser();
        JFrame myFrame = new JFrame();
        File folder = null;

        // In response to a button click:
        fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        int returnVal = fc.showDialog(myFrame, "Seleccione carpeta de salida");

        if (JFileChooser.APPROVE_OPTION == returnVal) {
            folder = fc.getSelectedFile();
        }

        if (folder != null) {
            this.jTFPathOut.setText(folder.getAbsolutePath());
        }
    }//GEN-LAST:event_jButtonOutActionPerformed

    private void jMIExitActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMIExitActionPerformed
        // TODO add your handling code here:
        this.dispatchEvent(new WindowEvent(this, WindowEvent.WINDOW_CLOSING));
    }//GEN-LAST:event_jMIExitActionPerformed

    private void jBStartActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jBStartActionPerformed
        // TODO add your handling code here:
        // Llamar a método unificador
        this.start();
    }//GEN-LAST:event_jBStartActionPerformed

    private void jSpinMinLengthStateChanged(javax.swing.event.ChangeEvent evt) {//GEN-FIRST:event_jSpinMinLengthStateChanged
        // TODO add your handling code here:
        if (new Integer(this.jSpinMinLength.getValue().toString()) < 3) {
            this.jSpinMinLength.setValue(3);
        }
    }//GEN-LAST:event_jSpinMinLengthStateChanged

    private void jSpinWoobleStateChanged(javax.swing.event.ChangeEvent evt) {//GEN-FIRST:event_jSpinWoobleStateChanged
        // TODO add your handling code here:
        // TODO add your handling code here:
        int aux = new Integer(this.jSpinWooble.getValue().toString());
        int max = new Integer(this.jSpinMaxLength.getValue().toString());

        if (aux < 0) {
            this.jSpinWooble.setValue(0);
        } else if (aux > max) {
            this.jSpinWooble.setValue(max);
        }
    }//GEN-LAST:event_jSpinWoobleStateChanged

    private void jSpinMismatchStateChanged(javax.swing.event.ChangeEvent evt) {//GEN-FIRST:event_jSpinMismatchStateChanged
        // TODO add your handling code here:
        int aux = new Integer(this.jSpinMismatch.getValue().toString());
        int max = new Integer(this.jSpinMaxLength.getValue().toString());

        if (aux < 0) {
            this.jSpinMismatch.setValue(0);
        } else if (aux > max) {
            this.jSpinMismatch.setValue(max);
        }
    }//GEN-LAST:event_jSpinMismatchStateChanged

    private void jSpinMaxLengthStateChanged(javax.swing.event.ChangeEvent evt) {//GEN-FIRST:event_jSpinMaxLengthStateChanged
        // TODO add your handling code here:
        int aux = new Integer(this.jSpinMaxLength.getValue().toString());
        int min = new Integer(this.jSpinMinLength.getValue().toString());

        if (aux < min) {
            this.jSpinMaxLength.setValue(min);
        }
    }//GEN-LAST:event_jSpinMaxLengthStateChanged

    private void jCBOriginItemStateChanged(java.awt.event.ItemEvent evt) {//GEN-FIRST:event_jCBOriginItemStateChanged
        // TODO add your handling code here:
    }//GEN-LAST:event_jCBOriginItemStateChanged

    private void jButtonInActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jButtonInActionPerformed
        // TODO add your handling code here:
        final JFileChooser fc = new JFileChooser();
        JFrame myFrame = new JFrame();
        File file = null;
        boolean isFolder = true;
        this.listOfFiles = null;

        // In response to a button click:
        fc.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
        int returnVal = fc.showDialog(myFrame, "Seleccione carpeta o archivo");

        if (returnVal == JFileChooser.APPROVE_OPTION) {
            file = fc.getSelectedFile();
            isFolder = file.isDirectory();
        } else
            return;

        if (file.exists()) {
            if(isFolder){
               this.listOfFiles = file.listFiles();
               this.jTFPathOut.setText(file.getAbsolutePath()); 
            } else {
                this.listOfFiles = new File[1];
                this.listOfFiles[0] = file;
                this.jTFPathOut.setText(file.getParent() );
            }
            this.jTFPathIn.setText(file.getAbsolutePath());
        }
    }//GEN-LAST:event_jButtonInActionPerformed

    private void formWindowClosing(java.awt.event.WindowEvent evt) {//GEN-FIRST:event_formWindowClosing
        // TODO add your handling code here:
        if(this.isRunning){
            MessageBox.show("No se puede salir de la aplicación mientras el " +
                    "proceso se está ejecutando", "Proceso en ejecución");
   
        } else
            System.exit(ReturnValue.SUCCESS.getReturnCode());
    }//GEN-LAST:event_formWindowClosing

    private void jCBOriginActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jCBOriginActionPerformed
        // TODO add your handling code here:
        int index = this.jCBOrigin.getSelectedIndex();
        setInputSequence(this.jCBOrigin.getItemAt(index));
    }//GEN-LAST:event_jCBOriginActionPerformed

    private void jcbExtendedActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jcbExtendedActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_jcbExtendedActionPerformed

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        /* Set the Nimbus look and feel */
        //<editor-fold defaultstate="collapsed" desc=" Look and feel setting code (optional) ">
        /* If Nimbus (introduced in Java SE 6) is not available, stay with the default look and feel.
         * For details see http://download.oracle.com/javase/tutorial/uiswing/lookandfeel/plaf.html 
         */
        try {
            for (javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels()) {
                if ("Nimbus".equals(info.getName())) {
                    javax.swing.UIManager.setLookAndFeel(info.getClassName());
                    break;
                }
            }
        } catch (ClassNotFoundException | InstantiationException | IllegalAccessException | javax.swing.UnsupportedLookAndFeelException ex) {
            java.util.logging.Logger.getLogger(GUIFrame.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>
        //</editor-fold>

        //</editor-fold>
        //</editor-fold>

        /* Create and display the form */
        java.awt.EventQueue.invokeLater(() -> {
            new GUIFrame().setVisible(true);
        });
    }

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton jBStart;
    private javax.swing.JButton jButtonIn;
    private javax.swing.JButton jButtonOut;
    private javax.swing.JComboBox<String> jCBOrigin;
    private javax.swing.JFormattedTextField jFTfieldSep;
    private javax.swing.JFormattedTextField jFTnumCols;
    private javax.swing.JFormattedTextField jFTpercMismatch;
    private javax.swing.JFormattedTextField jFTpercWooble;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel jLabel3;
    private javax.swing.JLabel jLabel4;
    private javax.swing.JLabel jLabel5;
    private javax.swing.JLabel jLabel6;
    private javax.swing.JLabel jLabel7;
    private javax.swing.JLabel jLblError;
    private javax.swing.JMenuItem jMIAbout;
    private javax.swing.JMenuItem jMIExit;
    private javax.swing.JMenuBar jMenuBar1;
    private javax.swing.JMenu jMenuFile;
    private javax.swing.JMenu jMenuHelp;
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JScrollPane jScrollPane2;
    private javax.swing.JSpinner jSpinMaxLength;
    private javax.swing.JSpinner jSpinMinLength;
    private javax.swing.JSpinner jSpinMismatch;
    private javax.swing.JSpinner jSpinWooble;
    private javax.swing.JTextArea jTAConsole;
    private javax.swing.JTextArea jTALoopPatterns;
    private javax.swing.JTextField jTFPathIn;
    private javax.swing.JTextField jTFPathOut;
    private javax.swing.JCheckBox jcbExtended;
    // End of variables declaration//GEN-END:variables

    protected LoopCatcher loopCatcher;
    protected File[] listOfFiles;
    protected boolean isRunning;

    public void setInputSequence(String value){
     InputSequence origin;
     String fs = "";
     String cols = "0";
     
        switch (value) {
            case "Ensembl":
                origin = InputSequence.ENSEMBL;
                fs = ";";
                cols = "10";
                break;
            case "FlyBase":
                origin = InputSequence.FLYBASE;
                fs = ",";
                cols = "8";
                break;
            case "BioMart":
                origin = InputSequence.BIOMART;
                fs = "|";
                cols = "6";
                break;
            default:
                origin = InputSequence.GENERIC;
                fs = this.jFTfieldSep.getText();
                cols = this.jFTnumCols.getText();
        }

        this.jFTfieldSep.setValue(fs);
        this.jFTnumCols.setValue(cols);
        this.loopCatcher.setInputType(origin);
    }
    
    public void setIsRunning(boolean value){
        this.isRunning = value;
        this.jBStart.setEnabled(!value);
        
        if(value){
            out.println("Comenzando proceso...");
            this.jBStart.setText("Espere...");
        } else {
            this.jBStart.setText("Comenzar");
        }
    }
    
    public boolean validateParameters(int min, int max, int wooble,
            int mismatch, String pathIn, String pathOut,
            ArrayList<String> loopList) {
        boolean isValid;
        isValid = true;
        String aux;

        // Validate numbers
        if (min > max || wooble > max || mismatch > max
                || min < 3 || wooble < 0 || mismatch < 0) {
            this.jLblError.setText("Error! Revise los parámetros.");
            isValid = false;
        }

        // Validate paths
        if (!new File(pathOut).exists()) {
            this.jLblError.setText("Error! Especificar carpeta de destino.");
            isValid = false;
        }
        if (!new File(pathIn).exists()) {
            this.jLblError.setText("Error! Especificar fuente.");
            isValid = false;
        }

        if (loopList.size() <= 0 || loopList.get(0).trim().length() == 0) {
            this.jLblError.setText("Error! Ingrese motivos a buscar.");
            isValid = false;
        }
        
        for(int i = 0; i < loopList.size() && isValid; i++){
            aux = loopList.get(i);
            aux = aux.trim();
            if(aux.length() <= 0){
                isValid = false;
                this.jLblError.setText("Error! Revise motivos a buscar.");
            }
        }

        return isValid;
    }
    
    public void start() {

        String inputValue;
        inputValue = this.jCBOrigin.getItemAt(jCBOrigin.getSelectedIndex());
        int max, min, wooble, mismatch;
        String pathOut = this.jTFPathOut.getText();
        String pathIn = this.jTFPathIn.getText();
        String[] loops = this.jTALoopPatterns.getText().split(",");
        ArrayList<String> loopList;
        loopList = new ArrayList<>();
        this.jLblError.setText("");

        if(this.isRunning)
            return;
        
        this.setIsRunning(true);
        this.setInputSequence(inputValue);
            
        // Values
        min = new Integer(this.jSpinMinLength.getValue().toString());
        max = new Integer(this.jSpinMaxLength.getValue().toString());
        wooble = new Integer(this.jSpinWooble.getValue().toString());
        mismatch = new Integer(this.jSpinMismatch.getValue().toString());

        // Loops
        loopList.addAll(Arrays.asList(loops));

        // 
        if (!validateParameters(min, max, wooble, mismatch, pathIn, pathOut,
                loopList)) {
            setIsRunning(false);
            return;
        }

        // Start process
        out.flush();
        loopCatcher.setLoopPatterns(loopList);
        loopCatcher.setMaxLength(max);
        loopCatcher.setMinLength(min);
        loopCatcher.setMaxMismatch(mismatch);
        loopCatcher.setMaxWooble(wooble);
        loopCatcher.setPathOut(pathOut);
        loopCatcher.setFileList(this.listOfFiles);
        loopCatcher.setIsExtendedMode(this.jcbExtended.isSelected() );
        
        Thread thread = new Thread((Runnable)loopCatcher);
        thread.start();
        
        try {
            thread.join();
        } catch (InterruptedException ex) {
            Logger.getLogger(GUIFrame.class.getName()).log(Level.SEVERE, null, ex);
        }
        setIsRunning(false);
    }
}
