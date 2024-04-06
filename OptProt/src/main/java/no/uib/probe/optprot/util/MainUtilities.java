/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package no.uib.probe.optprot.util;

import com.compomics.util.parameters.UtilitiesUserParameters;
import com.compomics.util.parameters.searchgui.OutputParameters;
import com.compomics.util.parameters.tools.ProcessingParameters;
import eu.isas.searchgui.SearchHandler;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.logging.Level;
import java.util.logging.Logger;
import no.uib.probe.optprot.configurations.Configurations;

/**
 *
 * @author yfa041
 */
public class MainUtilities {

    static {
        UtilitiesUserParameters userParameters = UtilitiesUserParameters.loadUserParameters();
        userParameters.setGzip(false);
        userParameters.setSearchGuiOutputParameters(OutputParameters.no_zip);
        userParameters.setRenameXTandemFile(true);
        UtilitiesUserParameters.saveUserParameters(userParameters);
        SearchHandler.setCloseProcessWhenDone(false);
        File resultsOutput = new File(Configurations.OUTPUT_FOLDER_PATH);
        if (resultsOutput.exists()) {
//            cleanFolder(resultsOutput);
        }
        resultsOutput.mkdir();
    }
    private static final ProcessingParameters Processing_Parameters = new ProcessingParameters();
    public static final OptProtWaitingHandler OptProt_Waiting_Handler = new OptProtWaitingHandler();

    public static ProcessingParameters getProcessingParameter() {
        if (Processing_Parameters == null) {
            Processing_Parameters.setnThreads(Runtime.getRuntime().availableProcessors());
            // Processing
            Processing_Parameters.setnThreads(15);

        }
        return Processing_Parameters;
    }

    public static void cleanOutputFolder() {
       
            File outputFolder = new File(Configurations.OUTPUT_FOLDER_PATH);
            deleteFolder(outputFolder);
            outputFolder.mkdir();
       
    }


    public  static void deleteFolder(File folder) {
        if (folder.exists() && folder.isDirectory()) {           
            for (File f : folder.listFiles()) {
                if (f.isDirectory()) {                    
                    deleteFolder(f);
                } else {
                    f.delete();
                }

            }

        }
            folder.delete();
     

    }

    public static double rundDouble(double args) {
        return Math.round(args * 100.0) / 100.0;

    }

}
