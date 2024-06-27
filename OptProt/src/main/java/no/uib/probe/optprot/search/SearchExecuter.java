/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package no.uib.probe.optprot.search;

import com.compomics.util.experiment.biology.modifications.ModificationFactory;
import com.compomics.util.experiment.identification.Advocate;
import com.compomics.util.parameters.identification.IdentificationParameters;
import com.compomics.util.parameters.identification.tool_specific.XtandemParameters;
import eu.isas.searchgui.SearchHandler;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import no.uib.probe.optprot.configurations.Configurations;
import no.uib.probe.optprot.model.SearchInputSetting;
import no.uib.probe.optprot.util.MainUtilities;

/**
 *
 * @author yfa041
 */
public class SearchExecuter {

    /**
     * The compomics PTM factory.
     */
    private static final ModificationFactory ptmFactory = ModificationFactory.getInstance();

    /**
     *
     * @param processId
     * @param searchInputSetting
     * @param msFile
     * @param fastaFile
     * @param tempIdParam
     * @param identificationParametersFile
     * @return results Folder
     */
    public synchronized static File executeSearch(String processId, SearchInputSetting searchInputSetting, File msFile, File fastaFile, IdentificationParameters tempIdParam, File identificationParametersFile) {

        if (searchInputSetting.isRunNovor() || searchInputSetting.isRunDirecTag()) {
            //remove terminal variable modifications andd add common before run novor         
            List<String> toRemoveMod = new ArrayList<>();
            for (String mod : tempIdParam.getSearchParameters().getModificationParameters().getVariableModifications()) {
                if (ptmFactory.getModification(mod).getModificationType().isNTerm() || ptmFactory.getModification(mod).getModificationType().isCTerm()) {
                    System.out.println("terminal : " + mod);
                    toRemoveMod.add(mod);
                }
            }
            for (String mod : toRemoveMod) {
                tempIdParam.getSearchParameters().getModificationParameters().removeVariableModification(mod);
            }
        }
//        if (searchInputSetting.isRunXTandem()) {
//            XtandemParameters xtandemParameters = (XtandemParameters) tempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
////            if (processId.contains("init_input_files")) {
////                xtandemParameters.setOutputResults("valid");//"valid"
////                xtandemParameters.setMaxEValue(0.01);
////                xtandemParameters.setProteinQuickAcetyl(false);
////                xtandemParameters.setQuickPyrolidone(false);
////                xtandemParameters.setStpBias(false);
////                xtandemParameters.setRefine(false);
////            } else {
////                xtandemParameters.setOutputResults("all");//"valid"
////            }
//////            xtandemParameters.setMaxEValue(0.01);
////            System.out.println("process id " + processId);
////            if (processId.contains("reference_run_")) {
////                xtandemParameters.setProteinQuickAcetyl(false);
////                xtandemParameters.setQuickPyrolidone(false);
////                xtandemParameters.setStpBias(false);
////                xtandemParameters.setRefine(false);
////            }
////            
//        }
        File resultOutput = new File(Configurations.OUTPUT_FOLDER_PATH, processId);
        resultOutput.mkdir();
        File tempSearchEngineFolder = new File(Configurations.OUTPUT_FOLDER_PATH, processId + "_temp");
        tempSearchEngineFolder.mkdir();

        ArrayList<File> msFileInList = new ArrayList<>();
        msFileInList.add(msFile);
        SearchHandler.setTempSearchEngineFolderPath(tempSearchEngineFolder.getAbsolutePath());
        SearchHandler.setTempFolderPath(null);
        SearchHandler searchHandler = new SearchHandler(tempIdParam, resultOutput, tempSearchEngineFolder, processId, msFileInList,
                fastaFile, new ArrayList<>(),
                identificationParametersFile,
                searchInputSetting.isRunOmssa(),
                searchInputSetting.isRunXTandem(),
                searchInputSetting.isRunMsgf(),
                searchInputSetting.isRunMsAmanda(),
                searchInputSetting.isRunMyriMatch(),
                searchInputSetting.isRunComet(),
                searchInputSetting.isRunTide(),
                searchInputSetting.isRunAndromeda(),
                searchInputSetting.isRunMetaMorpheus(),
                searchInputSetting.isRunSage(),
                searchInputSetting.isRunNovor(),
                searchInputSetting.isRunDirecTag(),
                searchInputSetting.getOmssaFolder(),
                searchInputSetting.getxTandemFolder(),
                searchInputSetting.getMsgfFolder(),
                searchInputSetting.getMsAmandaFolder(),
                searchInputSetting.getMyriMatchFolder(),
                searchInputSetting.getCometFolder(),
                searchInputSetting.getTideFolder(),
                searchInputSetting.getTideIndexLocation(),
                searchInputSetting.getAndromedaFolder(),
                searchInputSetting.getMetaMorpheusFolder(),
                searchInputSetting.getSageFolder(),
                searchInputSetting.getNovorFolder(),
                searchInputSetting.getDirecTagFolder(),
                searchInputSetting.getMakeblastdbFolder(),
                MainUtilities.getProcessingParameter()
        );

        try {
            searchHandler.startSearch(MainUtilities.OptProt_Waiting_Handler);
        } catch (InterruptedException ex) {
            ex.printStackTrace();
        }
        File resultsFile = searchHandler.getResultsFolder();
//        MainUtilities.deleteFolder(tempSearchEngineFolder);
//        System.exit(0);

        return resultsFile;

    }
}
