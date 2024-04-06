/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package no.uib.probe.optprot.search;

import com.compomics.util.experiment.biology.modifications.ModificationCategory;
import com.compomics.util.experiment.biology.modifications.ModificationFactory;
import com.compomics.util.parameters.identification.IdentificationParameters;
import eu.isas.searchgui.SearchHandler;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import no.uib.probe.optprot.configurations.Configurations;
import no.uib.probe.optprot.model.OptProtSearchParameters;
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
     * @param searchEngineParameters
     * @param msFile
     * @param fastaFile
     * @param tempIdParam
     * @param identificationParametersFile
     * @return results Folder
     */
    public synchronized static File executeSearch(String processId, OptProtSearchParameters searchEngineParameters, File msFile, File fastaFile, IdentificationParameters tempIdParam, File identificationParametersFile) {

        if (searchEngineParameters.isRunNovor()||searchEngineParameters.isRunDirecTag()) {
            //remove terminal variable modifications andd add common before run novor         
            List<String> toRemoveMod = new ArrayList<>();
            for (String mod : tempIdParam.getSearchParameters().getModificationParameters().getVariableModifications()) {
                if (ptmFactory.getModification(mod).getModificationType().isNTerm() || ptmFactory.getModification(mod).getModificationType().isCTerm()) {
                    System.out.println("terminal : "+mod);
                    toRemoveMod.add(mod);
                }
            }
            for (String mod : toRemoveMod) {
                tempIdParam.getSearchParameters().getModificationParameters().removeVariableModification(mod);
            }
        }
        File resultOutput = new File(Configurations.OUTPUT_FOLDER_PATH, processId);
        resultOutput.mkdir();
        File tempSearchEngineFolder = new File(Configurations.OUTPUT_FOLDER_PATH, processId + "_temp");
        tempSearchEngineFolder.mkdir();
        
        ArrayList<File> msFileInList = new ArrayList<>();
        msFileInList.add(msFile);
        SearchHandler.setTempSearchEngineFolderPath(tempSearchEngineFolder.getAbsolutePath());
        SearchHandler searchHandler = new SearchHandler(tempIdParam, resultOutput, tempSearchEngineFolder, processId, msFileInList,
                fastaFile, new ArrayList<>(),
                identificationParametersFile,
                searchEngineParameters.isRunOmssa(),
                searchEngineParameters.isRunXTandem(),
                searchEngineParameters.isRunMsgf(),
                searchEngineParameters.isRunMsAmanda(),
                searchEngineParameters.isRunMyriMatch(),
                searchEngineParameters.isRunComet(),
                searchEngineParameters.isRunTide(),
                searchEngineParameters.isRunAndromeda(),
                searchEngineParameters.isRunMetaMorpheus(),
                searchEngineParameters.isRunSage(),
                searchEngineParameters.isRunNovor(),
                searchEngineParameters.isRunDirecTag(),
                searchEngineParameters.getOmssaFolder(),
                searchEngineParameters.getxTandemFolder(),
                searchEngineParameters.getMsgfFolder(),
                searchEngineParameters.getMsAmandaFolder(),
                searchEngineParameters.getMyriMatchFolder(),
                searchEngineParameters.getCometFolder(),
                searchEngineParameters.getTideFolder(),
                searchEngineParameters.getTideIndexLocation(),
                searchEngineParameters.getAndromedaFolder(),
                searchEngineParameters.getMetaMorpheusFolder(),
                searchEngineParameters.getSageFolder(),
                searchEngineParameters.getNovorFolder(),
                searchEngineParameters.getDirecTagFolder(),
                searchEngineParameters.getMakeblastdbFolder(),
                MainUtilities.getProcessingParameter()
        );
        try {
         
            searchHandler.startSearch(MainUtilities.OptProt_Waiting_Handler);
        } catch (InterruptedException ex) {
            ex.printStackTrace();
        }
        File resultsFile = searchHandler.getResultsFolder();
        return resultsFile;

    }
}
