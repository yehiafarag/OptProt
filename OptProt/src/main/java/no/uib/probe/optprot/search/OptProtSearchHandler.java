package no.uib.probe.optprot.search;

import com.compomics.util.experiment.identification.Advocate;
import com.compomics.util.parameters.identification.IdentificationParameters;
import com.compomics.util.parameters.identification.tool_specific.MyriMatchParameters;
import com.compomics.util.parameters.identification.tool_specific.SageParameters;
import java.io.File;
import java.io.IOException;
import java.util.List;
import no.uib.probe.optprot.configurations.Configurations;
import no.uib.probe.optprot.dataset.model.SearchingSubDataset;
import no.uib.probe.optprot.model.SearchInputSetting;
import no.uib.probe.optprot.search.myrimatch.MyrimatchOptProtSearchOptimizer;
import no.uib.probe.optprot.search.sage.SageOptProtSearchOptimizer;
import no.uib.probe.optprot.search.xtandam.XTandemOptProtSearchOptimizer;

/**
 *
 * @author yfa041
 */
public class OptProtSearchHandler {

    public File optimizeSearchEngine(SearchingSubDataset searchingSubDataset, SearchInputSetting searchInputSetting, List<String> paramOrder) {
        try {
            final File generatedIdentificationParametersFile;
            IdentificationParameters identificationParameters = IdentificationParameters.getIdentificationParameters(searchingSubDataset.getSearchSettingsFile());
            if (searchInputSetting.isOptimizeAllParameters()) {
                generatedIdentificationParametersFile = new File(searchingSubDataset.getSubDataFolder(), Configurations.DEFAULT_RESULT_NAME + "_" + Configurations.DEFAULT_OPTPROT_SEARCH_SETTINGS_FILE_NAME.replace(".par", "_optAll.par"));
            } else {
                generatedIdentificationParametersFile = new File(searchingSubDataset.getSubDataFolder(), searchingSubDataset.getSearchSettingsFile().getName());
            }
            if (generatedIdentificationParametersFile.exists()) {
                generatedIdentificationParametersFile.delete();
            }
            generatedIdentificationParametersFile.createNewFile();
            IdentificationParameters.saveIdentificationParameters(identificationParameters, generatedIdentificationParametersFile);
            if (searchInputSetting.getSelectedSearchEngine().getIndex() == Advocate.xtandem.getIndex()) {
                XTandemOptProtSearchOptimizer xTandemOptProtSearchOptimizer = new XTandemOptProtSearchOptimizer(searchingSubDataset, searchInputSetting, generatedIdentificationParametersFile);
                xTandemOptProtSearchOptimizer.startProcess(paramOrder);
            } else if (searchInputSetting.getSelectedSearchEngine().getIndex() == Advocate.myriMatch.getIndex()) {
                MyriMatchParameters myriMatchParameters = (MyriMatchParameters) identificationParameters.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.myriMatch.getIndex());
               if(searchInputSetting.isOptimizeAllParameters()){
                myriMatchParameters.setMaxDynamicMods(4);
                myriMatchParameters.setNumberOfSpectrumMatches(1);
               }
                IdentificationParameters.saveIdentificationParameters(identificationParameters, generatedIdentificationParametersFile);
                MyrimatchOptProtSearchOptimizer myrimatchOptProtSearchOptimizer = new MyrimatchOptProtSearchOptimizer(searchingSubDataset, searchInputSetting, generatedIdentificationParametersFile);
                myrimatchOptProtSearchOptimizer.startProcess(paramOrder);
            } else if (searchInputSetting.getSelectedSearchEngine().getIndex() == Advocate.sage.getIndex()) {
                SageParameters myriMatchParameters = (SageParameters) identificationParameters.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());
               if(searchInputSetting.isOptimizeAllParameters()){
                myriMatchParameters.setMaxVariableMods(2);
               }
                IdentificationParameters.saveIdentificationParameters(identificationParameters, generatedIdentificationParametersFile);
                SageOptProtSearchOptimizer sageOptProtSearchOptimizer = new SageOptProtSearchOptimizer(searchingSubDataset, searchInputSetting, generatedIdentificationParametersFile);
                sageOptProtSearchOptimizer.startProcess(paramOrder);
                
            }
            return generatedIdentificationParametersFile;
        } catch (IOException ex) {
            ex.printStackTrace();
        }
        return null;

    }
}
