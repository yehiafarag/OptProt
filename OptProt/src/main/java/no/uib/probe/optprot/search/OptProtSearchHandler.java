package no.uib.probe.optprot.search;

import com.compomics.util.experiment.identification.Advocate;
import com.compomics.util.parameters.identification.IdentificationParameters;
import com.compomics.util.parameters.identification.tool_specific.MyriMatchParameters;
import com.compomics.util.parameters.identification.tool_specific.SageParameters;
import com.compomics.util.parameters.identification.tool_specific.XtandemParameters;
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
//            identificationParameters.getSearchParameters().getModificationParameters().clearFixedModifications();
//            identificationParameters.getSearchParameters().getModificationParameters().clearVariableModifications();
//            identificationParameters.getSearchParameters().getModificationParameters().clearRefinementModifications();
//            identificationParameters.getSearchParameters().getModificationParameters().getRefinementFixedModifications().clear();
//            IdentificationParameters.saveIdentificationParameters(identificationParameters, generatedIdentificationParametersFile);
            if (searchInputSetting.getSelectedSearchEngine().getIndex() == Advocate.xtandem.getIndex()) {
                XtandemParameters xtandemParameters = (XtandemParameters) identificationParameters.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
                xtandemParameters.setOutputResults("valid");//"valid"
                xtandemParameters.setMaxEValue(0.01);
                xtandemParameters.setProteinQuickAcetyl(false);
                xtandemParameters.setQuickPyrolidone(false);
                xtandemParameters.setStpBias(false);
                xtandemParameters.setRefine(false);
//                searchingSubDataset.setAcceptedIDRatioThreshold(0);
                IdentificationParameters.saveIdentificationParameters(identificationParameters, generatedIdentificationParametersFile);
                XTandemOptProtSearchOptimizer xTandemOptProtSearchOptimizer = new XTandemOptProtSearchOptimizer(searchingSubDataset, searchInputSetting, generatedIdentificationParametersFile);
                xTandemOptProtSearchOptimizer.startProcess(paramOrder);
            } else if (searchInputSetting.getSelectedSearchEngine().getIndex() == Advocate.myriMatch.getIndex()) {
                MyriMatchParameters myriMatchParameters = (MyriMatchParameters) identificationParameters.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.myriMatch.getIndex());
                myriMatchParameters.setMaxDynamicMods(4);
                myriMatchParameters.setNumberOfSpectrumMatches(1);
                searchingSubDataset.setAcceptedIDRatioThreshold(0);
                IdentificationParameters.saveIdentificationParameters(identificationParameters, generatedIdentificationParametersFile);
                MyrimatchOptProtSearchOptimizer myrimatchOptProtSearchOptimizer = new MyrimatchOptProtSearchOptimizer(searchingSubDataset, searchInputSetting, generatedIdentificationParametersFile);
                myrimatchOptProtSearchOptimizer.startProcess(paramOrder);
            } else if (searchInputSetting.getSelectedSearchEngine().getIndex() == Advocate.sage.getIndex()) {
                SageParameters myriMatchParameters = (SageParameters) identificationParameters.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());
                myriMatchParameters.setMaxVariableMods(4);
//                myriMatchParameters.setNumberOfSpectrumMatches(1);
                searchingSubDataset.setAcceptedIDRatioThreshold(0);
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
//    
//    public void printReport(){
//      System.out.println("-------------------------------------------------------------------------------------------");
//            System.out.println("Refrence score      :\t" + referenceIdRate);
//            System.out.println("Digestion           :\t" + generatedIdentificationParameters.getDigestionParameter());
//            System.out.println("Enzyme              :\t" + generatedIdentificationParameters.getEnzymeName());
//            System.out.println("Specificity         :\t" + generatedIdentificationParameters.getEnzymeSpecificity());
//            System.out.println("Max Missed Cleavages:\t" + generatedIdentificationParameters.getMaxMissedCleavage());
//            System.out.println("Fragment Ion Types  :\t" + generatedIdentificationParameters.getSelectedForwardIons() + "-" + generatedIdentificationParameters.getSelectedRewindIons());
//            System.out.println("Precursor Accuracy  :\t" + generatedIdentificationParameters.getPrecursorTolerance() + " ppm");
//            System.out.println("Fragment Accuracy   :\t" + generatedIdentificationParameters.getFragmentTolerance() + " Da");
//            System.out.println("PrecursorCharge     :\t" + generatedIdentificationParameters.getMinPrecursorCharge() + " - " + generatedIdentificationParameters.getMaxPrecursorCharge());
//            System.out.println("Isotops             :\t" + generatedIdentificationParameters.getMinIsotops() + " - " + generatedIdentificationParameters.getMaxIsotops());
//            System.out.println("default Variable mod:\t" + generatedIdentificationParameters.getVariableModifications() + "  Factor " + referenceFactor);
//            String fm = "";
//            if (generatedIdentificationParameters.getFixedModifications() != null) {
//                for (String fixedMod : generatedIdentificationParameters.getFixedModifications()) {
//                    fm += (fixedMod + "\n");
//                }
//            }
//            System.out.println("Fixed Mod mod:\n" + fm);
//            fm = "";
//            if (generatedIdentificationParameters.getSortedVariableModificationsMap() != null) {
//                for (int v : generatedIdentificationParameters.getSortedVariableModificationsMap().keySet()) {
//                    fm += (v + "% \t" + generatedIdentificationParameters.getSortedVariableModificationsMap().get(v).toString() + "\n");
//                }
//            }
//            System.out.println("Variable mod:\n" + fm);
//            if (optProtSearchParameters.isRunXTandem()) {
//                System.out.println("-----------------------------xtandem advanced-------------------------------------------");
//                System.out.println("Spectrum Dynamic Range:\t" + optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().getOptSpectrumDynamicRangeValue());
//                System.out.println("Number of Peaks       :\t" + optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().getOptSpectrumNumbPeaksValue());
//                System.out.println("MinimumFragmentMz     :\t" + optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().getOptSpectrumMinimumFragmentValue());
//                System.out.println("Minimum Peaks         :\t" + optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().getOptSpectrumMinimumPeaksValue());
//                System.out.println("Use NoiseSuppression  :\t" + optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().isOptNoiseSuppressionValue() + "  (" + optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().getOptMinimumPrecursorMassValue() + ")");
//                System.out.println("Use Parent isotop exp :\t" + optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().isOptParentIsotopeExpansionValue());
//                System.out.println("Use QuickAcetyl       :\t" + optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().isOptQuickAcetylValue());
//                System.out.println("Use QuickPyrolidone   :\t" + optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().isOptQuickPyrolidone());
//                System.out.println("Use stP Bias          :\t" + optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().isOptstPBiasValue());
//
//                System.out.println("Use Refinement        :\t" + optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().isOptUseRefine());
//                System.out.println("UnanticipatedCleavage :\t" + optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().isOptstPBiasValue());
//                System.out.println("SimiEnzymaticCleavage :\t" + optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().isOptstPBiasValue());
//                System.out.println("Potintial Modification:\t" + optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().isOptstPBiasValue());
//                System.out.println("Use PointMutations    :\t" + optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().isOptstPBiasValue());
//                System.out.println("Use SnAPs             :\t" + optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().isOptstPBiasValue());
//                System.out.println("Spectrum Synthesis    :\t" + optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().isOptstPBiasValue());
//
//                System.out.println("-------------------------------------------------------------------------------------------");
//                String rfm = "";
//                if (generatedIdentificationParameters.getRefinedFixedModifications() != null) {
//                    for (String fixedMod : generatedIdentificationParameters.getRefinedFixedModifications()) {
//                        rfm += (fixedMod + "\n");
//                    }
//                }
//                System.out.println("Refined Fixed Mod mod:\n" + rfm);
//                rfm = "";
//                if (generatedIdentificationParameters.getRefinedVariableModifications() != null) {
//                    for (String v : generatedIdentificationParameters.getRefinedVariableModifications()) {
//                        rfm += (v + "\n");
//                    }
//                }
//                System.out.println("Refined Variable mod:\n" + rfm);
//                System.out.println("-------------------------------------------------------------------------------------------");
//            }
//    }
}
