package no.uib.probe.optprot.arc;

import com.compomics.util.experiment.biology.enzymes.Enzyme;
import com.compomics.util.experiment.biology.enzymes.EnzymeFactory;
import com.compomics.util.experiment.biology.ions.impl.PeptideFragmentIon;
import com.compomics.util.experiment.biology.modifications.ModificationCategory;
import com.compomics.util.experiment.biology.modifications.ModificationFactory;
import com.compomics.util.experiment.biology.proteins.Peptide;
import com.compomics.util.experiment.identification.Advocate;
import com.compomics.util.experiment.identification.matches.SpectrumMatch;
import com.compomics.util.experiment.identification.modification.search_engine_mapping.ModificationLocalizationMapper;
import com.compomics.util.experiment.identification.protein_inference.fm_index.FMIndex;
import com.compomics.util.experiment.identification.spectrum_assumptions.PeptideAssumption;
import com.compomics.util.experiment.identification.spectrum_assumptions.TagAssumption;
import com.compomics.util.experiment.io.identification.IdfileReader;
import com.compomics.util.experiment.io.identification.IdfileReaderFactory;
import com.compomics.util.experiment.io.mass_spectrometry.MsFileHandler;
import com.compomics.util.io.IoUtil;
import com.compomics.util.parameters.identification.IdentificationParameters;
import com.compomics.util.parameters.identification.advanced.SequenceMatchingParameters;
import com.compomics.util.parameters.identification.search.DigestionParameters;
import com.compomics.util.parameters.identification.search.SearchParameters;
import com.compomics.util.parameters.identification.tool_specific.DirecTagParameters;
import com.compomics.util.parameters.identification.tool_specific.MyriMatchParameters;
import com.compomics.util.parameters.identification.tool_specific.XtandemParameters;
import eu.isas.searchgui.SearchHandler;
import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.concurrent.Future;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.xml.bind.JAXBException;
import javax.xml.stream.XMLStreamException;
import no.uib.probe.optprot.configurations.Configurations;
import no.uib.probe.optprot.model.OptimisedSearchResults;
import no.uib.probe.optprot.model.SearchInputSetting;
import no.uib.probe.optprot.search.SearchExecuter;
import no.uib.probe.optprot.util.MainUtilities;
import no.uib.probe.optprot.util.OptProtWaitingHandler;
import no.uib.probe.optprot.util.ReportExporter;
import no.uib.probe.optprot.util.SpectraUtilities;
import org.xmlpull.v1.XmlPullParserException;

/**
 * This class responsible for performing sub-searches in order to optimize
 * different search parameters
 *
 * @author Yehia Mokhtar Farag
 */
public class SearchOptimizerHandler {

    /**
     * The identification settings file.
     */
    private File identificationParametersFile;
    private File optimizedIdentificationParametersFile;
    private File subMsFile;
    private ArrayList<File> msFileInList;
    private File subFastaFile;
    private File resultsOutput;
//    private ExecutorService executor;
    /**
     * The identification file reader factory of compomics utilities.
     */
    private final IdfileReaderFactory readerFactory = IdfileReaderFactory.getInstance();
    private SearchInputSetting optProtSearchParameters;
    private final OptimisedSearchResults optimisedSearchParameter;
    /**
     * The compomics PTM factory.
     */
    private final ModificationFactory ptmFactory = ModificationFactory.getInstance();

    public SearchOptimizerHandler() {
//        executor = Executors.newFixedThreadPool(1);
        this.optimisedSearchParameter = new OptimisedSearchResults();
    }

    public void setSubMsFile(File msFile) {
        this.subMsFile = msFile;
        this.msFileInList = new ArrayList<>();
        this.msFileInList.add(msFile);
    }

    public void setSubFastaFile(File subFastaFile) {
        this.subFastaFile = subFastaFile;
    }

    public void setIdentificationParametersFile(File identificationParametersFile) {
        try {
            this.identificationParametersFile = identificationParametersFile;
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);

            optimizedIdentificationParametersFile = new File(identificationParametersFile.getParent(), Configurations.DEFAULT_RESULT_NAME + "_" + identificationParametersFile.getName());
            if (optimizedIdentificationParametersFile.exists()) {
                optimizedIdentificationParametersFile.delete();
            }
            optimizedIdentificationParametersFile.createNewFile();
            IdentificationParameters.saveIdentificationParameters(tempIdParam, optimizedIdentificationParametersFile);

            optimisedSearchParameter.setDigestionParameter(tempIdParam.getSearchParameters().getDigestionParameters().getCleavageParameter().name());
            if (optimisedSearchParameter.getDigestionParameter().equalsIgnoreCase(DigestionParameters.CleavageParameter.enzyme.name())) {
                optimisedSearchParameter.setEnzymeName(tempIdParam.getSearchParameters().getDigestionParameters().getEnzymes().get(0).getName());
                optimisedSearchParameter.setEnzymeSpecificity(tempIdParam.getSearchParameters().getDigestionParameters().getSpecificity(optimisedSearchParameter.getEnzymeName()).name());
                optimisedSearchParameter.setMaxMissedCleavage(tempIdParam.getSearchParameters().getDigestionParameters().getnMissedCleavages(optimisedSearchParameter.getEnzymeName()));
            }
            optimisedSearchParameter.setSelectedForwardIons(tempIdParam.getSearchParameters().getForwardIons());
            optimisedSearchParameter.setSelectedRewindIons(tempIdParam.getSearchParameters().getRewindIons());
            optimisedSearchParameter.setPrecursorTolerance(tempIdParam.getSearchParameters().getPrecursorAccuracy());
            optimisedSearchParameter.setFragmentTolerance(tempIdParam.getSearchParameters().getFragmentIonAccuracyInDaltons());

            optimisedSearchParameter.setMaxPrecursorCharge(tempIdParam.getSearchParameters().getMaxChargeSearched());
            optimisedSearchParameter.setMinPrecursorCharge(tempIdParam.getSearchParameters().getMinChargeSearched());

            optimisedSearchParameter.setMaxIsotops(tempIdParam.getSearchParameters().getMaxIsotopicCorrection());
            optimisedSearchParameter.setMinIsotops(tempIdParam.getSearchParameters().getMinIsotopicCorrection());

            for (String mod : tempIdParam.getSearchParameters().getModificationParameters().getVariableModifications()) {

                optimisedSearchParameter.addVariableModifications(mod, 0.0);
            }
            for (String mod : tempIdParam.getSearchParameters().getModificationParameters().getFixedModifications()) {
                optimisedSearchParameter.addFixedModifications(mod);
            }

        } catch (IOException ex) {
            Logger.getLogger(SearchOptimizerHandler.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    public void executeParametersOptimization() {
        try {

            IdentificationParameters optimizedIdentificationParameters = IdentificationParameters.getIdentificationParameters(optimizedIdentificationParametersFile);
            if (optProtSearchParameters.isRecalibrateSpectraParameter()) {
                File updatedMgfFile = this.runReferenceSearch(optProtSearchParameters.isRecalibrateSpectraParameter());
                System.out.println("before calibration: " + referenceIdRate);
                this.subMsFile = updatedMgfFile;
                this.runReferenceSearch(false);
                System.out.println("after calibration: " + referenceIdRate);
                optimizedIdentificationParameters.getSearchParameters().getDigestionParameters().setCleavageParameter(DigestionParameters.CleavageParameter.valueOf(optimisedSearchParameter.getDigestionParameter()));
//                IdentificationParameters.saveIdentificationParameters(optimizedIdentificationParameters, optimizedIdentificationParametersFile);
            }
//           
            if (optProtSearchParameters.isOptimizeDigestionParameter()) {
                this.optimizeDigestionParameter();
                optimizedIdentificationParameters.getSearchParameters().getDigestionParameters().setCleavageParameter(DigestionParameters.CleavageParameter.valueOf(optimisedSearchParameter.getDigestionParameter()));
//                IdentificationParameters.saveIdentificationParameters(optimizedIdentificationParameters, optimizedIdentificationParametersFile);
            }
//            if (optimisedSearchParameter.getDigestionParameter().equalsIgnoreCase(DigestionParameters.CleavageParameter.enzyme.name()) && (optProtSearchParameters.isOptimizeDigestionParameter() || optProtSearchParameters.isOptimizeEnzymeParameter())) {
//                this.optimizeEnzymeParameter();
//                optimizedIdentificationParameters.getSearchParameters().getDigestionParameters().getEnzymes().clear();
//                optimizedIdentificationParameters.getSearchParameters().getDigestionParameters().getEnzymes().add(EnzymeFactory.getInstance().getEnzyme(optimisedSearchParameter.getEnzymeName()));
//                optimizedIdentificationParameters.getSearchParameters().getDigestionParameters().setnMissedCleavages(optimisedSearchParameter.getEnzymeName(), optimisedSearchParameter.getMaxMissedCleavage());
////                IdentificationParameters.saveIdentificationParameters(optimizedIdentificationParameters, optimizedIdentificationParametersFile);
//            }
//            if (optimisedSearchParameter.getDigestionParameter().equalsIgnoreCase(DigestionParameters.CleavageParameter.enzyme.name()) && (optProtSearchParameters.isOptimizeDigestionParameter() || optProtSearchParameters.isOptimizeSpecificityParameter())) {
//                this.optimizeSpecificityParameter();
//                optimizedIdentificationParameters.getSearchParameters().getDigestionParameters().setSpecificity(optimisedSearchParameter.getEnzymeName(), DigestionParameters.Specificity.valueOf(optimisedSearchParameter.getEnzymeSpecificity()));
////                IdentificationParameters.saveIdentificationParameters(optimizedIdentificationParameters, optimizedIdentificationParametersFile);
//            }
//            if ((optimisedSearchParameter.getDigestionParameter().equalsIgnoreCase(DigestionParameters.CleavageParameter.enzyme.name()) && (optProtSearchParameters.isOptimizeDigestionParameter()) || optProtSearchParameters.isOptimizeMaxMissCleavagesParameter())) {
//                this.optimizeMaxMissCleavagesParameter();
//                optimizedIdentificationParameters.getSearchParameters().getDigestionParameters().setnMissedCleavages(optimisedSearchParameter.getEnzymeName(), optimisedSearchParameter.getMaxMissedCleavage());
//
////                IdentificationParameters.saveIdentificationParameters(optimizedIdentificationParameters, optimizedIdentificationParametersFile);
//            }
            if (optProtSearchParameters.isOptimizeFragmentIonTypesParameter()) {
                this.optimizeFragmentIonTypesParameter();
                optimizedIdentificationParameters.getSearchParameters().setForwardIons(optimisedSearchParameter.getSelectedForwardIons());
                optimizedIdentificationParameters.getSearchParameters().setRewindIons(optimisedSearchParameter.getSelectedRewindIons());
//                IdentificationParameters.saveIdentificationParameters(optimizedIdentificationParameters, optimizedIdentificationParametersFile);
            }
            if (this.optProtSearchParameters.isOptimizePrecursorToleranceParameter()) {
                this.optimizePrecursorToleranceParameter();
                optimizedIdentificationParameters.getSearchParameters().setPrecursorAccuracy(optimisedSearchParameter.getPrecursorTolerance());
                optimizedIdentificationParameters.getSearchParameters().setPrecursorAccuracyType(SearchParameters.MassAccuracyType.PPM);
//                IdentificationParameters.saveIdentificationParameters(optimizedIdentificationParameters, optimizedIdentificationParametersFile);
            }
            if (this.optProtSearchParameters.isOptimizeFragmentToleranceParameter()) {
                this.optimizeFragmentToleranceParameter();
                optimizedIdentificationParameters.getSearchParameters().setFragmentIonAccuracy(optimisedSearchParameter.getFragmentTolerance());
                optimizedIdentificationParameters.getSearchParameters().setFragmentAccuracyType(SearchParameters.MassAccuracyType.DA);
//                IdentificationParameters.saveIdentificationParameters(optimizedIdentificationParameters, optimizedIdentificationParametersFile);
            }

            if (this.optProtSearchParameters.isOptimizePrecursorChargeParameter()) {
                this.optimizePrecursorChargeParameter();
                optimizedIdentificationParameters.getSearchParameters().setMinChargeSearched(optimisedSearchParameter.getMinPrecursorCharge());
                optimizedIdentificationParameters.getSearchParameters().setMaxChargeSearched(optimisedSearchParameter.getMaxPrecursorCharge());
//                IdentificationParameters.saveIdentificationParameters(optimizedIdentificationParameters, optimizedIdentificationParametersFile);
            }
            if (this.optProtSearchParameters.isOptimizeIsotopsParameter()) {
                this.optimizeIsotopParameter();
                optimizedIdentificationParameters.getSearchParameters().setMaxIsotopicCorrection(optimisedSearchParameter.getMinIsotops());
                optimizedIdentificationParameters.getSearchParameters().setMaxIsotopicCorrection(optimisedSearchParameter.getMaxIsotops());
//                IdentificationParameters.saveIdentificationParameters(optimizedIdentificationParameters, optimizedIdentificationParametersFile);
            }
            if (this.optProtSearchParameters.isOptimizeModificationParameter()) {
                //create reference data
//                this.optimizeModificationsParameter();
//            
                this.optimizeModificationsParameter();
                optimizedIdentificationParameters = IdentificationParameters.getIdentificationParameters(optimizedIdentificationParametersFile);
                optimizedIdentificationParameters.getSearchParameters().getModificationParameters().clearFixedModifications();
                optimizedIdentificationParameters.getSearchParameters().getModificationParameters().clearVariableModifications();
                optimizedIdentificationParameters.getSearchParameters().getModificationParameters().clearRefinementModifications();
                optimizedIdentificationParameters.getSearchParameters().getModificationParameters().getRefinementFixedModifications().clear();
                optimizedIdentificationParameters.getSearchParameters().getModificationParameters().getRefinementVariableModifications().clear();
                for (String modId : optimisedSearchParameter.getSortedVariableModificationsMap().firstEntry().getValue()) {
                    optimizedIdentificationParameters.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(modId));
                }
                for (String modId : optimisedSearchParameter.getFixedModifications()) {
                    optimizedIdentificationParameters.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(modId));
                }
                for (String modId : optimisedSearchParameter.getRefinedFixedModifications()) {
                    optimizedIdentificationParameters.getSearchParameters().getModificationParameters().addRefinementFixedModification(ptmFactory.getModification(modId));
                }
                for (String modId : optimisedSearchParameter.getRefinedVariableModifications()) {
                    optimizedIdentificationParameters.getSearchParameters().getModificationParameters().addRefinementVariableModification(ptmFactory.getModification(modId));
                }
                IdentificationParameters.saveIdentificationParameters(optimizedIdentificationParameters, optimizedIdentificationParametersFile);
            }

            if (optProtSearchParameters.isRunXTandem() && optProtSearchParameters.isOptimizeXtandemAdvancedParameter()) {

                XtandemAdvancedSearchOptimizerHandler xtendemHandler = new XtandemAdvancedSearchOptimizerHandler(subMsFile, subFastaFile, optimizedIdentificationParametersFile, optProtSearchParameters, referenceIdRate);
                if (optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().isOptSpectrumDynamicRange()) {
                    optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().setOptSpectrumDynamicRangeValue(xtendemHandler.optimizeSpectrumDynamicRange());

                }
                if (optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().isOptSpectrumNumbPeaks()) {
                    optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().setOptSpectrumNumbPeaksValue(xtendemHandler.optimizeSpectrumPeaksNumber());

                }
                if (optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().isOptSpectrumMinimumFragment()) {
                    optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().setOptSpectrumMinimumFragmentValue(xtendemHandler.optimizeMinimumFragmentMz());

                }
                if (optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().isOptSpectrumMinimumPeaks()) {
                    optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().setOptSpectrumMinimumPeaksValue(xtendemHandler.optimizeMinimumPeaks());

                }
                if (optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().isOptNoiseSuppression()) {
                    double minPrecursorMass = xtendemHandler.optimizeNoiseSuppression();
                    optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().setOptNoiseSuppressionValue(minPrecursorMass == -1);
                    if (minPrecursorMass != -1) {
                        optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().setOptMinimumPrecursorMassValue(xtendemHandler.optimizeNoiseSuppression());
                    }

                }
                if (optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().isOptParentIsotopeExpansion()) {
                    optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().setOptParentIsotopeExpansionValue(xtendemHandler.optimizeParentIsotopExpansion());
                }
                if (optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().isOptQuickAcetyl()) {
                    optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().setOptQuickAcetylValue(xtendemHandler.optimizeParentIsotopExpansion());
                }
                if (optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().isOptQuickPyrolidone()) {
                    optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().setOptQuickPyrolidoneValue(xtendemHandler.optimizeParentIsotopExpansion());
                }
                if (optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().isOptstPBias()) {
                    optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().setOptstPBiasValue(xtendemHandler.optimizeParentIsotopExpansion());
                }

//                if (xtendemHandler.optimizeUseRefine()) {
                if (optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().isOptRefineUnanticipatedCleavage()) {
                    optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().setOptRefineUnanticipatedCleavage(xtendemHandler.optimizeRefineUnanticipatedCleavage());
                }
                if (optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().isOptRefineSimiEnzymaticCleavage()) {
                    optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().setOptRefineSimiEnzymaticCleavageValue(xtendemHandler.optimizeRefineSimiEnzymaticCleavage());
                }
                if (optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().isOptRefinePotintialModification()) {
                    optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().setOptRefinePotintialModificationValue(xtendemHandler.optimizeRefinePotintialModification());
                }
                if (optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().isOptRefinePointMutations()) {
                    optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().setOptRefinePointMutationsValue(xtendemHandler.optimizeRefinePointMutations());
                }
                if (optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().isOptRefineSnAPs()) {
                    optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().setOptRefineSnAPsValue(xtendemHandler.optimizeRefineSnAPs());
                }
                if (optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().isOptRefineSpectrumSynthesis()) {
                    optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().setOptRefineSpectrumSynthesisValue(xtendemHandler.optimizeRefineSpectrumSynthesis());
                }
                if (optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().isOptUseRefine()) {
                    optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().setOptUseRefine(xtendemHandler.optimizeUseRefine());

                }
//                }

            }
        } catch (IOException ex) {
            ex.printStackTrace();
        } finally {

////            executor.shutdown();
//            while (!executor.isTerminated()) {
//            }
//            MainUtilities.cleanFolder(resultsOutput);

            System.out.println("-------------------------------------------------------------------------------------------");
            System.out.println("Refrence score      :\t" + referenceIdRate);
            System.out.println("Digestion           :\t" + optimisedSearchParameter.getDigestionParameter());
            System.out.println("Enzyme              :\t" + optimisedSearchParameter.getEnzymeName());
            System.out.println("Specificity         :\t" + optimisedSearchParameter.getEnzymeSpecificity());
            System.out.println("Max Missed Cleavages:\t" + optimisedSearchParameter.getMaxMissedCleavage());
            System.out.println("Fragment Ion Types  :\t" + optimisedSearchParameter.getSelectedForwardIons() + "-" + optimisedSearchParameter.getSelectedRewindIons());
            System.out.println("Precursor Accuracy  :\t" + optimisedSearchParameter.getPrecursorTolerance() + " ppm");
            System.out.println("Fragment Accuracy   :\t" + optimisedSearchParameter.getFragmentTolerance() + " Da");
            System.out.println("PrecursorCharge     :\t" + optimisedSearchParameter.getMinPrecursorCharge() + " - " + optimisedSearchParameter.getMaxPrecursorCharge());
            System.out.println("Isotops             :\t" + optimisedSearchParameter.getMinIsotops() + " - " + optimisedSearchParameter.getMaxIsotops());
            System.out.println("default Variable mod:\t" + optimisedSearchParameter.getVariableModifications() + "  Factor " + referenceFactor);
            String fm = "";
            if (optimisedSearchParameter.getFixedModifications() != null) {
                for (String fixedMod : optimisedSearchParameter.getFixedModifications()) {
                    fm += (fixedMod + "\n");
                }
            }
            System.out.println("Fixed Mod mod:\n" + fm);
            fm = "";
            if (optimisedSearchParameter.getSortedVariableModificationsMap() != null) {
                for (int v : optimisedSearchParameter.getSortedVariableModificationsMap().keySet()) {
                    fm += (v + "% \t" + optimisedSearchParameter.getSortedVariableModificationsMap().get(v).toString() + "\n");
                }
            }
            System.out.println("Variable mod:\n" + fm);
            if (optProtSearchParameters.isRunXTandem()) {
                System.out.println("-----------------------------xtandem advanced-------------------------------------------");
                System.out.println("Spectrum Dynamic Range:\t" + optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().getOptSpectrumDynamicRangeValue());
                System.out.println("Number of Peaks       :\t" + optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().getOptSpectrumNumbPeaksValue());
                System.out.println("MinimumFragmentMz     :\t" + optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().getOptSpectrumMinimumFragmentValue());
                System.out.println("Minimum Peaks         :\t" + optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().getOptSpectrumMinimumPeaksValue());
                System.out.println("Use NoiseSuppression  :\t" + optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().isOptNoiseSuppressionValue() + "  (" + optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().getOptMinimumPrecursorMassValue() + ")");
                System.out.println("Use Parent isotop exp :\t" + optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().isOptParentIsotopeExpansionValue());
                System.out.println("Use QuickAcetyl       :\t" + optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().isOptQuickAcetylValue());
                System.out.println("Use QuickPyrolidone   :\t" + optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().isOptQuickPyrolidone());
                System.out.println("Use stP Bias          :\t" + optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().isOptstPBiasValue());

                System.out.println("Use Refinement        :\t" + optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().isOptUseRefine());
                System.out.println("UnanticipatedCleavage :\t" + optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().isOptstPBiasValue());
                System.out.println("SimiEnzymaticCleavage :\t" + optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().isOptstPBiasValue());
                System.out.println("Potintial Modification:\t" + optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().isOptstPBiasValue());
                System.out.println("Use PointMutations    :\t" + optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().isOptstPBiasValue());
                System.out.println("Use SnAPs             :\t" + optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().isOptstPBiasValue());
                System.out.println("Spectrum Synthesis    :\t" + optProtSearchParameters.getXtandemOptProtAdvancedSearchParameters().isOptstPBiasValue());

                System.out.println("-------------------------------------------------------------------------------------------");
                String rfm = "";
                if (optimisedSearchParameter.getRefinedFixedModifications() != null) {
                    for (String fixedMod : optimisedSearchParameter.getRefinedFixedModifications()) {
                        rfm += (fixedMod + "\n");
                    }
                }
                System.out.println("Refined Fixed Mod mod:\n" + rfm);
                rfm = "";
                if (optimisedSearchParameter.getRefinedVariableModifications() != null) {
                    for (String v : optimisedSearchParameter.getRefinedVariableModifications()) {
                        rfm += (v + "\n");
                    }
                }
                System.out.println("Refined Variable mod:\n" + rfm);
                System.out.println("-------------------------------------------------------------------------------------------");
            }
        }

    }

    public void setOptProtSearchParameters(SearchInputSetting optProtSearchParameters) {
        this.optProtSearchParameters = optProtSearchParameters;
    }
    private double referenceIdRate;
    private double referenceFactor;

    public synchronized double getReferenceIdRate() {

        return referenceIdRate;
    }

    public File runReferenceSearch(boolean recalibratSpectra) throws IOException {
        IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        final String option = "reference_run";
        final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option;
        ArrayList<SpectrumMatch> spectrumMatches = new ArrayList<>();

        Future f = MainUtilities.getExecutorService().submit(() -> {
            spectrumMatches.addAll(excuteSearch(updatedName, "Reference Run", option, tempIdParam, recalibratSpectra));

        });
        while (!f.isDone()) {
        }

        referenceIdRate = spectrumMatches.size();
        if (recalibratSpectra) {
            return SpectraUtilities.recalibrateSpectra(subMsFile, subFastaFile, spectrumMatches, tempIdParam);
        }

//        List<String> mods = new ArrayList<>(tempIdParam.getSearchParameters().getModificationParameters().getVariableModifications());
//        for (String mod : mods) {
//            tempIdParam.getSearchParameters().getModificationParameters().removeVariableModification(mod);
//        }
//        f = executor.submit(() -> {
//            referenceFactor = 0.20 * (referenceIdRate - countSpectra(excuteSearch(updatedName + "_factor", "Reference Factor Run", option + "_factor", tempIdParam)));
//        });
//        while (!f.isDone()) {
//        }
//        System.out.println("job is done- reference factor " + referenceFactor);
//        for (String mod : mods) {
//            System.out.println("param v ptm  " + mod);
//            optimisedSearchParameter.addVariableModifications(mod, referenceIdRate);
//        }
        return null;
    }

    private void optimizeDigestionParameter() throws IOException {
        long start1 = System.currentTimeMillis();
        Map<String, Double> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        double idRate = -1;
        String selectedOption = "";
        String[] cleavageParameters = new String[]{"enzyme", "wholeProtein", "unSpecific"};
        // Creates a pool of 3 threads  
        for (int i = 0; i < cleavageParameters.length; i++) {
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            tempIdParam.getSearchParameters().getDigestionParameters().setCleavageParameter(DigestionParameters.CleavageParameter.valueOf(cleavageParameters[i]));
            final String option = cleavageParameters[i];
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option;

            Future f = MainUtilities.getExecutorService().submit(() -> {
                resultsMap.put(option, countSpectra(excuteSearch(updatedName, "Cleavage_Parameters", option, tempIdParam, false)));
            });
            while (!f.isDone()) {
            }
            if (resultsMap.get(option) > 20.0) {
                break;
            }
        }
//        
        for (String option : resultsMap.keySet()) {
            if (resultsMap.get(option) > idRate) {
                idRate = resultsMap.get(option);
                selectedOption = option;
            }
        }
        if (idRate > referenceIdRate) {
            referenceIdRate = idRate;
            System.out.println("------------->> I work is done " + resultsMap);
            System.out.println("param updated " + selectedOption);
            optimisedSearchParameter.setDigestionParameter(selectedOption);
        }
        long end1 = System.currentTimeMillis();
        double total = (end1 - start1) / 1000.0;
        ReportExporter.addElementToReport(Configurations.getDataset_Id(), "Cleavage_Parameters", optimisedSearchParameter.getDigestionParameter(), resultsMap.get(optimisedSearchParameter.getDigestionParameter()), total);

    }

    private void optimizeEnzymeParameter() throws IOException {
        // Processing        
        long start1 = System.currentTimeMillis();
        Map<String, Double> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        double idRate = -1;
        String selectedOption = "";
        //optimise enzyme  
        Enzyme[] enzymes = new Enzyme[EnzymeFactory.getInstance().getEnzymes().size()];
        enzymes[0] = EnzymeFactory.getInstance().getEnzyme("Trypsin");
        int i = 1;
        for (Enzyme enzyme : EnzymeFactory.getInstance().getEnzymes()) {
            if (enzyme.getName().equals("Trypsin")) {
                continue;
            }
            enzymes[i] = enzyme;
            i++;
        }
        for (Enzyme enzyme : enzymes) {
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            tempIdParam.getSearchParameters().getDigestionParameters().clearEnzymes();
            tempIdParam.getSearchParameters().getDigestionParameters().addEnzyme(enzyme);
            final String option = enzyme.getName();
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option;
            Future f = MainUtilities.getExecutorService().submit(() -> {
                resultsMap.put(option, countSpectra(excuteSearch(updatedName, "Enzyme_Parameters", option, tempIdParam, false)));
            });
            while (!f.isDone()) {
            }
            if (enzyme.getName().equals("Trypsin") && resultsMap.get(option) > 20.0) {
                break;
            }

        }

        for (String option : resultsMap.keySet()) {
            if (resultsMap.get(option) > idRate) {
                idRate = resultsMap.get(option);
                selectedOption = option;
            }
        }
        if (idRate > referenceIdRate) {
            referenceIdRate = idRate;
            System.out.println("------------->> II work is done " + resultsMap);
            System.out.println("param updated " + selectedOption);
            optimisedSearchParameter.setEnzymeName(selectedOption);
        }
        long end1 = System.currentTimeMillis();
        double total = (end1 - start1) / 1000.0;

        ReportExporter.addElementToReport(Configurations.getDataset_Id(), "Enzyme_Parameters", optimisedSearchParameter.getEnzymeName(), resultsMap.get(optimisedSearchParameter.getEnzymeName()), total);

    }

    private void optimizeSpecificityParameter() throws IOException {
        Map<String, Double> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        double idRate = -1;
        String selectedOption = "";
        for (int i = 0; i < DigestionParameters.Specificity.values().length; i++) {
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            tempIdParam.getSearchParameters().getDigestionParameters().setSpecificity(optimisedSearchParameter.getEnzymeName(), DigestionParameters.Specificity.getSpecificity(i));
            final String option = DigestionParameters.Specificity.getSpecificity(i).name();
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option;

            Future f = MainUtilities.getExecutorService().submit(() -> {
                resultsMap.put(option, countSpectra(excuteSearch(updatedName, "Specific_Parameters", option, tempIdParam, false)));
            });
            while (!f.isDone()) {
            }
        }
        for (String option : resultsMap.keySet()) {
            if (resultsMap.get(option) > idRate) {
                idRate = resultsMap.get(option);
                selectedOption = option;
            }
        }
        if (idRate > referenceIdRate) {
            referenceIdRate = idRate;
            System.out.println("------------->> III work is done " + resultsMap);
//            System.out.println("param updated " + selectedOption);
            optimisedSearchParameter.setEnzymeSpecificity(selectedOption);
        }

    }

    private void optimizeMaxMissCleavagesParameter() throws IOException {
        Map<String, Double> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        double idRate = -1;
        String selectedOption = "";
        for (int i = 0; i < 7; i++) {
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            tempIdParam.getSearchParameters().getDigestionParameters().setnMissedCleavages(optimisedSearchParameter.getEnzymeName(), i);
            final String option = "missedCleavages_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option;
            
            Future future = MainUtilities.getExecutorService().submit(() -> {
                resultsMap.put(option, countSpectra(excuteSearch(updatedName, "missCleav_Parameters", option, tempIdParam, false)));
            });
            while (!future.isDone()) {
            }
        }
        for (String option : resultsMap.keySet()) {
            if (resultsMap.get(option) > idRate) {
                idRate = resultsMap.get(option);
                selectedOption = option;

            }
        }

        if (idRate > referenceIdRate) {
            int selectedNum = Integer.parseInt(selectedOption.split("_")[1]);
            referenceIdRate = idRate;
            System.out.println("------------->> IV work is done " + resultsMap + "  " + selectedNum);
            System.out.println("param updated " + selectedOption);
            optimisedSearchParameter.setMaxMissedCleavage(selectedNum);
        }

    }

    private void optimizeFragmentIonTypesParameter() throws IOException {
        Map<String, Double> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        double idRate = -1;
        IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        SearchParameters tempSearchParameters = tempIdParam.getSearchParameters();
        ArrayList<Integer> selectedForwardIons = tempSearchParameters.getForwardIons();
        String[] forwardIons = new String[]{"a", "b", "c"};
        String[] rewindIons = new String[]{"x", "y", "z"};

        for (String forwardIon : forwardIons) {
            selectedForwardIons.clear();
            Integer forwardIonType = PeptideFragmentIon.getIonType(forwardIon);
            selectedForwardIons.add(forwardIonType);
            for (String rewindIon : rewindIons) {
                Integer rewindIonType = PeptideFragmentIon.getIonType(rewindIon);
                ArrayList<Integer> selectedRewindIons = new ArrayList<>();
                selectedRewindIons.add(rewindIonType);
                tempSearchParameters.setRewindIons(selectedRewindIons);
                String option = selectedForwardIons + "-" + selectedRewindIons;
                final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option;

                Future f = MainUtilities.getExecutorService().submit(() -> {
                    resultsMap.put(option, countSpectra(excuteSearch(updatedName, "Ions_Parameters", option, tempIdParam, false)));
                });
                while (!f.isDone()) {
                }
            }
        }

        String selectedOption = "";
        for (String option : resultsMap.keySet()) {
            if (resultsMap.get(option) > idRate) {
                idRate = resultsMap.get(option);
                selectedOption = option;
            }
        }
        selectedOption = selectedOption.replace("[", "").replace("]", "");
        int forward = Integer.parseInt(selectedOption.split("-")[0]);
        int rewind = Integer.parseInt(selectedOption.split("-")[1]);
        ArrayList<Integer> rewindIonsList = new ArrayList<>();
        ArrayList<Integer> forwardIonsList = new ArrayList<>();
        forwardIonsList.add(forward);
        rewindIonsList.add(rewind);
        if (idRate > referenceIdRate) {
            referenceIdRate = idRate;
            System.out.println("------------->> V work is done " + forwardIonsList + "  " + rewindIonsList);
            System.out.println("param updated " + selectedOption);
            optimisedSearchParameter.setSelectedForwardIons(forwardIonsList);
            optimisedSearchParameter.setSelectedRewindIons(rewindIonsList);
        }

    }

    private void optimizePrecursorToleranceParameter() throws IOException {
        Map<String, Double> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        double idRate = -1;
        String selectedOption = "";
        for (double i = 5; i < 30;) {
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            tempIdParam.getSearchParameters().setPrecursorAccuracy(i);
            tempIdParam.getSearchParameters().setPrecursorAccuracyType(SearchParameters.MassAccuracyType.PPM);
            final String option = "precursorAccuracy_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option;

            Future future = MainUtilities.getExecutorService().submit(() -> {
                resultsMap.put(option, countSpectra(excuteSearch(updatedName, "precursorAccuracy_", option, tempIdParam, false)));
            });
            while (!future.isDone()) {
            }
            i += 5;
        }
        for (String option : resultsMap.keySet()) {
            if (resultsMap.get(option) > idRate) {
                idRate = resultsMap.get(option);
                selectedOption = option;

            }
        }
        if (idRate > referenceIdRate) {
            referenceIdRate = idRate;
            System.out.println("------------->> VI work is done " + resultsMap + "  " + selectedOption);
            System.out.println("param updated " + selectedOption);
            double selectedNum = Double.parseDouble(selectedOption.split("_")[1]);
            optimisedSearchParameter.setPrecursorTolerance(selectedNum);
        }

    }

    private void optimizeFragmentToleranceParameter() throws IOException {
        Map<String, Double> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        double idRate = -1;
        String selectedOption = "";
        double[] values = new double[]{0.02, 0.05, 0.07, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5};
        for (double i : values) {
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            tempIdParam.getSearchParameters().setFragmentIonAccuracy(i);
            tempIdParam.getSearchParameters().setFragmentAccuracyType(SearchParameters.MassAccuracyType.DA);
            final String option = "fragmentAccuracy_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option;

            Future future = MainUtilities.getExecutorService().submit(() -> {
                resultsMap.put(option, countSpectra(excuteSearch(updatedName, "fragment_accuracy", option, tempIdParam, false)));
            });
            while (!future.isDone()) {
            }
        }
        for (String option : resultsMap.keySet()) {
            if (resultsMap.get(option) > idRate) {
                idRate = resultsMap.get(option);
                selectedOption = option;

            }
        }

        if (idRate > referenceIdRate) {
            double selectedNum = Double.parseDouble(selectedOption.split("_")[1]);
            referenceIdRate = idRate;
            System.out.println("------------->> VII work is done " + resultsMap + "  " + selectedNum);
            System.out.println("param updated " + selectedOption);
            optimisedSearchParameter.setFragmentTolerance(selectedNum);
        }

    }

    private void optimizePrecursorChargeParameter() throws IOException {
        Map<String, Double> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        double idRate = -1;
        String selectedOption = "";
        int selectedNum = 100;
        for (int i = 2; i < 6; i++) {
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            tempIdParam.getSearchParameters().setMaxChargeSearched(i);
//            tempIdParam.getSearchParameters().setMinChargeSearched(0);
            final String option = "maxCharge_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option;
            Future future = MainUtilities.getExecutorService().submit(() -> {
                resultsMap.put(option, countSpectra(excuteSearch(updatedName, "charge_Parameters", option, tempIdParam, false)));
            });
            while (!future.isDone()) {
            }
        }
        for (String option : resultsMap.keySet()) {
            if (resultsMap.get(option) > idRate) {
                idRate = resultsMap.get(option);
                selectedOption = option;

            }
        }
        selectedNum = Integer.parseInt(selectedOption.split("_")[1]);
        if (idRate > referenceIdRate) {
            referenceIdRate = idRate;
            System.out.println("------------->> VIII work is done " + resultsMap + "  " + selectedNum);
            System.out.println("param updated " + selectedOption);
            optimisedSearchParameter.setMaxPrecursorCharge(selectedNum);
        }

        resultsMap.clear();
        for (int i = optimisedSearchParameter.getMaxPrecursorCharge() - 1; i > 0; i--) {
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//            tempIdParam.getSearchParameters().setMaxChargeSearched(optimisedSearchParameter.getMaxPrecursorCharge());
            tempIdParam.getSearchParameters().setMinChargeSearched(i);
            final String option = "minCharge_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option;
            Future future = MainUtilities.getExecutorService().submit(() -> {
                resultsMap.put(option, countSpectra(excuteSearch(updatedName, "charge_Parameters", option, tempIdParam, false)));
            });
            while (!future.isDone()) {
            }
        }

        idRate = -1;
        for (String option : resultsMap.keySet()) {
            int valueNum = Integer.parseInt(option.split("_")[1]);
            if (resultsMap.get(option) > idRate) {
                idRate = resultsMap.get(option);
                selectedOption = option;
                selectedNum = valueNum;
            } else if (resultsMap.get(option) == idRate && valueNum < selectedNum) {
                selectedNum = valueNum;
                selectedOption = option;
            }
        }

        if (idRate > referenceIdRate) {
            selectedNum = Integer.parseInt(selectedOption.split("_")[1]);
            referenceIdRate = idRate;
            System.out.println("------------->> VIIII work is done " + resultsMap + "  " + selectedNum);
            System.out.println("param updated " + selectedOption);
            optimisedSearchParameter.setMinPrecursorCharge(selectedNum);
        }

    }

    private void optimizeIsotopParameter() throws IOException {
        Map<String, Double> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        double idRate = -1;
        String selectedOption = "";
        for (int i = 7; i > -8; i--) {
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//            tempIdParam.getSearchParameters().setMaxIsotopicCorrection(8);
            tempIdParam.getSearchParameters().setMinChargeSearched(i);
            final String option = "minIsotop_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option;

            Future future = MainUtilities.getExecutorService().submit(() -> {
                resultsMap.put(option, countSpectra(excuteSearch(updatedName, "isotop_Parameters", option, tempIdParam, false)));
            });
            while (!future.isDone()) {
            }
        }
        for (String option : resultsMap.keySet()) {
            if (resultsMap.get(option) > idRate) {
                idRate = resultsMap.get(option);
                selectedOption = option;

            }
        }

        int selectedNum = 100;
        for (String option : resultsMap.keySet()) {
            int valueNum = Math.abs(Integer.parseInt(option.split("_")[1]));
            if (resultsMap.get(option) > idRate) {
                idRate = resultsMap.get(option);
                selectedOption = option;
                selectedNum = valueNum;
            } else if (resultsMap.get(option) == idRate && valueNum < selectedNum) {
                selectedNum = valueNum;
                selectedOption = option;
            }
        }

        if (idRate > referenceIdRate) {
            selectedNum = Integer.parseInt(selectedOption.split("_")[1]);
            referenceIdRate = idRate;
            System.out.println("------------->> X work is done " + resultsMap + "  " + selectedNum);
            System.out.println("param updated " + selectedOption);
            optimisedSearchParameter.setMinIsotops(selectedNum);
        }
        resultsMap.clear();
        idRate = -1;
        for (int i = optimisedSearchParameter.getMinIsotops() + 1; i < 9; i++) {
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            tempIdParam.getSearchParameters().setMaxIsotopicCorrection(i);
//            tempIdParam.getSearchParameters().setMinIsotopicCorrection(optimisedSearchParameter.getMinIsotops());
            final String option = "maxIsotop_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option;

            Future future = MainUtilities.getExecutorService().submit(() -> {
                resultsMap.put(option, countSpectra(excuteSearch(updatedName, "isotop_Parameters", option, tempIdParam, false)));
            });
            while (!future.isDone()) {
            }
        }
        for (String option : resultsMap.keySet()) {
            if (resultsMap.get(option) > idRate) {
                idRate = resultsMap.get(option);
                selectedOption = option;

            }
        }
        selectedNum = Integer.parseInt(selectedOption.split("_")[1]);
        if (idRate > referenceIdRate) {
            referenceIdRate = idRate;
            System.out.println("------------->> XI work is done " + resultsMap + "  " + selectedNum);
            System.out.println("param updated " + selectedOption);
            optimisedSearchParameter.setMaxIsotops(selectedNum);
        }

    }
    int refScore = 0;

    private void optimizeModificationsParameter() throws IOException {
        List<String> mods = new ArrayList<>();
        mods.addAll(ptmFactory.getModifications(ModificationCategory.Common));
        mods.addAll(ptmFactory.getModifications(ModificationCategory.Common_Biological));
        mods.addAll(ptmFactory.getModifications(ModificationCategory.Common_Artifact));
        Map<Integer, Set<String>> orderedFixedModificationScore = new TreeMap<>(Collections.reverseOrder());
        Map<String, Integer> fixedModificationScore = new LinkedHashMap<>();
        Set<String> defaultFixedModifications = new HashSet<>();
        String ref_mod = "reference_mod";
        IdentificationParameters idParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        idParam.getSearchParameters().getModificationParameters().clearFixedModifications();
        idParam.getSearchParameters().getModificationParameters().clearVariableModifications();
        idParam.getSearchParameters().getModificationParameters().clearRefinementModifications();
        idParam.getSearchParameters().getModificationParameters().getRefinementFixedModifications().clear();
        final String referenc_mod_Name = Configurations.DEFAULT_RESULT_NAME + "_" + ref_mod;
        ArrayList<SpectrumMatch> sms = testModifications(referenc_mod_Name, ref_mod, defaultFixedModifications, new HashSet<>(), false);
        refScore = sms.size();

        for (String modId : mods) {
            final String option = modId;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "f_" + option;
            Set<String> fModifications = new HashSet<>();
            fModifications.add(modId);
            sms = testModifications(updatedName, option, defaultFixedModifications, fModifications, false);
            if (sms.size() > refScore) {
                if (!orderedFixedModificationScore.containsKey(sms.size())) {
                    orderedFixedModificationScore.put(sms.size(), new LinkedHashSet<>());
                }
                orderedFixedModificationScore.get(sms.size()).add(modId);
            }
        }

        //select only modification with unique target 
        Set<String> targetedAAToModList = new LinkedHashSet<>();
        for (int score : orderedFixedModificationScore.keySet()) {
            Set<String> fmods = orderedFixedModificationScore.get(score);
            for (String mod : fmods) {
                if (!targetedAAToModList.contains(ptmFactory.getModification(mod).getPattern().toString())) {
                    fixedModificationScore.put(mod, score);
                    targetedAAToModList.add(ptmFactory.getModification(mod).getPattern().toString());
                }

            }
        }
        MainUtilities.cleanOutputFolder();
        final String tempOption = fixedModificationScore.keySet().toString().replace(",", "_").replace("[", "").replace("]", "");
        final String tempUpdatedName = Configurations.DEFAULT_RESULT_NAME + "f_" + tempOption;
        Set<String> fModifications = new HashSet<>(fixedModificationScore.keySet());
        sms = testModifications(tempUpdatedName, tempOption, defaultFixedModifications, fModifications, false);
        refScore = sms.size();
        defaultFixedModifications.addAll(fixedModificationScore.keySet());
        ///check refined modification 
        for (String fixedMod : fixedModificationScore.keySet()) {
            idParam.getSearchParameters().getModificationParameters().addRefinementFixedModification(ptmFactory.getModification(fixedMod));
        }
        optimisedSearchParameter.getFixedModifications().clear();
        optimisedSearchParameter.getFixedModifications().addAll(fixedModificationScore.keySet());
        optimisedSearchParameter.setRefinedFixedModifications(fixedModificationScore.keySet());

        //start variable modification check
        TreeMap<Integer, ArrayList<String>> topVScoreMap = new TreeMap<>(Collections.reverseOrder());
        topVScoreMap.put(refScore, new ArrayList<>());
        for (String modId : mods) {
            if (defaultFixedModifications.contains(modId)) {
                continue;
            }
            final String option = modId;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "v_" + option;
            Set<String> vModifications = new HashSet<>();
            vModifications.add(modId);
            sms = testModifications(updatedName, option, defaultFixedModifications, vModifications, true);

            if (sms.size() > refScore) {
                if (!topVScoreMap.containsKey(sms.size())) {
                    topVScoreMap.put(sms.size(), new ArrayList<>());
                }
                topVScoreMap.get(sms.size()).add(modId);
                System.out.println("at variable v1 mod " + option + "  " + sms.size());
            }
        }
        Map<String, Integer> allFilteredVariableModifications = new LinkedHashMap<>();
        for (int key : topVScoreMap.keySet()) {
            ArrayList<String> vmods = topVScoreMap.get(key);
            for (String mod : vmods) {
                allFilteredVariableModifications.put(mod, key);
                if (allFilteredVariableModifications.size() >= 5) {
                    break;
                }
                if (allFilteredVariableModifications.size() >= 5) {
                    break;
                }
            }

        }
        int topScore = topVScoreMap.firstKey();
        topVScoreMap.clear();
        topVScoreMap.put(topScore, new ArrayList<>());
        int countI = 0;
        for (String modI : allFilteredVariableModifications.keySet()) {
            int countII = 0;
            for (String modII : allFilteredVariableModifications.keySet()) {
                if (countII <= countI) {
                    countII++;
                    continue;
                }
                final String option = modI + "_" + modII;
                final String updatedName = Configurations.DEFAULT_RESULT_NAME + "v2_" + option;
                Set<String> vModifications = new HashSet<>();
                vModifications.add(modI);
                vModifications.add(modII);
                sms = testModifications(updatedName, option, defaultFixedModifications, vModifications, true);
                System.out.println("at variable v2  mod " + option + "  " + sms.size());
                if (topVScoreMap.firstKey() < sms.size()) {
                    topVScoreMap.put(sms.size(), new ArrayList<>());
                    topVScoreMap.get(sms.size()).add(modI);
                    topVScoreMap.get(sms.size()).add(modII);
                }
                countII++;
            }
            countI++;
        }

        TreeMap<Integer, ArrayList<String>> finalVariableMdificationsMap = new TreeMap<>(Collections.reverseOrder());
        finalVariableMdificationsMap.putAll(topVScoreMap);
        for (String modI : allFilteredVariableModifications.keySet()) {
            for (ArrayList<String> modIIs : topVScoreMap.values()) {
                if (modIIs.contains(modI)) {
                    continue;
                }
                final String option = modI + "_" + modIIs.toString();
                final String updatedName = Configurations.DEFAULT_RESULT_NAME + "v3_" + option;
                Set<String> vModifications = new HashSet<>(modIIs);
                vModifications.add(modI);
                sms = testModifications(updatedName, option, defaultFixedModifications, vModifications, true);
                if (finalVariableMdificationsMap.firstKey() < sms.size()) {
                    finalVariableMdificationsMap.put(sms.size(), new ArrayList<>(modIIs));
                    finalVariableMdificationsMap.get(sms.size()).add(modI);
                    System.out.println("at variable v3 mod " + option + "  " + sms.size());
                }
            }
        }

        optimisedSearchParameter.setSortedVariableModificationsMap(finalVariableMdificationsMap);
        optimisedSearchParameter.setRefinedVariableModifications(allFilteredVariableModifications.keySet());
        for (int i : finalVariableMdificationsMap.keySet()) {
            System.out.println("at vm " + i + " " + finalVariableMdificationsMap.get(i));
        }
        for (String i : allFilteredVariableModifications.keySet()) {
            System.out.println("at vm " + allFilteredVariableModifications.get(i) + " " + i);
        }

//        System.exit(0);
    }

    private ArrayList<SpectrumMatch> testModifications(String updatedName, String option, Set<String> defaultFixedModifications, Set<String> modifications, boolean variableModificatins) throws IOException {

        ArrayList<SpectrumMatch> spectrumMaches = new ArrayList<>();
        IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        tempIdParam.getSearchParameters().getModificationParameters().clearFixedModifications();
        tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
//        tempIdParam.getSearchParameters().getModificationParameters().getRefinementFixedModifications().clear();
//        tempIdParam.getSearchParameters().getModificationParameters().getRefinementVariableModifications().clear();
        if (optProtSearchParameters.isRunXTandem()) {
            boolean onlyTerminal = !variableModificatins;
            for (String mod : modifications) {
                if (!ptmFactory.getModification(mod).getModificationType().isNTerm() && !ptmFactory.getModification(mod).getModificationType().isCTerm()) {
                    onlyTerminal = false;
                    break;
                }
            }

            XtandemParameters xtandemParameters = (XtandemParameters) tempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
            if (xtandemParameters.isQuickPyrolidone() && xtandemParameters.isProteinQuickAcetyl() && onlyTerminal && !modifications.isEmpty()) {
                return new ArrayList<>();

            }
        }
        for (String modId : defaultFixedModifications) {
            tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(modId));
        }
        for (String modId : modifications) {
            if (variableModificatins) {
                tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(modId));
            } else {
                tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(modId));
            }
        }
        Future future = MainUtilities.getExecutorService().submit(() -> {
            ArrayList<SpectrumMatch> sms = excuteSearch(updatedName, "ref_modification__Parameters", option, tempIdParam, false);
            spectrumMaches.addAll(sms);

        });
        while (!future.isDone()) {
        }
        return spectrumMaches;

    }

    private double countSpectra(ArrayList<SpectrumMatch> matches) {
        try {
            if (matches.isEmpty()) {
                return -1;
            }
            MsFileHandler msFileHandler = new MsFileHandler();
            msFileHandler.register(subMsFile, MainUtilities.OptProt_Waiting_Handler);
            double leng = msFileHandler.getSpectrumTitles(IoUtil.removeExtension(subMsFile.getName())).length;

            return MainUtilities.rundDouble((matches.size() * 100.0) / leng);
        } catch (IOException ex) {
            Logger.getLogger(SearchOptimizerHandler.class.getName()).log(Level.SEVERE, null, ex);
        }
        return -1;
    }

    public ArrayList<SpectrumMatch> excuteSearch(String defaultOutputFileName, String paramName, String paramOption, IdentificationParameters tempIdParam, boolean addPeptideMass) {
        if (optProtSearchParameters.isRunXTandem()) {
            return excuteXTandom(defaultOutputFileName, paramOption, tempIdParam, addPeptideMass);
        } else if (optProtSearchParameters.isRunMyriMatch()) {
            return excuteMyriMatch(defaultOutputFileName, paramOption, tempIdParam, addPeptideMass);
        }
        return new ArrayList<>();
    }

    public synchronized Set<String> excuteNovor(SearchInputSetting searchEngineParameters) {

        try {
            String spectraFileName = IoUtil.removeExtension(subMsFile.getName());
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + spectraFileName;
            this.optProtSearchParameters = searchEngineParameters;
            File configFolder = new File(resultsOutput, updatedName + "_temp");
            System.out.println("at output for novor " + configFolder.getAbsolutePath());
            configFolder.mkdir();
            File resultOutput = new File(resultsOutput, updatedName);
            resultOutput.mkdir();
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            SearchParameters searchParameters = tempIdParam.getSearchParameters();
            String error = SearchHandler.loadModifications(searchParameters);
            if (error != null) {
                System.out.println(error);
            }
            OptProtWaitingHandler waitingHandlerCLIImpl = new OptProtWaitingHandler();
            searchEngineParameters.setRunNovor(true);
            SearchHandler searchHandler = new SearchHandler(tempIdParam, resultOutput, configFolder, spectraFileName, msFileInList,
                    subFastaFile, new ArrayList<>(),
                    identificationParametersFile,
                    searchEngineParameters.isRunOmssa(),
                    searchEngineParameters.isRunXTandem(),
                    searchEngineParameters.isRunMsgf(),
                    searchEngineParameters.isRunMsAmanda(),
                    this.optProtSearchParameters.isRunMyriMatch(),
                    this.optProtSearchParameters.isRunComet(),
                    this.optProtSearchParameters.isRunTide(),
                    this.optProtSearchParameters.isRunAndromeda(),
                    this.optProtSearchParameters.isRunMetaMorpheus(),
                    this.optProtSearchParameters.isRunSage(),
                    this.optProtSearchParameters.isRunNovor(),
                    this.optProtSearchParameters.isRunDirecTag(),
                    this.optProtSearchParameters.getOmssaFolder(),
                    this.optProtSearchParameters.getxTandemFolder(),
                    this.optProtSearchParameters.getMsgfFolder(),
                    this.optProtSearchParameters.getMsAmandaFolder(),
                    this.optProtSearchParameters.getMyriMatchFolder(),
                    this.optProtSearchParameters.getCometFolder(),
                    this.optProtSearchParameters.getTideFolder(),
                    this.optProtSearchParameters.getTideIndexLocation(),
                    this.optProtSearchParameters.getAndromedaFolder(),
                    this.optProtSearchParameters.getMetaMorpheusFolder(),
                    this.optProtSearchParameters.getSageFolder(),
                    this.optProtSearchParameters.getNovorFolder(),
                    this.optProtSearchParameters.getDirecTagFolder(),
                    this.optProtSearchParameters.getMakeblastdbFolder(),
                    null
            );

            searchHandler.startSearch(waitingHandlerCLIImpl);
            File resultsFile = searchHandler.getResultsFolder();
            File NovorFile = new File(resultsFile, spectraFileName + ".novor.csv");
            Set<String> sequences = SpectraUtilities.getSequences(NovorFile);
            return sequences;
        } catch (IOException | InterruptedException ex) {
            ex.printStackTrace();
        }
        return new HashSet<>();
    }

    public synchronized Map<String, String> excuteDirecTag(SearchInputSetting searchEngineParameters) {

        try {
            String spectraFileName = IoUtil.removeExtension(subMsFile.getName());
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + spectraFileName;
            this.optProtSearchParameters = searchEngineParameters;
            File configFolder = new File(resultsOutput, updatedName + "_temp");
            configFolder.mkdir();
            File resultOutput = new File(resultsOutput, updatedName);
            resultOutput.mkdir();
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            SearchParameters searchParameters = tempIdParam.getSearchParameters();
            DirecTagParameters direcTagParameters = (DirecTagParameters) searchParameters.getIdentificationAlgorithmParameter(Advocate.direcTag.getIndex());
            direcTagParameters.setMaxTagCount(1);
            direcTagParameters.setTagLength(4);
            direcTagParameters.setNumChargeStates(4);
            OptProtWaitingHandler waitingHandlerCLIImpl = new OptProtWaitingHandler();
            this.optProtSearchParameters.setRunXTandem(false);
            this.optProtSearchParameters.setRunDirecTag(true);
            SearchHandler searchHandler = new SearchHandler(tempIdParam, resultOutput, configFolder, spectraFileName, msFileInList,
                    subFastaFile, new ArrayList<>(),
                    identificationParametersFile,
                    this.optProtSearchParameters.isRunOmssa(),
                    this.optProtSearchParameters.isRunXTandem(),
                    this.optProtSearchParameters.isRunMsgf(),
                    this.optProtSearchParameters.isRunMsAmanda(),
                    this.optProtSearchParameters.isRunMyriMatch(),
                    this.optProtSearchParameters.isRunComet(),
                    this.optProtSearchParameters.isRunTide(),
                    this.optProtSearchParameters.isRunAndromeda(),
                    this.optProtSearchParameters.isRunMetaMorpheus(),
                    this.optProtSearchParameters.isRunSage(),
                    this.optProtSearchParameters.isRunNovor(),
                    this.optProtSearchParameters.isRunDirecTag(),
                    this.optProtSearchParameters.getOmssaFolder(),
                    this.optProtSearchParameters.getxTandemFolder(),
                    this.optProtSearchParameters.getMsgfFolder(),
                    this.optProtSearchParameters.getMsAmandaFolder(),
                    this.optProtSearchParameters.getMyriMatchFolder(),
                    this.optProtSearchParameters.getCometFolder(),
                    this.optProtSearchParameters.getTideFolder(),
                    this.optProtSearchParameters.getTideIndexLocation(),
                    this.optProtSearchParameters.getAndromedaFolder(),
                    this.optProtSearchParameters.getMetaMorpheusFolder(),
                    this.optProtSearchParameters.getSageFolder(),
                    this.optProtSearchParameters.getNovorFolder(),
                    this.optProtSearchParameters.getDirecTagFolder(),
                    this.optProtSearchParameters.getMakeblastdbFolder(),
                    null
            );
//            MainUtilities.getExecutorService() = MainUtilities.getExecutorService()s.newFixedThreadPool(2);
            Future future = MainUtilities.getExecutorService().submit(() -> {
                try {
                    searchHandler.startSearch(waitingHandlerCLIImpl);
                } catch (InterruptedException ex) {
                    Logger.getLogger(SearchOptimizerHandler.class.getName()).log(Level.SEVERE, null, ex);
                }
            });
            while (!future.isDone()) {
            }
//            MainUtilities.getExecutorService().shutdown();
            File resultsFile = searchHandler.getResultsFolder();
            File NovorFile = new File(resultsFile, spectraFileName + ".tags");
            IdfileReader idReader = readerFactory.getFileReader(NovorFile);
            MsFileHandler msFileHandler = new MsFileHandler();
            msFileHandler.register(subMsFile, waitingHandlerCLIImpl);
            try {
                ArrayList<SpectrumMatch> matches = idReader.getAllSpectrumMatches(msFileHandler, waitingHandlerCLIImpl, searchParameters);
                Map<String, String> specTagMap = new LinkedHashMap<>();
                for (SpectrumMatch sm : matches) {
                    TagAssumption tag = sm.getAllTagAssumptions().toList().get(0);
                    if (tag.getScore() <= 0.01) {
                        specTagMap.put(sm.getSpectrumTitle(), tag.getTag().getContent().get(1).asSequence());
                    }
                }
                return specTagMap;
            } catch (SQLException | ClassNotFoundException | JAXBException | XmlPullParserException | XMLStreamException ex) {
                Logger.getLogger(SearchOptimizerHandler.class.getName()).log(Level.SEVERE, null, ex);
            }

        } catch (IOException | InterruptedException ex) {
            ex.printStackTrace();
        }
        return new LinkedHashMap<>();
    }

    private synchronized ArrayList<SpectrumMatch> excuteXTandom(String defaultOutputFileName, String paramOption, IdentificationParameters tempIdParam, boolean addPeptideMasses) {

        try {

            if (!optProtSearchParameters.getXTandemEnabledParameters().getParamsToOptimize().isEnabledParam(paramOption.split("_")[0])) {
                System.out.println("param " + paramOption + " is not supported " + paramOption);
                return new ArrayList<>();
            }
            SearchParameters searchParameters = tempIdParam.getSearchParameters();
            XtandemParameters xtandemParameters = (XtandemParameters) searchParameters.getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
            xtandemParameters.setQuickPyrolidone(true);
            xtandemParameters.setProteinQuickAcetyl(true);
            xtandemParameters.setStpBias(true);
            xtandemParameters.setRefine(true);

            File resultOutput = SearchExecuter.executeSearch(defaultOutputFileName, optProtSearchParameters, subMsFile, subFastaFile, tempIdParam, optimizedIdentificationParametersFile);
            File xTandemFile = new File(resultOutput, SearchHandler.getXTandemFileName(IoUtil.removeExtension(subMsFile.getName())));
            IdfileReader idReader = readerFactory.getFileReader(xTandemFile);
            MsFileHandler msFileHandler = new MsFileHandler();
            msFileHandler.register(subMsFile, MainUtilities.OptProt_Waiting_Handler);
            final List<SpectrumMatch> matches = Collections.synchronizedList(idReader.getAllSpectrumMatches(msFileHandler, MainUtilities.OptProt_Waiting_Handler, searchParameters));

            if (addPeptideMasses) {
                SequenceMatchingParameters modificationSequenceMatchingParameters = tempIdParam.getModificationLocalizationParameters().getSequenceMatchingParameters();
                FMIndex sequenceProvider = new FMIndex(subFastaFile, null, new OptProtWaitingHandler(), false, tempIdParam);
//                MainUtilities.getExecutorService() = MainUtilities.getExecutorService()s.newFixedThreadPool(2);
                Future future = MainUtilities.getExecutorService().submit(() -> {
                    for (SpectrumMatch sm : matches) {
                        for (PeptideAssumption pepAss : sm.getAllPeptideAssumptions().toList()) {
                            Peptide pep = pepAss.getPeptide();
                            ModificationLocalizationMapper.modificationLocalization(
                                    pep,
                                    tempIdParam,
                                    idReader,
                                    ptmFactory,
                                    sequenceProvider
                            );
                            pepAss.getPeptide().setMass(pep.getMass(searchParameters.getModificationParameters(), sequenceProvider, modificationSequenceMatchingParameters));
                            sm.setBestPeptideAssumption(pepAss);
                        }
                    }
                });
                while (!future.isDone()) {
                }
//                MainUtilities.getExecutorService().shutdown();
            }
            return new ArrayList<>(matches);
        } catch (IOException | InterruptedException | SQLException | ClassNotFoundException | JAXBException | XMLStreamException | XmlPullParserException ex) {
            ex.printStackTrace();
        }
        return new ArrayList<>();
    }

    private synchronized ArrayList<SpectrumMatch> excuteMyriMatch(String defaultOutputFileName, String paramOption, IdentificationParameters tempIdParam, boolean addPeptideMasses) {

        try {

            if (!optProtSearchParameters.getMyriMatchEnabledParameters().getParamsToOptimize().isEnabledParam(paramOption.split("_")[0])) {
                System.out.println("param " + paramOption + " is not supported " + paramOption);
                return new ArrayList<>();
            }
            SearchParameters searchParameters = tempIdParam.getSearchParameters();
            MyriMatchParameters myriMatchParameters = (MyriMatchParameters) searchParameters.getAlgorithmSpecificParameters().get(Advocate.myriMatch.getIndex());
            File resultOutput = SearchExecuter.executeSearch(defaultOutputFileName, optProtSearchParameters, subMsFile, subFastaFile, tempIdParam, optimizedIdentificationParametersFile);
            File myriMatchFile = new File(resultOutput, SearchHandler.getMyriMatchFileName(IoUtil.removeExtension(subMsFile.getName()), myriMatchParameters));
            IdfileReader idReader = readerFactory.getFileReader(myriMatchFile);
            MsFileHandler msFileHandler = new MsFileHandler();
            msFileHandler.register(subMsFile, MainUtilities.OptProt_Waiting_Handler);
            final List<SpectrumMatch> matches = Collections.synchronizedList(idReader.getAllSpectrumMatches(msFileHandler, MainUtilities.OptProt_Waiting_Handler, searchParameters));

            if (addPeptideMasses) {
                SequenceMatchingParameters modificationSequenceMatchingParameters = tempIdParam.getModificationLocalizationParameters().getSequenceMatchingParameters();
                FMIndex sequenceProvider = new FMIndex(subFastaFile, null, new OptProtWaitingHandler(), false, tempIdParam);
//                MainUtilities.getExecutorService() = MainUtilities.getExecutorService()s.newFixedThreadPool(2);
                Future future = MainUtilities.getExecutorService().submit(() -> {
                    for (SpectrumMatch sm : matches) {
                        for (PeptideAssumption pepAss : sm.getAllPeptideAssumptions().toList()) {
                            Peptide pep = pepAss.getPeptide();
                            ModificationLocalizationMapper.modificationLocalization(
                                    pep,
                                    tempIdParam,
                                    idReader,
                                    ptmFactory,
                                    sequenceProvider
                            );
                            pepAss.getPeptide().setMass(pep.getMass(searchParameters.getModificationParameters(), sequenceProvider, modificationSequenceMatchingParameters));
                            sm.setBestPeptideAssumption(pepAss);
                        }
                    }
                });
                while (!future.isDone()) {
                }
//                MainUtilities.getExecutorService().shutdown();
            }
            return new ArrayList<>(matches);
        } catch (IOException | InterruptedException | SQLException | ClassNotFoundException | JAXBException | XMLStreamException | XmlPullParserException ex) {
            ex.printStackTrace();
        }
        return new ArrayList<>();
    }

}
