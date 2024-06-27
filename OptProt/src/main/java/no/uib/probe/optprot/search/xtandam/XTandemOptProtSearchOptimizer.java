/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package no.uib.probe.optprot.search.xtandam;

import com.compomics.util.experiment.biology.enzymes.EnzymeFactory;
import com.compomics.util.experiment.biology.modifications.ModificationFactory;
import com.compomics.util.experiment.identification.Advocate;
import com.compomics.util.experiment.identification.matches.SpectrumMatch;
import com.compomics.util.experiment.io.identification.IdfileReaderFactory;
import com.compomics.util.parameters.identification.IdentificationParameters;
import com.compomics.util.parameters.identification.search.DigestionParameters;
import com.compomics.util.parameters.identification.search.SearchParameters;
import com.compomics.util.parameters.identification.tool_specific.XtandemParameters;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import no.uib.probe.optprot.dataset.model.SearchingSubDataset;
import no.uib.probe.optprot.model.OptimisedSearchResults;
import no.uib.probe.optprot.model.ParameterScoreModel;
import no.uib.probe.optprot.model.RawScoreModel;
import no.uib.probe.optprot.model.SearchInputSetting;
import no.uib.probe.optprot.search.DefaultOptProtSearchOptimizer;
import no.uib.probe.optprot.search.SearchExecuter;
import no.uib.probe.optprot.util.MainUtilities;
import no.uib.probe.optprot.util.SpectraFileUtilities;

/**
 *
 * @author yfa041
 */
public class XTandemOptProtSearchOptimizer extends DefaultOptProtSearchOptimizer {

    /**
     * The identification file reader factory of compomics utilities.
     */
    private final IdfileReaderFactory readerFactory = IdfileReaderFactory.getInstance();
    /**
     * The compomics PTM factory.
     */
    private final ModificationFactory ptmFactory = ModificationFactory.getInstance();

    private final SearchingSubDataset optProtDataset;
    private final SearchInputSetting searchInputSetting;
    private final File identificationParametersFile;
    private final OptimisedSearchResults optimisedSearchResults;
    private final IdentificationParameters identificationParameters;
    private final Map<String, TreeSet<ParameterScoreModel>> parameterScoreMap;

    public XTandemOptProtSearchOptimizer(SearchingSubDataset optProtDataset, SearchInputSetting searchInputSetting, File identificationParametersFile) throws IOException {

        this.optProtDataset = optProtDataset;
        this.searchInputSetting = searchInputSetting;
        this.identificationParametersFile = identificationParametersFile;
        this.identificationParameters = IdentificationParameters.getIdentificationParameters(identificationParametersFile);

        System.out.println("identificationParametersFile " + identificationParametersFile.getAbsolutePath());
        System.out.println("identificationParametersFile " + identificationParameters.getSearchParameters().getModificationParameters().getAllModifications());
        this.optimisedSearchResults = new OptimisedSearchResults();
        this.parameterScoreMap = new LinkedHashMap<>();
        MainUtilities.cleanOutputFolder();
        parameterScoreMap.put("DigestionParameter", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("EnzymeParameter", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("SpecificityParameter", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("MaxMissCleavagesParameter", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("FragmentIonTypesParameter", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("PrecursorToleranceParameter", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("FragmentToleranceParameter", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("PrecursorChargeParameter", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("IsotopParameter", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("ModificationsParameter", new TreeSet<>(Collections.reverseOrder()));

        parameterScoreMap.put("XtandemSpectrumDynamicRange", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("XtandemNumberOfPeaks", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("XtandemMinimumFragmentMz", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("XtandemMinimumPeaks", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("XtandemNoiseSuppression", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("XtandemParentIsotopExpansion", new TreeSet<>(Collections.reverseOrder()));

        parameterScoreMap.put("XtandemQuickAcetyl", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("XtandemQuickPyrolidone", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("XtandemStPBias", new TreeSet<>(Collections.reverseOrder()));

        parameterScoreMap.put("XtandemUseRefine", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("XtandemUnanticipatedCleavage", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("XtandemRefineSimiEnzymaticCleavage", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("XtandemPotintialModification", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("XtandemPointMutations", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("XtandemSnAPs", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("XtandemSpectrumSynthesis", new TreeSet<>(Collections.reverseOrder()));

//        XtandemParameters xtandemParameters = (XtandemParameters) identificationParameters.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
//        xtandemParameters.setQuickPyrolidone(false);
//        xtandemParameters.setProteinQuickAcetyl(false);
//        xtandemParameters.setStpBias(false);
//        IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
    }

    private String digestionParameterOpt;

    public void startProcess(List<String> paramOrder) throws IOException {
        System.out.println("thershuld for id is " + optProtDataset.getAcceptedIDRatioThreshold());
        if (searchInputSetting.isOptimizeAllParameters()) {
//            advancedParam = false;

        }
        digestionParameterOpt = identificationParameters.getSearchParameters().getDigestionParameters().getCleavageParameter().name();
        for (String param : paramOrder) {

//        System.exit(0);
            System.out.println("-------------------------------------------param " + param + "-------------------------------------------");
            if (param.equalsIgnoreCase("DigestionParameter_1") && searchInputSetting.isOptimizeDigestionParameter()) {
                optimisedSearchResults.setDigestionParameter("enzyme");

                String value = this.optimizeEnzymeParameter(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("EnzymeParameter"));
                System.out.println("optimizeEnzymeParameter " + value + "------------------------------------------------------------------------->>> 2 id rate " + optProtDataset.getActiveIdentificationNum());
                optimisedSearchResults.setEnzymeName(value);

                if (!value.equalsIgnoreCase(identificationParameters.getSearchParameters().getDigestionParameters().getEnzymes().get(0).getName())) {
                    int nMissesCleavages = identificationParameters.getSearchParameters().getDigestionParameters().getnMissedCleavages(value);
                    identificationParameters.getSearchParameters().getDigestionParameters().clearEnzymes();
                    identificationParameters.getSearchParameters().getDigestionParameters().addEnzyme(EnzymeFactory.getInstance().getEnzyme(value));
                    identificationParameters.getSearchParameters().getDigestionParameters().setnMissedCleavages(value, nMissesCleavages);
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
                }
                value = this.optimizeSpecificityParameter(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("SpecificityParameter"));
                if (!value.equalsIgnoreCase(identificationParameters.getSearchParameters().getDigestionParameters().getSpecificity(identificationParameters.getSearchParameters().getDigestionParameters().getEnzymes().get(0).getName()).name())) {
                    identificationParameters.getSearchParameters().getDigestionParameters().setSpecificity(identificationParameters.getSearchParameters().getDigestionParameters().getEnzymes().get(0).getName(), DigestionParameters.Specificity.valueOf(value));
                    System.out.println("optimizeSpecificityParameter " + value + "------------------------------------------------------------------------->>> 3 id rate " + optProtDataset.getActiveIdentificationNum());
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
                }
                continue;

            }

            if (param.equalsIgnoreCase("FragmentIonTypesParameter") && searchInputSetting.isOptimizeFragmentIonTypesParameter()) {
                String value = this.optimizeFragmentIonTypesParameter(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("FragmentIonTypesParameter"));
                int forward = Integer.parseInt(value.split("-")[0]);
                int rewind = Integer.parseInt(value.split("-")[1]);
                boolean update = false;
                if (!identificationParameters.getSearchParameters().getForwardIons().contains(forward)) {
                    ArrayList<Integer> forwardIonsList = new ArrayList<>();
                    forwardIonsList.add(forward);
                    identificationParameters.getSearchParameters().setForwardIons(forwardIonsList);
                    update = true;
                }
                if (!identificationParameters.getSearchParameters().getRewindIons().contains(rewind)) {
                    ArrayList<Integer> rewindIonsList = new ArrayList<>();
                    rewindIonsList.add(rewind);
                    identificationParameters.getSearchParameters().setRewindIons(rewindIonsList);
                    update = true;
                }
                if (update) {
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
                }
                System.out.println("OptimizeFragmentIonTypesParameter" + value + "------------------------------------------------------------------------->>> 4 id rate " + optProtDataset.getActiveIdentificationNum());
                continue;

            }

//confusing param
            if (param.equalsIgnoreCase("DigestionParameter_2") && searchInputSetting.isOptimizeDigestionParameter()) {
                int value = this.optimizeMaxMissCleavagesParameter(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("MaxMissCleavagesParameter"));
                if (value != identificationParameters.getSearchParameters().getDigestionParameters().getnMissedCleavages(identificationParameters.getSearchParameters().getDigestionParameters().getEnzymes().get(0).getName())) {
                    identificationParameters.getSearchParameters().getDigestionParameters().setnMissedCleavages(identificationParameters.getSearchParameters().getDigestionParameters().getEnzymes().get(0).getName(), value);
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
                    System.out.println("optimizeMaxMissCleavagesParameter " + value + "------------------------------------------------------------------------->>> 6 id rate " + optProtDataset.getActiveIdentificationNum());
                }

                digestionParameterOpt = this.optimizeDigestionParameter(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("DigestionParameter"));

                continue;

            }

            if (param.equalsIgnoreCase("FragmentToleranceParameter") && searchInputSetting.isOptimizeFragmentToleranceParameter()) {
                double value = this.optimizeFragmentToleranceParameter(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("FragmentToleranceParameter"));
                if (value != identificationParameters.getSearchParameters().getFragmentIonAccuracy()) {
                    identificationParameters.getSearchParameters().setFragmentIonAccuracy(value);
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
                }
                System.out.println("optimizeFragmentToleranceParameter" + value + "------------------------------------------------------------------------->>> 8 id rate " + optProtDataset.getActiveIdentificationNum());
                continue;
            }

            if (param.equalsIgnoreCase("PrecursorChargeParameter") && searchInputSetting.isOptimizePrecursorChargeParameter()) {

                int[] values = this.optimizePrecursorChargeParameter(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("PrecursorChargeParameter"));
                if (values[1] != identificationParameters.getSearchParameters().getMaxChargeSearched()) {
                    identificationParameters.getSearchParameters().setMinChargeSearched(values[0]);
                    identificationParameters.getSearchParameters().setMaxChargeSearched(values[1]);
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
                    System.out.println("optimizePrecursorChargeParameter" + values[0] + "--" + values[1] + "------------------------------------------------------------------------->>> 9 id rate " + optProtDataset.getActiveIdentificationNum());

                }
                continue;
            }

//            if (param.equalsIgnoreCase("ModificationParameter") && searchInputSetting.isOptimizeModificationParameter()) {
//
//                Map<String, Set<String>> modificationsResults = this.optimizeModificationsParameter(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("ModificationsParameter"));
//
//                identificationParameters.getSearchParameters().getModificationParameters().clearFixedModifications();
//                identificationParameters.getSearchParameters().getModificationParameters().clearVariableModifications();
//                identificationParameters.getSearchParameters().getModificationParameters().clearRefinementModifications();
//                identificationParameters.getSearchParameters().getModificationParameters().getRefinementFixedModifications().clear();
//                for (String fixedMod : modificationsResults.get("fixedModifications")) {
//                    if (ptmFactory.getModification(fixedMod) != null) {
//                        identificationParameters.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(fixedMod));
//
//                        identificationParameters.getSearchParameters().getModificationParameters().addRefinementFixedModification(ptmFactory.getModification(fixedMod));
//                    }
//                }
//                for (String variableMod : modificationsResults.get("variableModifications")) {
//                    if (ptmFactory.getModification(variableMod) != null) {
//                        identificationParameters.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(variableMod));
//                    }
//                }
//                for (String refinmentVariableMod : modificationsResults.get("refinmentVariableModifications")) {
//                    if (ptmFactory.getModification(refinmentVariableMod) != null) {
//                        identificationParameters.getSearchParameters().getModificationParameters().addRefinementVariableModification(ptmFactory.getModification(refinmentVariableMod));
//                    }
//                }
//                IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
//                continue;
//            }
            if (param.equalsIgnoreCase("ModificationParameter") && searchInputSetting.isOptimizeModificationParameter()) {
                Map<String, Set<String>> modificationsResults = this.optimizeModificationsParameter(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("ModificationsParameter"));
                identificationParameters.getSearchParameters().getModificationParameters().clearFixedModifications();
                identificationParameters.getSearchParameters().getModificationParameters().clearVariableModifications();
                identificationParameters.getSearchParameters().getModificationParameters().clearRefinementModifications();
                identificationParameters.getSearchParameters().getModificationParameters().getRefinementFixedModifications().clear();
                for (String fixedMod : modificationsResults.get("fixedModifications")) {
                    if (ptmFactory.getModification(fixedMod) != null) {
                        identificationParameters.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(fixedMod));
                        identificationParameters.getSearchParameters().getModificationParameters().addRefinementFixedModification(ptmFactory.getModification(fixedMod));
                    }
                }
                for (String variableMod : modificationsResults.get("variableModifications")) {
                    if (ptmFactory.getModification(variableMod) != null) {
                        identificationParameters.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(variableMod));
                    }
                }
                IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
                continue;
            }

            if (param.equalsIgnoreCase("XtandemAdvancedParameter") && searchInputSetting.isOptimizeXtandemAdvancedParameter()) {
                advancedParam = true;
                XtandemParameters xtandemParameters = (XtandemParameters) identificationParameters.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
                xtandemParameters.setProteinQuickAcetyl(false);
                xtandemParameters.setQuickPyrolidone(false);
                xtandemParameters.setStpBias(false);
//                if (searchInputSetting.getXtandemOptProtAdvancedSearchParameters().isOptQuickAcetyl() || searchInputSetting.isOptimizeAllParameters()) {
//                boolean bvalue = optimizeQuickAcetyl(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("XtandemQuickAcetyl"));
//                if (bvalue) {
//                    xtandemParameters.setProteinQuickAcetyl(bvalue);
//                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
//                }
////                }
////                if (searchInputSetting.getXtandemOptProtAdvancedSearchParameters().isOptQuickPyrolidone() || searchInputSetting.isOptimizeAllParameters()) {
//
//                if (bvalue) {
//                    bvalue = optimizeQuickPyrolidone(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("XtandemQuickPyrolidone"));
//                    xtandemParameters.setQuickPyrolidone(bvalue);
//                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
//                }
////                if (searchInputSetting.getXtandemOptProtAdvancedSearchParameters().isOptstPBias() || searchInputSetting.isOptimizeAllParameters()) {
//                bvalue = optimizeStPBias(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("XtandemStPBias"));
//                if (bvalue) {
//                    xtandemParameters.setStpBias(bvalue);
//                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
//                }
////                if (searchInputSetting.getXtandemOptProtAdvancedSearchParameters().isOptSpectrumDynamicRange() || searchInputSetting.isOptimizeAllParameters()) {
//                double dvalue = optimizeSpectrumDynamicRange(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("XtandemSpectrumDynamicRange"));
//                if (dvalue != xtandemParameters.getDynamicRange()) {
//                    xtandemParameters.setDynamicRange(dvalue);
//                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
//                }
////                if (searchInputSetting.getXtandemOptProtAdvancedSearchParameters().isOptSpectrumNumbPeaks() || searchInputSetting.isOptimizeAllParameters()) {
//                int ivalue = optimizeSpectrumPeaksNumber(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("XtandemNumberOfPeaks"));
//                if (ivalue != xtandemParameters.getnPeaks()) {
//                    xtandemParameters.setnPeaks(ivalue);
//                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
//                }
////                if (searchInputSetting.getXtandemOptProtAdvancedSearchParameters().isOptSpectrumMinimumFragment() || searchInputSetting.isOptimizeAllParameters()) {
//                dvalue = optimizeMinimumFragmentMz(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("XtandemMinimumFragmentMz"));
//                if (dvalue != xtandemParameters.getMinFragmentMz()) {
//                    xtandemParameters.setMinFragmentMz(dvalue);
//                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
//                }
////                if (searchInputSetting.getXtandemOptProtAdvancedSearchParameters().isOptSpectrumMinimumPeaks() || searchInputSetting.isOptimizeAllParameters()) {
//                ivalue = optimizeMinimumPeaks(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("XtandemMinimumPeaks"));
//                if (ivalue != xtandemParameters.getMinPeaksPerSpectrum()) {
//                    xtandemParameters.setMinPeaksPerSpectrum(ivalue);
//                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
//                }
////                if (searchInputSetting.getXtandemOptProtAdvancedSearchParameters().isOptNoiseSuppression() || searchInputSetting.isOptimizeAllParameters()) {
//                double minPrecursorMass = optimizeNoiseSuppression(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("XtandemNoiseSuppression"));
//                bvalue = (minPrecursorMass == 0.0);
//                if (bvalue != xtandemParameters.isUseNoiseSuppression()) {
//                    xtandemParameters.setUseNoiseSuppression(bvalue);
//                }
//                if (!bvalue) {
//                    xtandemParameters.setMinPrecursorMass(minPrecursorMass);
//                }
//                IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
////                }
////                if (xtandemParameters.isRefine() || searchInputSetting.getXtandemOptProtAdvancedSearchParameters().isOptUseRefine() || searchInputSetting.isOptimizeAllParameters()) {
//
////                    if (searchInputSetting.getXtandemOptProtAdvancedSearchParameters().isOptRefineSimiEnzymaticCleavage() || searchInputSetting.isOptimizeAllParameters()) {
//                bvalue = optimizeRefineSimiEnzymaticCleavage(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("XtandemRefineSimiEnzymaticCleavage"));
//                if (bvalue != xtandemParameters.isRefineSemi()) {
//                    xtandemParameters.setRefineSemi(bvalue);
//                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
//                }
////                    if (searchInputSetting.getXtandemOptProtAdvancedSearchParameters().isOptRefinePotintialModification() || searchInputSetting.isOptimizeAllParameters()) {
//                bvalue = optimizePotintialModification(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("XtandemPotintialModification"));
//                if (bvalue != xtandemParameters.isPotentialModificationsForFullRefinment()) {
//                    xtandemParameters.setPotentialModificationsForFullRefinment(bvalue);
//                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
//                }
////                    if (searchInputSetting.getXtandemOptProtAdvancedSearchParameters().isOptRefinePointMutations() || searchInputSetting.isOptimizeAllParameters()) {
//                bvalue = optimizeRefinePointMutations(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("XtandemPointMutations"));
//                if (bvalue != xtandemParameters.isRefinePointMutations()) {
//                    xtandemParameters.setRefinePointMutations(bvalue);
//                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
//                }
////                    if (searchInputSetting.getXtandemOptProtAdvancedSearchParameters().isOptRefineSnAPs() || searchInputSetting.isOptimizeAllParameters()) {
//                bvalue = optimizeRefineSnAPs(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("XtandemSnAPs"));
//                if (bvalue != xtandemParameters.isRefineSnaps()) {
//                    xtandemParameters.setRefineSnaps(bvalue);
//                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
//                }
////                    if (searchInputSetting.getXtandemOptProtAdvancedSearchParameters().isOptRefineSpectrumSynthesis() || searchInputSetting.isOptimizeAllParameters()) {
//                bvalue = optimizeRefineSpectrumSynthesis(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("XtandemSpectrumSynthesis"));
//                if (bvalue != xtandemParameters.isRefineSpectrumSynthesis()) {
//                    xtandemParameters.setRefineSpectrumSynthesis(bvalue);
//                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
//                }
////                    if (searchInputSetting.getXtandemOptProtAdvancedSearchParameters().isOptRefineUnanticipatedCleavage() || searchInputSetting.isOptimizeAllParameters()) {
//                bvalue = optimizeRefineUnanticipatedCleavage(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("XtandemUnanticipatedCleavage"));
//                if (bvalue != xtandemParameters.isRefineUnanticipatedCleavages()) {
//                    xtandemParameters.setRefineUnanticipatedCleavages(bvalue);
//                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
//                }
////                }
//
////                if (searchInputSetting.getXtandemOptProtAdvancedSearchParameters().isOptUseRefine() || searchInputSetting.isOptimizeAllParameters()) {
//                bvalue = optimizeUseRefine(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("XtandemUseRefine"));
//                if (bvalue != xtandemParameters.isRefine()) {
//                    xtandemParameters.setRefine(bvalue);
//                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
//                }
////                if (searchInputSetting.getXtandemOptProtAdvancedSearchParameters().isOptParentIsotopeExpansion() || searchInputSetting.isOptimizeAllParameters()) {
//                bvalue = optimizeParentIsotopExpansion(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("XtandemParentIsotopExpansion"));
//                if (bvalue != xtandemParameters.getParentMonoisotopicMassIsotopeError()) {
//                    xtandemParameters.setParentMonoisotopicMassIsotopeError(bvalue);
//                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
//                }
                continue;
            }
            if (param.equalsIgnoreCase("PrecursorToleranceParameter") && searchInputSetting.isOptimizePrecursorToleranceParameter()) {

                double value = this.optimizePrecursorToleranceParameter(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("PrecursorToleranceParameter"));
                if (value != identificationParameters.getSearchParameters().getPrecursorAccuracy()) {
                    identificationParameters.getSearchParameters().setPrecursorAccuracy(value);
                    if (value > 1) {
                        identificationParameters.getSearchParameters().setPrecursorAccuracyType(SearchParameters.MassAccuracyType.PPM);
                    } else {
                        identificationParameters.getSearchParameters().setPrecursorAccuracyType(SearchParameters.MassAccuracyType.DA);
                    }
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
                    System.out.println("OptimizePrecursorToleranceParameter" + value + "------------------------------------------------------------------------->>> 5 id rate " + optProtDataset.getActiveIdentificationNum());
                }
            }

        }
        if (!digestionParameterOpt.equalsIgnoreCase(identificationParameters.getSearchParameters().getDigestionParameters().getCleavageParameter().name())) {
            System.out.println("optimizeDigestionParameter " + digestionParameterOpt + "------------------------------------------------------------------------->>> 7 id rate " + optProtDataset.getActiveIdentificationNum() + " results:  " + optProtDataset.getTempIdentificationNum());
            optimisedSearchResults.setDigestionParameter(digestionParameterOpt);
            identificationParameters.getSearchParameters().getDigestionParameters().setCleavageParameter(DigestionParameters.CleavageParameter.valueOf(digestionParameterOpt));
            IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
//                digestionParameterOpt = "";
        }
//        ReportExporter.printFullReport(identificationParameters, optProtDataset, Advocate.xtandem);
//        for (String key : parameterScoreMap.keySet()) {
//            System.out.println("at param map " + key);
//            for (ParameterScoreModel m : parameterScoreMap.get(key)) {
//                System.out.println(" " + m);
//            }
//            System.out.println("---------------------------------------");
//        }
//        System.exit(0);

    }
    private boolean advancedParam = true;

    @Override
    public synchronized RawScoreModel excuteSearch(SearchingSubDataset optProtDataset, String defaultOutputFileName, String paramOption, IdentificationParameters tempIdParam, boolean addSpectraList, SearchInputSetting optProtSearchSettings, File identificationParametersFile,boolean pairData) {
//        try {
        if (!optProtSearchSettings.getXTandemEnabledParameters().getParamsToOptimize().isEnabledParam(paramOption.split("_")[0])) {
            System.out.println("param " + paramOption + " is not supported " + paramOption);
            return new RawScoreModel();
        }
        if (!advancedParam && tempIdParam.getSearchParameters().getModificationParameters().getFixedModifications().size() == 1) {
            boolean terminalMod = ptmFactory.getModification(tempIdParam.getSearchParameters().getModificationParameters().getFixedModifications().get(0)).getModificationType().isCTerm() || ptmFactory.getModification(tempIdParam.getSearchParameters().getModificationParameters().getFixedModifications().get(0)).getModificationType().isNTerm();
            if (terminalMod) {
                System.out.println("terminal mod only ");
                return new RawScoreModel();
            }

        }
        if (defaultOutputFileName.contains("_resultsf_Carbamilation of protein N-term")||defaultOutputFileName.contains("resultsf_Acetylation of protein N-term")) {
            System.out.println("param " + paramOption + " is not supported " + paramOption);
            return new RawScoreModel();
        }
        
        
        
        SearchParameters searchParameters = tempIdParam.getSearchParameters();
        if (!advancedParam) {
            XtandemParameters xtandemParameters = (XtandemParameters) searchParameters.getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
            xtandemParameters.setQuickPyrolidone(false);
            xtandemParameters.setProteinQuickAcetyl(false);
            xtandemParameters.setStpBias(false);
            xtandemParameters.setRefine(false);
        }
        File resultOutput = SearchExecuter.executeSearch(defaultOutputFileName, optProtSearchSettings, optProtDataset.getSubMsFile(), optProtDataset.getSubFastaFile(), tempIdParam, identificationParametersFile);
        final List<SpectrumMatch> validatedMaches = SpectraFileUtilities.getValidatedIdentificationResults(resultOutput, optProtDataset.getSubMsFile(), Advocate.xtandem, tempIdParam);
        RawScoreModel rawScore = SpectraFileUtilities.getComparableRawScore(optProtDataset, validatedMaches, Advocate.xtandem,pairData);//(optProtDataset, resultOutput, optProtDataset.getSubMsFile(), Advocate.sage, tempIdParam, updateDataReference);

        MainUtilities.deleteFolder(resultOutput);
        if (addSpectraList && rawScore.isSignificatChange()) {
            rawScore.setSpectrumMatchResult(validatedMaches);
        }
        return (rawScore);

//            
//            if (addPeptideMasses) {
//                File xTandemFile = new File(resultOutput, SearchHandler.getXTandemFileName(IoUtil.removeExtension(optProtDataset.getSubMsFile().getName())));
//                IdfileReader idReader = readerFactory.getFileReader(xTandemFile);
//                SequenceMatchingParameters modificationSequenceMatchingParameters = tempIdParam.getModificationLocalizationParameters().getSequenceMatchingParameters();
//                FMIndex sequenceProvider = new FMIndex(optProtDataset.getSubFastaFile(), null, new OptProtWaitingHandler(), false, tempIdParam);
//                for (SpectrumMatch sm : validatedMaches) {
//                    Thread t = new Thread(() -> {
//                        for (PeptideAssumption pepAss : sm.getAllPeptideAssumptions().toList()) {
//                            try {
//                                Peptide pep = pepAss.getPeptide();
//                                ModificationLocalizationMapper.modificationLocalization(
//                                        pep,
//                                        tempIdParam,
//                                        idReader,
//                                        ptmFactory,
//                                        sequenceProvider
//                                );
//                                pepAss.getPeptide().setMass(pep.getMass(searchParameters.getModificationParameters(), sequenceProvider, modificationSequenceMatchingParameters));
//                                sm.setBestPeptideAssumption(pepAss);
//                            } catch (Exception e) {
//                            }
//                        }
//                    });
//                    t.start();
//                    while (t.isAlive()) {
//                        Thread.sleep(10);
//                    }
//
//                }
//            }
//            MainUtilities.deleteFolder(resultOutput);
//            return (validatedMaches);
//        } catch (IOException | InterruptedException ex) {
//            identificationParametersFile.delete();
//            System.out.println("error thrown here ----------------- " + paramOption);
//            ex.printStackTrace();
//            System.exit(0);
//        }
//        return new ArrayList<>();
    }
//
//    public double optimizeSpectrumDynamicRange(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
//        Map<Double, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        int idRate = optProtDataset.getActiveIdentificationNum();
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
//        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
//        double selectedOption = xtandemParameters.getDynamicRange();
//
//        for (double i = 100.0; i > 0;) {
//            final String option = "spectrumDR_" + i;
//            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//            xtandemParameters = (XtandemParameters) tempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
//            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
//            xtandemParameters.setDynamicRange(i);
//            final double j = i;
//            final ParameterScoreModel paramScore = new ParameterScoreModel();
//            paramScore.setParamId("spectrumDR");
//
//            Future future = executor.submit(() -> {
//                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
//                paramScore.setScore(resultsMap.get(j));
//                paramScore.setParamValue(option);
//                parameterScoreSet.add(paramScore);
//            });
//            while (!future.isDone()) {
//            }
//            i -= 20;
//        }
//
//        int localId = -1;
//        double localSelection = 0;
//        for (double dRangeScore : resultsMap.keySet()) {
////            System.out.println(" option 1 " + dRangeScore + "  " + resultsMap.get(dRangeScore) + " > " + localId + "  " + idRate);
//            if (resultsMap.get(dRangeScore) > localId) {
//                localId = resultsMap.get(dRangeScore);
//                localSelection = dRangeScore;
//
//            }
//        }
//        if (localId >= idRate) {
//            selectedOption = localSelection;
//            optProtDataset.setActiveIdentificationNum(localId);
//        }
//        return selectedOption;
//    }
//
//    public int optimizeSpectrumPeaksNumber(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
//
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
//        Integer selectedOption = xtandemParameters.getnPeaks();
//        int idRate = optProtDataset.getActiveIdentificationNum();
//        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
//
//        for (int i = 30; i <= 70;) {
//            final String option = "peaksNum_" + i;
//            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//            xtandemParameters = (XtandemParameters) tempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
//            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
//            xtandemParameters.setnPeaks(i);
//            final int j = i;
//            final ParameterScoreModel paramScore = new ParameterScoreModel();
//            paramScore.setParamId("peaksNum");
//
//            Future future = executor.submit(() -> {
//                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
//                paramScore.setScore(resultsMap.get(j));
//                paramScore.setParamValue(option);
//                parameterScoreSet.add(paramScore);
//            });
//            while (!future.isDone()) {
//            }
//            i += 10;
//        }
//        int localId = -1;
//        int localSelection = 0;
//        for (int option : resultsMap.keySet()) {
////            System.out.println(" option #peaks " + option + "  " + resultsMap.get(option) + " > " + localId + "  " + idRate);
//            if (resultsMap.get(option) > localId) {
//                localId = resultsMap.get(option);
//                localSelection = option;
//
//            }
//        }
//        if (localId >= idRate) {
//            selectedOption = localSelection;
//            optProtDataset.setActiveIdentificationNum(localId);
//        }
//        return selectedOption;
//    }
//
//    public double optimizeMinimumFragmentMz(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
//        Map<Double, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        int idRate = optProtDataset.getActiveIdentificationNum();
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
//        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
//        double selectedOption = xtandemParameters.getMinFragmentMz();
//        for (int i = 100; i <= 300;) {
//            final String option = "minimumFragmentMz_" + i;
//            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//            xtandemParameters = (XtandemParameters) tempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
//            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
//            xtandemParameters.setMinFragmentMz(i);
//            final double j = i;
//            final ParameterScoreModel paramScore = new ParameterScoreModel();
//            paramScore.setParamId("minimumFragmentMz");
//
//            Future future = executor.submit(() -> {
//                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
//                paramScore.setScore(resultsMap.get(j));
//                paramScore.setParamValue(option);
//                parameterScoreSet.add(paramScore);
//            });
//            while (!future.isDone()) {
//            }
//            i += 50;
//        }
//        int localId = -1;
//        double localSelection = 0;
//        for (double dRangeScore : resultsMap.keySet()) {
//            System.out.println(" option 1 " + dRangeScore + "  " + resultsMap.get(dRangeScore) + " > " + localId + "  " + idRate);
//            if (resultsMap.get(dRangeScore) > localId) {
//                localId = resultsMap.get(dRangeScore);
//                localSelection = dRangeScore;
//
//            }
//        }
//        if (localId >= idRate) {
//            selectedOption = localSelection;
//            optProtDataset.setActiveIdentificationNum(localId);
//        }
//
//        return selectedOption;
//    }
//
//    public int optimizeMinimumPeaks(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
//        Integer selectedOption = xtandemParameters.getMinPeaksPerSpectrum();
//        int idRate = optProtDataset.getActiveIdentificationNum();
//        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
//        int lastValue = 0;
//        for (int i = 5; i <= 100;) {
//            final String option = "minpeaksNum_" + i;
//            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//            xtandemParameters = (XtandemParameters) tempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
//            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
//            xtandemParameters.setMinPeaksPerSpectrum(i);
//            final int j = i;
//            final ParameterScoreModel paramScore = new ParameterScoreModel();
//            paramScore.setParamId("minpeaksNum");
//
//            Future future = executor.submit(() -> {
//                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
//                paramScore.setScore(resultsMap.get(j));
//                paramScore.setParamValue(option);
//                parameterScoreSet.add(paramScore);
//            });
//            while (!future.isDone()) {
//            }
//            if (lastValue > resultsMap.get(i)) {
//                break;
//            }
//            lastValue = resultsMap.get(i);
//            i += 10;
//
//        }
//        int localId = -1;
//        int localSelection = 0;
//        for (int option : resultsMap.keySet()) {
//            if (resultsMap.get(option) > localId) {
//                localId = resultsMap.get(option);
//                localSelection = option;
//
//            }
//        }
//        if (localId >= idRate) {
//            selectedOption = localSelection;
//            optProtDataset.setActiveIdentificationNum(localId);
//        }
//        return selectedOption;
//    }
//
//    public double optimizeNoiseSuppression(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
//        Map<Double, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        int idRate = optProtDataset.getActiveIdentificationNum();
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
//        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
//        double selectedOption = xtandemParameters.getMinPrecursorMass();
//        final String option = "noiseSupression_" + false;
//        final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
//        xtandemParameters.setUseNoiseSuppression(false);
//        final ParameterScoreModel paramScore = new ParameterScoreModel();
//        paramScore.setParamId("noiseSupression");
//
//        Future future = executor.submit(() -> {
//            resultsMap.put(0.0, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
//            paramScore.setScore(resultsMap.get(0.0));
//            paramScore.setParamValue(option);
//            parameterScoreSet.add(paramScore);
//        });
//        while (!future.isDone()) {
//        }
//        xtandemParameters.setUseNoiseSuppression(true);
//        for (double j = 500; j < 1600;) {
//            final String suboption = "noiseSupression_" + true + "_" + j;
//            final String subupdatedName = Configurations.DEFAULT_RESULT_NAME + "_" + suboption + "_" + msFileName;
//            final double i = j;
//            xtandemParameters.setMinPrecursorMass(j);
//            final ParameterScoreModel paramScore2 = new ParameterScoreModel();
//            paramScore2.setParamId("noiseSupression");
//            future = executor.submit(() -> {
//                resultsMap.put(i, excuteSearch(optProtDataset, subupdatedName, suboption, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
//                paramScore2.setScore(resultsMap.get(i));
//                paramScore2.setParamValue(option);
//                parameterScoreSet.add(paramScore2);
//            });
//            while (!future.isDone()) {
//            }
//            j += 350;
//
//        }
//
//        int localId = -1;
//        double localSelection = 0;
//        for (double dRangeScore : resultsMap.keySet()) {
////            System.out.println(" option 1 " + option + "  " + resultsMap.get(option) + " > " + localId + "  " + idRate);
//            if (resultsMap.get(dRangeScore) > localId) {
//                localId = resultsMap.get(dRangeScore);
//                localSelection = dRangeScore;
//
//            }
//        }
//        if (localId >= idRate) {
//            selectedOption = localSelection;
//            optProtDataset.setActiveIdentificationNum(localId);
//        }
//
//        return selectedOption;
//    }
////
//
//    public boolean optimizeParentIsotopExpansion(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
//        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        int idRate = optProtDataset.getActiveIdentificationNum();
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
//        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
//        boolean selectedOption = xtandemParameters.getParentMonoisotopicMassIsotopeError();
//
//        for (int i = 0; i < 2; i++) {
//            final String option = "parentMonoisotopicMassIsotopeError_" + (i == 1);
//            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
//            xtandemParameters.setParentMonoisotopicMassIsotopeError(i == 1);
//            final int j = i;
//            final ParameterScoreModel paramScore = new ParameterScoreModel();
//            paramScore.setParamId("parentMonoisotopicMassIsotopeError");
//
//            Future future = executor.submit(() -> {
//                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
//                paramScore.setScore(resultsMap.get(j));
//                paramScore.setParamValue(option);
//                parameterScoreSet.add(paramScore);
//            });
//            while (!future.isDone()) {
//            }
////            System.out.println("at parent istop " + option + " " + resultsMap.get(i));
//        }
//
//        int localId = -1;
//        int localSelection = 0;
//        for (int dRangeScore : resultsMap.keySet()) {
//            if (resultsMap.get(dRangeScore) > localId) {
//                localId = resultsMap.get(dRangeScore);
//                localSelection = dRangeScore;
//
//            }
//        }
//        if (localId >= idRate) {
//            selectedOption = localSelection == 1;
//            optProtDataset.setActiveIdentificationNum(localId);
//        }
//        return selectedOption;
//    }
////
//
//    public boolean optimizeQuickAcetyl(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
//        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        int idRate = optProtDataset.getActiveIdentificationNum();
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
//        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
//        boolean selectedOption = xtandemParameters.isProteinQuickAcetyl();
//
//        for (int i = 0; i < 2; i++) {
//            boolean useQuickAcetyl = (i == 1);
//            final String option = "useQuickAcetyl_" + useQuickAcetyl;
//            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
//            xtandemParameters.setProteinQuickAcetyl(useQuickAcetyl);
//            final int j = i;
//            final ParameterScoreModel paramScore = new ParameterScoreModel();
//            paramScore.setParamId("QuickAcetyl");
//            Future future = executor.submit(() -> {
//                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
//                paramScore.setScore(resultsMap.get(j));
//                paramScore.setParamValue(option);
//                parameterScoreSet.add(paramScore);
//            });
//            while (!future.isDone()) {
//            }
////            System.out.println("at quick acetyle " + useQuickAcetyl + " " + resultsMap.get(j));
//        }
//
//        int localId = -1;
//        int localSelection = 0;
//        for (int dRangeScore : resultsMap.keySet()) {
//            if (resultsMap.get(dRangeScore) > localId) {
//                localId = resultsMap.get(dRangeScore);
//                localSelection = dRangeScore;
//
//            }
//        }
//        if (localId >= idRate) {
//            selectedOption = localSelection == 1;
//            optProtDataset.setActiveIdentificationNum(localId);
//        }
//        return selectedOption;
//    }
//
//    public boolean optimizeQuickPyrolidone(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
//        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        int idRate = optProtDataset.getActiveIdentificationNum();
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
//        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
//        boolean selectedOption = xtandemParameters.isQuickPyrolidone();
//
//        for (int i = 0; i < 2; i++) {
//            boolean useQuickPyrolidone = (i == 1);
//            final String option = "useQuickPyrolidone_" + useQuickPyrolidone;
//            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
//            xtandemParameters.setQuickPyrolidone(useQuickPyrolidone);
//            final int j = i;
//            final ParameterScoreModel paramScore = new ParameterScoreModel();
//            paramScore.setParamId("QuickPyrolidone");
//
//            Future future = executor.submit(() -> {
//                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
//                paramScore.setScore(resultsMap.get(j));
//                paramScore.setParamValue(option);
//                parameterScoreSet.add(paramScore);
//            });
//            while (!future.isDone()) {
//            }
//            System.out.println("at quick prolyien  " + useQuickPyrolidone + " " + resultsMap.get(j));
//        }
//
//        int localId = -1;
//        int localSelection = 0;
//        for (int dRangeScore : resultsMap.keySet()) {
//            if (resultsMap.get(dRangeScore) > localId) {
//                localId = resultsMap.get(dRangeScore);
//                localSelection = dRangeScore;
//
//            }
//        }
//        if (localId >= idRate) {
//            selectedOption = localSelection == 1;
//            optProtDataset.setActiveIdentificationNum(localId);
//        }
//
//        return selectedOption;
//    }
////
//
//    public boolean optimizeStPBias(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
//        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        int idRate = optProtDataset.getActiveIdentificationNum();
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
//        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
//        boolean selectedOption = xtandemParameters.isStpBias();
//
//        for (int i = 0; i < 2; i++) {
//            boolean useStpBias = (i == 1);
//            final String option = "useStpBias_" + useStpBias;
//            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
//            xtandemParameters.setStpBias(useStpBias);
//            final int j = i;
//            final ParameterScoreModel paramScore = new ParameterScoreModel();
//            paramScore.setParamId("StpBias");
//
//            Future future = executor.submit(() -> {
//                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
//                paramScore.setScore(resultsMap.get(j));
//                paramScore.setParamValue(option);
//                parameterScoreSet.add(paramScore);
//            });
//            while (!future.isDone()) {
//            }
//
////            System.out.println("at useStpBias  " + useStpBias + " " + resultsMap.get(j));
//        }
//
//        int localId = -1;
//        int localSelection = 0;
//        for (int dRangeScore : resultsMap.keySet()) {
//            if (resultsMap.get(dRangeScore) > localId) {
//                localId = resultsMap.get(dRangeScore);
//                localSelection = dRangeScore;
//
//            }
//        }
//        if (localId >= idRate) {
//            selectedOption = localSelection == 1;
//            optProtDataset.setActiveIdentificationNum(localId);
//        }
//
//        return selectedOption;
//    }
////
//
//    public boolean optimizeUseRefine(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
//        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        int idRate = optProtDataset.getActiveIdentificationNum();
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
//        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
//        boolean selectedOption = xtandemParameters.isRefine();
//
//        for (int i = 0; i < 2; i++) {
//            boolean useRefine = (i == 1);
//            final String option = "useRefine_" + useRefine;
//            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
//            xtandemParameters.setRefine(useRefine);
//            final int j = i;
//            final ParameterScoreModel paramScore = new ParameterScoreModel();
//            paramScore.setParamId("useRefineuseRefine");
//
//            Future future = executor.submit(() -> {
//                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
//                paramScore.setScore(resultsMap.get(j));
//                paramScore.setParamValue(option);
//                parameterScoreSet.add(paramScore);
//            });
//            while (!future.isDone()) {
//            }
//        }
//
//        int localId = -1;
//        int localSelection = 0;
//        for (int dRangeScore : resultsMap.keySet()) {
//            if (resultsMap.get(dRangeScore) > localId) {
//                localId = resultsMap.get(dRangeScore);
//                localSelection = dRangeScore;
//
//            }
//        }
//        if (localId >= idRate) {
//            selectedOption = localSelection == 1;
//            optProtDataset.setActiveIdentificationNum(localId);
//        }
//
//        return selectedOption;
//
//    }
////
//
//    public boolean optimizeRefineUnanticipatedCleavage(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
//        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        int idRate = optProtDataset.getActiveIdentificationNum();
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
//        xtandemParameters.setRefine(true);
//        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
//        boolean selectedOption = xtandemParameters.isRefineUnanticipatedCleavages();
//
//        for (int i = 0; i < 2; i++) {
//            boolean useRefineUnanticipatedCleavages = (i == 1);
//            final String option = "useRefineUnanticipatedCleavages_" + useRefineUnanticipatedCleavages;
//            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
//            xtandemParameters.setRefineUnanticipatedCleavages(useRefineUnanticipatedCleavages);
//            final int j = i;
//            final ParameterScoreModel paramScore = new ParameterScoreModel();
//            paramScore.setParamId("UnanticipatedCleavages");
//
//            Future future = executor.submit(() -> {
//                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
//                paramScore.setScore(resultsMap.get(j));
//                paramScore.setParamValue(option);
//                parameterScoreSet.add(paramScore);
//            });
//            while (!future.isDone()) {
//            }
//        }
//
//        int localId = -1;
//        int localSelection = 0;
//        for (int dRangeScore : resultsMap.keySet()) {
//            if (resultsMap.get(dRangeScore) > localId) {
//                localId = resultsMap.get(dRangeScore);
//                localSelection = dRangeScore;
//
//            }
//        }
//        if (localId >= idRate) {
//            selectedOption = localSelection == 1;
//            optProtDataset.setActiveIdentificationNum(localId);
//        }
//
//        return selectedOption;
//
//    }
////
//
//    public boolean optimizeRefineSimiEnzymaticCleavage(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
//        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        int idRate = optProtDataset.getActiveIdentificationNum();
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
//        xtandemParameters.setRefine(true);
//        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
//        boolean selectedOption = xtandemParameters.isRefineSemi();
//
//        for (int i = 0; i < 2; i++) {
//            boolean useRefineSimiEnzymaticCleavage = (i == 1);
//            final String option = "useRefineSimiEnzymaticCleavage_" + useRefineSimiEnzymaticCleavage;
//            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
//            xtandemParameters.setRefineSemi(useRefineSimiEnzymaticCleavage);
//            final int j = i;
//            final ParameterScoreModel paramScore = new ParameterScoreModel();
//            paramScore.setParamId("SimiEnzymaticCleavage");
//
//            Future future = executor.submit(() -> {
//                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
//                paramScore.setScore(resultsMap.get(j));
//                paramScore.setParamValue(option);
//                parameterScoreSet.add(paramScore);
//            });
//            while (!future.isDone()) {
//            }
//        }
//
//        int localId = -1;
//        int localSelection = 0;
//        for (int dRangeScore : resultsMap.keySet()) {
//            if (resultsMap.get(dRangeScore) > localId) {
//                localId = resultsMap.get(dRangeScore);
//                localSelection = dRangeScore;
//
//            }
//        }
//        if (localId >= idRate) {
//            selectedOption = localSelection == 1;
//            optProtDataset.setActiveIdentificationNum(localId);
//        }
//
//        return selectedOption;
//
//    }
////
//
//    public boolean optimizePotintialModification(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
//        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        int idRate = optProtDataset.getActiveIdentificationNum();
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
//        xtandemParameters.setRefine(true);
//        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
//        boolean selectedOption = xtandemParameters.isPotentialModificationsForFullRefinment();
//
//        for (int i = 0; i < 2; i++) {
//            boolean usePotintialModification = (i == 1);
//            final String option = "usePotintialModification_" + usePotintialModification;
//            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
//            xtandemParameters.setPotentialModificationsForFullRefinment(usePotintialModification);
//            final int j = i;
//            final ParameterScoreModel paramScore = new ParameterScoreModel();
//            paramScore.setParamId("PotintialModification");
//
//            Future future = executor.submit(() -> {
//                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
//                paramScore.setScore(resultsMap.get(j));
//                paramScore.setParamValue(option);
//                parameterScoreSet.add(paramScore);
//            });
//            while (!future.isDone()) {
//            }
//        }
//
//        int localId = -1;
//        int localSelection = 0;
//        for (int dRangeScore : resultsMap.keySet()) {
//            if (resultsMap.get(dRangeScore) > localId) {
//                localId = resultsMap.get(dRangeScore);
//                localSelection = dRangeScore;
//
//            }
//        }
//        if (localId >= idRate) {
//            selectedOption = localSelection == 1;
//            optProtDataset.setActiveIdentificationNum(localId);
//        }
//
//        return selectedOption;
//
//    }
////
//
//    public boolean optimizeRefinePointMutations(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
//        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        int idRate = optProtDataset.getActiveIdentificationNum();
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
//        xtandemParameters.setRefine(true);
//        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
//        boolean selectedOption = xtandemParameters.isRefinePointMutations();
//
//        for (int i = 0; i < 2; i++) {
//            boolean useRefinePointMutations = (i == 1);
//            final String option = "useRefinePointMutations_" + useRefinePointMutations;
//            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
//            xtandemParameters.setRefinePointMutations(useRefinePointMutations);
//            final int j = i;
//            final ParameterScoreModel paramScore = new ParameterScoreModel();
//            paramScore.setParamId("PointMutations");
//
//            Future future = executor.submit(() -> {
//                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
//                paramScore.setScore(resultsMap.get(j));
//                paramScore.setParamValue(option);
//                parameterScoreSet.add(paramScore);
//            });
//            while (!future.isDone()) {
//            }
//        }
//
//        int localId = -1;
//        int localSelection = 0;
//        for (int dRangeScore : resultsMap.keySet()) {
//            if (resultsMap.get(dRangeScore) > localId) {
//                localId = resultsMap.get(dRangeScore);
//                localSelection = dRangeScore;
//
//            }
//        }
//        if (localId >= idRate) {
//            selectedOption = localSelection == 1;
//            optProtDataset.setActiveIdentificationNum(localId);
//        }
//
//        return selectedOption;
//
//    }
////
////
//
//    public boolean optimizeRefineSnAPs(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
//        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        int idRate = optProtDataset.getActiveIdentificationNum();
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
//        xtandemParameters.setRefine(true);
//        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
//        boolean selectedOption = xtandemParameters.isRefineSnaps();
//
//        for (int i = 0; i < 2; i++) {
//            boolean useRefineSnAPs = (i == 1);
//            final String option = "useRefineSnAPs_" + useRefineSnAPs;
//            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
//            xtandemParameters.setRefineSnaps(useRefineSnAPs);
//            final int j = i;
//            final ParameterScoreModel paramScore = new ParameterScoreModel();
//            paramScore.setParamId("SnAPs");
//
//            Future future = executor.submit(() -> {
//                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
//                paramScore.setScore(resultsMap.get(j));
//                paramScore.setParamValue(option);
//                parameterScoreSet.add(paramScore);
//            });
//            while (!future.isDone()) {
//            }
//        }
//
//        int localId = -1;
//        int localSelection = 0;
//        for (int dRangeScore : resultsMap.keySet()) {
//            if (resultsMap.get(dRangeScore) > localId) {
//                localId = resultsMap.get(dRangeScore);
//                localSelection = dRangeScore;
//
//            }
//        }
//        if (localId >= idRate) {
//            selectedOption = localSelection == 1;
//            optProtDataset.setActiveIdentificationNum(localId);
//        }
//
//        return selectedOption;
//
//    }
//
//    public boolean optimizeRefineSpectrumSynthesis(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
//        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        int idRate = optProtDataset.getActiveIdentificationNum();
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
//        xtandemParameters.setRefine(true);
//        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
//        boolean selectedOption = xtandemParameters.isRefineSpectrumSynthesis();
//
//        for (int i = 0; i < 2; i++) {
//            boolean useRefineSpectrumSynthesis = (i == 1);
//            final String option = "useRefineSpectrumSynthesis_" + useRefineSpectrumSynthesis;
//            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
//            xtandemParameters.setRefineSpectrumSynthesis(useRefineSpectrumSynthesis);
//            final int j = i;
//            final ParameterScoreModel paramScore = new ParameterScoreModel();
//            paramScore.setParamId("SpectrumSynthesis");
//
//            Future future = executor.submit(() -> {
//                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
//                paramScore.setScore(resultsMap.get(j));
//                paramScore.setParamValue(option);
//                parameterScoreSet.add(paramScore);
//            });
//            while (!future.isDone()) {
//            }
//        }
//
//        int localId = -1;
//        int localSelection = 0;
//        for (int dRangeScore : resultsMap.keySet()) {
//            if (resultsMap.get(dRangeScore) > localId) {
//                localId = resultsMap.get(dRangeScore);
//                localSelection = dRangeScore;
//
//            }
//        }
//        if (localId >= idRate) {
//            selectedOption = localSelection == 1;
//            optProtDataset.setActiveIdentificationNum(localId);
//        }
//
//        return selectedOption;
//
//    }

}
