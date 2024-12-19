package no.uib.probe.optprot.search;

import com.compomics.util.experiment.biology.enzymes.Enzyme;
import com.compomics.util.experiment.biology.enzymes.EnzymeFactory;
import com.compomics.util.experiment.biology.ions.impl.PeptideFragmentIon;
import com.compomics.util.experiment.biology.modifications.Modification;
import com.compomics.util.experiment.biology.modifications.ModificationCategory;
import com.compomics.util.experiment.biology.modifications.ModificationFactory;
import com.compomics.util.experiment.identification.Advocate;
import com.compomics.util.experiment.identification.matches.SpectrumMatch;
import com.compomics.util.io.IoUtil;
import com.compomics.util.parameters.identification.IdentificationParameters;
import com.compomics.util.parameters.identification.search.DigestionParameters;
import com.compomics.util.parameters.identification.search.SearchParameters;
import com.compomics.util.parameters.identification.tool_specific.DirecTagParameters;
import com.compomics.util.parameters.identification.tool_specific.NovorParameters;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import no.uib.probe.optprot.configurations.Configurations;
import no.uib.probe.optprot.dataset.model.SearchingSubDataset;
import no.uib.probe.optprot.model.ParameterScoreModel;
import no.uib.probe.optprot.model.RawScoreModel;
import no.uib.probe.optprot.model.SearchInputSetting;
import no.uib.probe.optprot.model.SortedPTMs;
import no.uib.probe.optprot.util.MainUtilities;
import no.uib.probe.optprot.util.SpectraUtilities;
import org.apache.batik.svggen.font.table.Table;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

/**
 *
 * @author yfa041
 */
public abstract class DefaultOptProtSearchOptimizer {

    /**
     * The compomics PTM factory.
     */
    private final ModificationFactory ptmFactory = ModificationFactory.getInstance();

    public String optimizeDigestionCleavageParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("CleavageParameter");
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        String selectedOption = oreginaltempIdParam.getSearchParameters().getDigestionParameters().getCleavageParameter().name();
        int idRate = optProtDataset.getActiveIdentificationNum();
        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        RawScoreModel oreginalScore = new RawScoreModel("CleavageParameter");
        oreginalScore.setIdPSMNumber(idRate);
        String[] cleavageParameters = new String[]{"wholeProtein", "unSpecific"};

        resultsMap.put(selectedOption, oreginalScore);
        int spectraCounter = optProtDataset.getActiveIdentificationNum();
        for (String cleavageParameter : cleavageParameters) {
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            tempIdParam.getSearchParameters().getDigestionParameters().setCleavageParameter(DigestionParameters.CleavageParameter.valueOf(cleavageParameter));
            if (cleavageParameter.equalsIgnoreCase("enzyme")) {
                tempIdParam.getSearchParameters().getDigestionParameters().addEnzyme(EnzymeFactory.getInstance().getEnzyme("Trypsin"));
            } else {
                tempIdParam.getSearchParameters().getDigestionParameters().clearEnzymes();
            }
            final String option = cleavageParameter;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, true, optimisedSearchParameter, identificationParametersFile, false);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.getFinalScore() > 1) {
                    spectraCounter = Math.max(spectraCounter, scoreModel.getSpectrumMatchResult().size());
                    resultsMap.put(option, f.get());

                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }

        int total = optProtDataset.getActiveIdentificationNum();
        if (!resultsMap.isEmpty()) {
            String bestScore = SpectraUtilities.compareScoresSet(resultsMap);
            selectedOption = bestScore;
            if (selectedOption.equalsIgnoreCase("unSpecific")) {
                double impact = Math.round((double) (resultsMap.get(selectedOption).getSpectrumMatchResult().size() - optProtDataset.getActiveIdentificationNum()) * 100.0 / (double) optProtDataset.getActiveIdentificationNum()) / 100.0;
                paramScore.setImpact(impact);
                total = resultsMap.get(selectedOption).getSpectrumMatchResult().size();
                paramScore.setComments("Extremely slow processing");
            }
        }
        paramScore.setScore(total);
        paramScore.setParamValue(selectedOption);
        parameterScoreSet.add(paramScore);
        return selectedOption;

    }

    /**
     * Optimize digestion enzyme
     *
     * @param optProtDataset
     * @param identificationParametersFile
     * @param optimisedSearchParameter
     * @param parameterScoreSet
     * @return
     * @throws IOException
     */
    public String[] optimizeEnzymeParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("Enzyme");
        String[] values = new String[3];
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//        boolean optimizeOnlyMissedCleaveNum = false;

        values[0] = oreginaltempIdParam.getSearchParameters().getDigestionParameters().getEnzymes().get(0).getName();
        values[1] = oreginaltempIdParam.getSearchParameters().getDigestionParameters().getSpecificity(values[0]).name();
        int missedClavageNumb = oreginaltempIdParam.getSearchParameters().getDigestionParameters().getnMissedCleavages(values[0]);
        values[2] = "" + missedClavageNumb;
//        if(optimisedSearchParameter.isOptimizeEnzymeParameter()){
//         oreginaltempIdParam.getSearchParameters().getDigestionParameters().setCleavageParameter(DigestionParameters.CleavageParameter.valueOf("enzyme"));
//            compareBetweenEnzymes = true;
//            values[1] = "specific";
//        }
//       else if ((optimisedSearchParameter.getDigestionParameterOpt().equalsIgnoreCase("enzyme")) || (optimisedSearchParameter.isOptimizeMaxMissCleavagesParameter()&&!optimisedSearchParameter.isOptimizeEnzymeParameter())) {//            
//            values[0] = oreginaltempIdParam.getSearchParameters().getDigestionParameters().getEnzymes().get(0).getName();
//            values[1] = oreginaltempIdParam.getSearchParameters().getDigestionParameters().getSpecificity(values[0]).name();
//            missedClavageNumb = oreginaltempIdParam.getSearchParameters().getDigestionParameters().getnMissedCleavages(values[0]);
//            System.out.println("max missed cleave number is " + missedClavageNumb);
//            optimizeOnlyMissedCleaveNum = true;
//        } else {
//            oreginaltempIdParam.getSearchParameters().getDigestionParameters().setCleavageParameter(DigestionParameters.CleavageParameter.valueOf("enzyme"));
//            compareBetweenEnzymes = true;
//            values[1] = "specific";
//        }
        Map<String, RawScoreModel> resultsMapI = Collections.synchronizedMap(new LinkedHashMap<>());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        double threshold = 1;
        //optimise enzyme  

        //rerernce search with    Chymotrypsin (no P rule)   delete after get the data
//        {
//            Enzyme enzyme = EnzymeFactory.getInstance().getEnzyme("CNBr");
//            oreginaltempIdParam.getSearchParameters().getDigestionParameters().clearEnzymes();
//            oreginaltempIdParam.getSearchParameters().getDigestionParameters().addEnzyme(enzyme);
//            oreginaltempIdParam.getSearchParameters().getDigestionParameters().setnMissedCleavages(enzyme.getName(), missedClavageNumb);
//            final String option = enzyme.getName();
//            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
//            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
//                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, true, optimisedSearchParameter, identificationParametersFile, false);
//                return scoreModel;
//            });
//            try {
//                RawScoreModel scoreModel = f.get();
//                System.out.println("Enzyme: " + enzyme.getName() + "    Score: " + scoreModel.getFinalScore() + "    #PSMs: " + scoreModel.getIdPSMNumber() + "   " + scoreModel + "   ");
//                optProtDataset.setActiveScoreModel(scoreModel);
//            } catch (InterruptedException | ExecutionException ex) {
//                ex.printStackTrace();
//            }
//        }

        if (optimisedSearchParameter.isOptimizeEnzymeParameter()) {
            for (Enzyme enzyme : EnzymeFactory.getInstance().getEnzymes()) {
                if (enzyme.getName().replace(" ", "").equalsIgnoreCase("Trypsin(noPrule)")) {
                    continue;
                }
                oreginaltempIdParam.getSearchParameters().getDigestionParameters().clearEnzymes();
                oreginaltempIdParam.getSearchParameters().getDigestionParameters().addEnzyme(enzyme);
                oreginaltempIdParam.getSearchParameters().getDigestionParameters().setnMissedCleavages(enzyme.getName(), missedClavageNumb);
                final String option = enzyme.getName();
                final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
                Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                    RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, true, optimisedSearchParameter, identificationParametersFile, false);
                    return scoreModel;
                });
                try {
                    RawScoreModel scoreModel = f.get();
                    System.out.println("Enzyme: " + enzyme.getName() + "    Score: " + scoreModel.getFinalScore() + "    #PSMs: " + scoreModel.getIdPSMNumber() + "   " + scoreModel + "   ");

                    if (scoreModel.isSensitiveChange()) {
                        resultsMapI.put(option, scoreModel);
                    }
                } catch (InterruptedException | ExecutionException ex) {
                    ex.printStackTrace();
                }
            }
            if (!resultsMapI.isEmpty()) {

                String enzymeName = SpectraUtilities.compareScoresSet(resultsMapI);
                values[0] = enzymeName;
                optProtDataset.setActiveScoreModel(resultsMapI.get(enzymeName));
            } else {
                System.out.println("case 3");
//                values[0] = "Trypsin";
            }

            //optimize specifty 
            oreginaltempIdParam.getSearchParameters().getDigestionParameters().clearEnzymes();
            oreginaltempIdParam.getSearchParameters().getDigestionParameters().addEnzyme(EnzymeFactory.getInstance().getEnzyme(values[0]));
            oreginaltempIdParam.getSearchParameters().getDigestionParameters().setnMissedCleavages(values[0], missedClavageNumb);
            oreginaltempIdParam.getSearchParameters().getDigestionParameters().setSpecificity(values[0], DigestionParameters.Specificity.valueOf(values[1]));
        }
        resultsMapI.clear();
        if (optimisedSearchParameter.isOptimizeSpecificityParameter()) {
            for (int i = 0; i < DigestionParameters.Specificity.values().length; i++) {
                final String option = DigestionParameters.Specificity.getSpecificity(i).name();
                oreginaltempIdParam.getSearchParameters().getDigestionParameters().setSpecificity(values[0], DigestionParameters.Specificity.getSpecificity(i));
                final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
                Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                    RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, true, optimisedSearchParameter, identificationParametersFile, false);
                    return scoreModel;
                });
                try {
                    RawScoreModel scoreModel = f.get();
                    if (scoreModel.getFinalScore() > threshold) {
                        resultsMapI.put(option, scoreModel);
                    }

                } catch (ExecutionException | InterruptedException ex) {
                    ex.printStackTrace();
                }
            }
            if (!resultsMapI.isEmpty()) {
                String specifty = SpectraUtilities.compareScoresSet(resultsMapI);
                values[1] = specifty;
                double impact = Math.round((double) (resultsMapI.get(specifty).getSpectrumMatchResult().size() - optProtDataset.getActiveIdentificationNum()) * 100.0 / (double) optProtDataset.getActiveIdentificationNum()) / 100.0;
                paramScore.setImpact(impact);
                optProtDataset.setActiveScoreModel(resultsMapI.get(specifty));

            }
            oreginaltempIdParam.getSearchParameters().getDigestionParameters().setSpecificity(values[0], DigestionParameters.Specificity.valueOf("specific"));
        }

///number op missed cleavage
        resultsMapI.clear();
        if (optimisedSearchParameter.isOptimizeMaxMissCleavagesParameter()) {
            for (int i = 0; i < 5; i++) {
                oreginaltempIdParam.getSearchParameters().getDigestionParameters().setnMissedCleavages(values[0], i);
                final String option = "missedCleavages_" + i;
                final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;

                Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                    RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, true, optimisedSearchParameter, identificationParametersFile, false);
                    return scoreModel;
                });
                try {
                    RawScoreModel scoreModel = f.get();
                    System.out.println("mmissed " + i + "--->" + scoreModel + "    +\\t  " + optProtDataset.getCurrentScoreModel());
//                System.exit(0);
                    if (i < missedClavageNumb && scoreModel.getFinalScore() > 0 && scoreModel.getSharedDataSize() == optProtDataset.getCurrentScoreModel().getIdPSMNumber()) {
                        resultsMapI.put(i + "", scoreModel);
                    } else if (i > missedClavageNumb && scoreModel.getFinalScore() > 0 && scoreModel.getSharedDataSize() == optProtDataset.getCurrentScoreModel().getIdPSMNumber() && scoreModel.getIdPSMNumber() >= 1.05 * optProtDataset.getCurrentScoreModel().getIdPSMNumber()) {
                        resultsMapI.put(i + "", scoreModel);
                    } else if (i > missedClavageNumb) {
                        break;
                    }
                } catch (ExecutionException | InterruptedException ex) {
                    ex.printStackTrace();
                }
            }

            String numbOfMissedCleavage = missedClavageNumb + "";
            if (!resultsMapI.isEmpty()) {
                numbOfMissedCleavage = (SpectraUtilities.compareScoresSet(resultsMapI));
                double impact = Math.round((double) (resultsMapI.get(numbOfMissedCleavage).getSpectrumMatchResult().size() - optProtDataset.getActiveIdentificationNum()) * 100.0 / (double) optProtDataset.getActiveIdentificationNum()) / 100.0;
                paramScore.setImpact(impact);
                optProtDataset.setActiveScoreModel(resultsMapI.get(numbOfMissedCleavage));
            }
            values[2] = numbOfMissedCleavage;
        }
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(Arrays.asList(values).toString());
        parameterScoreSet.add(paramScore);
        return values;

    }

    /**
     *
     * @param optProtDataset
     * @param identificationParametersFile
     * @param optimisedSearchParameter
     * @param parameterScoreSet
     * @return
     * @throws IOException
     */
    public String optimizeFragmentIonTypesParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("FragmentIons");
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        String selectedOption = oreginaltempIdParam.getSearchParameters().getForwardIons() + "-" + oreginaltempIdParam.getSearchParameters().getRewindIons();
        IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        SearchParameters tempSearchParameters = tempIdParam.getSearchParameters();
        ArrayList<Integer> selectedForwardIons = tempSearchParameters.getForwardIons();
        String[] forwardIons = new String[]{"b", "a", "c"};
        String[] rewindIons = new String[]{"y", "x", "z"};

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
                if (option.equalsIgnoreCase(selectedOption)) {
                    continue;
                }

                final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;

                Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                    RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile, false);
                    return scoreModel;
                });
                try {
                    RawScoreModel scoreModel = f.get();
                    if (scoreModel.isAcceptedChange()) {
                        resultsMap.put(option, scoreModel);
                    }
                } catch (ExecutionException | InterruptedException ex) {
                    ex.printStackTrace();
                }
            }
        }

        if (!resultsMap.isEmpty()) {
            String bestScore = SpectraUtilities.compareScoresSet(resultsMap);
            selectedOption = bestScore;
            double impact = Math.round((double) (resultsMap.get(selectedOption + "").getSpectrumMatchResult().size() - optProtDataset.getActiveIdentificationNum()) * 100.0 / (double) optProtDataset.getActiveIdentificationNum()) / 100.0;
            paramScore.setImpact(impact);
            optProtDataset.setActiveScoreModel(resultsMap.get(bestScore));
        }

        selectedOption = selectedOption.replace("[", "").replace("]", "");
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(selectedOption);
        parameterScoreSet.add(paramScore);
        return selectedOption;
    }

    public Integer optimizeMaxMissCleavagesParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("missedCleavages");
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        if (!oreginaltempIdParam.getSearchParameters().getDigestionParameters().getCleavageParameter().name().equalsIgnoreCase("enzyme")) {
            return -1;
        }
        String enzymeName = oreginaltempIdParam.getSearchParameters().getDigestionParameters().getEnzymes().get(0).getName();
        Integer selectedOption = oreginaltempIdParam.getSearchParameters().getDigestionParameters().getnMissedCleavages(enzymeName);
        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());

        for (int i = 0; i < 5; i++) {
            if (i == selectedOption) {
                continue;
            }
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            tempIdParam.getSearchParameters().getDigestionParameters().setnMissedCleavages(enzymeName, i);
            final String option = "missedCleavages_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile, false);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isAcceptedChange()) {
                    resultsMap.put(i + "", scoreModel);
                } else if (i > selectedOption && !scoreModel.isSensitiveChange()) {
                    break;
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }

        if (!resultsMap.isEmpty()) {
            String bestScore = SpectraUtilities.compareScoresSet(resultsMap);
            selectedOption = Integer.valueOf(bestScore);
            double impact = Math.round((double) (resultsMap.get(selectedOption + "").getSpectrumMatchResult().size() - optProtDataset.getActiveIdentificationNum()) * 100.0 / (double) optProtDataset.getActiveIdentificationNum()) / 100.0;
            paramScore.setImpact(impact);
            optProtDataset.setActiveScoreModel(resultsMap.get(bestScore));
        }
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(selectedOption + "");
        parameterScoreSet.add(paramScore);

        return selectedOption;

    }

    public double optimizeFragmentToleranceParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("fragmentAccuracy");
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        double selectedOption = oreginaltempIdParam.getSearchParameters().getFragmentIonAccuracy();
        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        double[] values = new double[]{0.01, 0.02, 0.05, 0.1, 0.2, 0.5};
        double threshold = 1.5 + optProtDataset.getComparisonsThreshold();
        for (double i : values) {
            if (selectedOption == i) {
                continue;
            }
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            tempIdParam.getSearchParameters().setFragmentIonAccuracy(i);
            tempIdParam.getSearchParameters().setFragmentAccuracyType(SearchParameters.MassAccuracyType.DA);
            final String option = "fragmentAccuracy_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile, false);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                System.out.println("Fragment accurcy " + i + "  " + scoreModel + "   first " + (i < selectedOption && scoreModel.isSensitiveChange()) + "   second " + (scoreModel.getFinalScore() > threshold && scoreModel.getIdPSMNumber() >= optProtDataset.getCurrentScoreModel().getIdPSMNumber()) + "   " + threshold);
                if ((i < selectedOption && scoreModel.isSensitiveChange()) || (scoreModel.getFinalScore() > threshold && scoreModel.getIdPSMNumber() >= optProtDataset.getCurrentScoreModel().getIdPSMNumber())) {
                    resultsMap.put(i + "", scoreModel);
                    threshold++;
                } else if (i > selectedOption) {
                    break;
                }

            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }

        }
        if (!resultsMap.isEmpty()) {
            selectedOption = Double.parseDouble(SpectraUtilities.compareScoresSet(resultsMap));
            double impact = Math.round((double) (resultsMap.get(selectedOption + "").getSpectrumMatchResult().size() - optProtDataset.getActiveIdentificationNum()) * 100.0 / (double) optProtDataset.getActiveIdentificationNum()) / 100.0;
            paramScore.setImpact(impact);
            optProtDataset.setActiveScoreModel(resultsMap.get(selectedOption + ""));

        }
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(selectedOption + "");
        parameterScoreSet.add(paramScore);
        return selectedOption;

    }

    public int[] optimizePrecursorChargeParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("charge");
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        int selectedMaxChargeOption = oreginaltempIdParam.getSearchParameters().getMaxChargeSearched();
        int selectedMinChargeOption = oreginaltempIdParam.getSearchParameters().getMinChargeSearched();
        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        int spectraCounter = optProtDataset.getActiveIdentificationNum();
        for (int i = 1; i < 5; i++) {

            for (int j = 2; j <= 5; j++) {
                if (j <= i) {
                    continue;
                }
                tempIdParam.getSearchParameters().setMinChargeSearched(i);
                tempIdParam.getSearchParameters().setMaxChargeSearched(j);

                final String option = "charge-" + i + "," + j;
                final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;

                Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {

                    RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile, false);
                    return scoreModel;
                });
                try {

                    RawScoreModel scoreModel = f.get();
                    if (scoreModel.isAcceptedChange()) {
                        if (scoreModel.getSpectrumMatchResult().size() < spectraCounter) {
                            continue;
                        }
                        spectraCounter = Math.max(spectraCounter, scoreModel.getSpectrumMatchResult().size());
                        resultsMap.put(option, scoreModel);

                    }
                } catch (ExecutionException | InterruptedException ex) {
                    ex.printStackTrace();
                }
            }

        }

        if (!resultsMap.isEmpty()) {
            String bestScore = SpectraUtilities.compareScoresSet(resultsMap);
            double impact = Math.round((double) (resultsMap.get(bestScore + "").getSpectrumMatchResult().size() - optProtDataset.getActiveIdentificationNum()) * 100.0 / (double) optProtDataset.getActiveIdentificationNum()) / 100.0;
            paramScore.setImpact(impact);
            optProtDataset.setActiveScoreModel(resultsMap.get(bestScore));
            String[] topOption = bestScore.split("-")[1].split(",");
            selectedMinChargeOption = Integer.parseInt(topOption[0]);
            selectedMaxChargeOption = Integer.parseInt(topOption[1]);
        }
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(selectedMinChargeOption + "," + selectedMaxChargeOption);
        parameterScoreSet.add(paramScore);
        return new int[]{selectedMinChargeOption, selectedMaxChargeOption};

    }

    public int[] optimizeIsotopParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {

        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("isotop_");
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        int selectedMaxIsotopicCorrectionOption = oreginaltempIdParam.getSearchParameters().getMaxIsotopicCorrection();
        int selectedMinIsotopicCorrectioneOption = oreginaltempIdParam.getSearchParameters().getMinIsotopicCorrection();
        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);

        for (int i = -2; i < 2; i++) {
            for (int j = -1; j <= 2; j++) {
                if (j <= i) {
                    continue;
                }
                tempIdParam.getSearchParameters().setMinIsotopicCorrection(i);
                tempIdParam.getSearchParameters().setMaxIsotopicCorrection(j);
                final String option = "isotop_" + i + "," + j;
                final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;

                Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                    RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile, false);
                    return scoreModel;
                });
                try {
                    RawScoreModel scoreModel = f.get();
                    if (scoreModel.isAcceptedChange()) {
                        resultsMap.put(option, scoreModel);
                    }
                } catch (ExecutionException | InterruptedException ex) {
                    ex.printStackTrace();
                }
            }
        }

        if (!resultsMap.isEmpty()) {
            String bestScore = SpectraUtilities.compareScoresSet(resultsMap);
            double impact = Math.round((double) (resultsMap.get(bestScore + "").getSpectrumMatchResult().size() - optProtDataset.getActiveIdentificationNum()) * 100.0 / (double) optProtDataset.getActiveIdentificationNum()) / 100.0;
            paramScore.setImpact(impact);
            optProtDataset.setActiveScoreModel(resultsMap.get(bestScore));
            String[] topOption = bestScore.split("_")[1].split(",");
            selectedMinIsotopicCorrectioneOption = Integer.parseInt(topOption[0]);
            selectedMaxIsotopicCorrectionOption = Integer.parseInt(topOption[1]);
        }
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(selectedMinIsotopicCorrectioneOption + "," + selectedMaxIsotopicCorrectionOption);
        parameterScoreSet.add(paramScore);

        return new int[]{selectedMinIsotopicCorrectioneOption, selectedMaxIsotopicCorrectionOption};
    }

    public double optimizePrecursorToleranceParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("PrecursorAccuracy");
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        double selectedOption = oreginaltempIdParam.getSearchParameters().getPrecursorAccuracy();
        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        double[] iValues = new double[]{5, 10, 15, 20, 25};
        boolean toEnd = false;
        double threshold = 1.5;
        for (double i : iValues) {
            if (i == selectedOption) {
                continue;
            }
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            tempIdParam.getSearchParameters().setPrecursorAccuracy(i);
            tempIdParam.getSearchParameters().setPrecursorAccuracyType(SearchParameters.MassAccuracyType.PPM);
            final String option = "precursorAccuracy_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile, false);
                return scoreModel;
            });
            try {

                RawScoreModel scoreModel = f.get();
                System.out.println("PT " + i + "  vs " + selectedOption + "  " + scoreModel);
                if (scoreModel.getFinalScore() > threshold) {
                    threshold += 0.5;
                    resultsMap.put(i + "", scoreModel);
                } else if (i > selectedOption) {
                    toEnd = true;
                    break;

                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }

        }

        if (!toEnd) {
            iValues = new double[]{0.1, 0.3, 0.5, 0.7, 0.9};
            if (!optProtDataset.isHighResolutionMassSpectrometers()) {
                for (double i : iValues) {
                    IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
                    tempIdParam.getSearchParameters().setPrecursorAccuracy(i);
                    tempIdParam.getSearchParameters().setPrecursorAccuracyType(SearchParameters.MassAccuracyType.DA);
                    final String option = "precursorAccuracy_Da" + i;
                    final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
                    Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                        RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile, false);
                        return scoreModel;
                    });
                    try {
                        RawScoreModel scoreModel = f.get();
                        if (scoreModel.getFinalScore() > threshold) {
                            threshold += 0.5;
                            resultsMap.put(i + "", scoreModel);
                        }
                    } catch (ExecutionException | InterruptedException ex) {
                        ex.printStackTrace();
                    }

                }
            }
        }
        if (!resultsMap.isEmpty()) {
            String bestScore = SpectraUtilities.compareScoresSet(resultsMap);
            double impact = Math.round((double) (resultsMap.get(bestScore + "").getSpectrumMatchResult().size() - optProtDataset.getActiveIdentificationNum()) * 100.0 / (double) optProtDataset.getActiveIdentificationNum()) / 100.0;
            paramScore.setImpact(impact);
            optProtDataset.setActiveScoreModel(resultsMap.get(bestScore));
            selectedOption = Double.parseDouble(bestScore);
        }

        paramScore.setScore(optProtDataset.getActiveIdentificationNum());

        parameterScoreSet.add(paramScore);
        if (selectedOption >= 5) {
            paramScore.setParamValue(selectedOption + "PPM");
            paramScore.setComments("High-Resolution Mass Spectrometers: Instruments like Orbitrap or Fourier Transform Ion Cyclotron Resonance (FT-ICR)");
        } else {
            paramScore.setParamValue(selectedOption + "Da");
            paramScore.setComments("Low-Resolution Mass Spectrometers: Quadrupole and ion trap mass spectrometers have lower mass accuracy");
        }
        return selectedOption;
    }

    /**
     * This function responsible for selection best combination of
     * (fixed,variable,refinement ) modifications-PTMs
     *
     * @param optProtDataset input data-set
     * @param identificationParametersFile identification file
     * @param searchInputSetting
     * @param parameterScoreSet set to store the parameter information object
     * @return modification map
     * @throws IOException
     */
    public Map<String, Set<String>> optimizeModificationsParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting searchInputSetting, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {

        Set<String> preservedMods = new HashSet<>();
        preservedMods.add("Deamidation of N");
        preservedMods.add("Dimethylation of K");
        preservedMods.add("Methylation of K");
        preservedMods.add("Formylation of K");

        final ParameterScoreModel fixedModParamScore = new ParameterScoreModel();
        fixedModParamScore.setParamId("FixedModifications");

        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        ArrayList<String> selectedFixedModificationOption = new ArrayList<>();
        ArrayList<String> selectedVariableModificationOption = new ArrayList<>();
        Map<String, Set<String>> modificationsResults = new HashMap<>();
        List<String> mods = new ArrayList<>();
        mods.addAll(ptmFactory.getModifications(ModificationCategory.Common));
        mods.addAll(ptmFactory.getModifications(ModificationCategory.Common_Biological));
        mods.addAll(ptmFactory.getModifications(ModificationCategory.Common_Artifact));
        IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);

        tempIdParam.getSearchParameters().getModificationParameters().clearFixedModifications();
        tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
        tempIdParam.getSearchParameters().getModificationParameters().clearRefinementModifications();
        tempIdParam.getSearchParameters().getModificationParameters().getRefinementFixedModifications().clear();
        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        Set<String> potintialMods = new LinkedHashSet<>();
        String commonFixedMod = "Carbamidomethylation of C";
        potintialMods.add(commonFixedMod);
        String commonVariableMod = "Oxidation of M";
        Map<String, RawScoreModel> targtedFixedModificationScore = new TreeMap<>();
        Map<String, RawScoreModel> fullFixedModificationScore = new LinkedHashMap<>();
        //first stage common fixed modification
        String prefix = "f_";
        resultsMap.putAll(this.checkModificationsScores(selectedFixedModificationOption, selectedVariableModificationOption, potintialMods, true, msFileName, tempIdParam, optProtDataset, identificationParametersFile, searchInputSetting, prefix, true));
        if (!resultsMap.isEmpty()) {
            String bestMod = SpectraUtilities.compareScoresSet(resultsMap, true);
            if (resultsMap.get(bestMod).getFinalScore() > 0 || (resultsMap.get(bestMod).getFinalScore() < 0.0 && (resultsMap.get(bestMod).getFinalScore() * -1.0) < 0.1)) {
                selectedFixedModificationOption.add(bestMod);
                optProtDataset.setActiveScoreModel(resultsMap.get(bestMod));
                potintialMods.clear();
                targtedFixedModificationScore.put("C", resultsMap.get(bestMod));
                MainUtilities.cleanOutputFolder();
                resultsMap.clear();
            }

        }
        //try variable coomon mod first
        // stage 2 test for common variable modification   
        final ParameterScoreModel variableModParamScore = new ParameterScoreModel();
        variableModParamScore.setParamId("VariableModifications");
        potintialMods.add(commonVariableMod);
        prefix = "v_";
        resultsMap.putAll(this.checkModificationsScores(selectedFixedModificationOption, selectedVariableModificationOption, potintialMods, false, msFileName, tempIdParam, optProtDataset, identificationParametersFile, searchInputSetting, prefix, false));
        if (!resultsMap.isEmpty()) {
            String bestMod = SpectraUtilities.compareScoresSet(resultsMap, true);
            selectedVariableModificationOption.add(bestMod);
            optProtDataset.setActiveScoreModel(resultsMap.get(bestMod));
            potintialMods.clear();
            MainUtilities.cleanOutputFolder();
            resultsMap.clear();
        }
        MainUtilities.cleanOutputFolder();
        //process variable modifications
        mods.removeAll(selectedFixedModificationOption);
        mods.removeAll(selectedVariableModificationOption);
        //second stage fixed modifications
        for (String modId : mods) {
            String modPattern = ptmFactory.getModification(modId).getPattern().toString();
            if (modPattern.equals("") || ptmFactory.getModification(modId).getModificationType().isNTerm() || ptmFactory.getModification(modId).getModificationType().isCTerm()) {
                modPattern = modPattern + "-" + ptmFactory.getModification(modId).getModificationType().isNTerm() + "-" + ptmFactory.getModification(modId).getModificationType().isCTerm();
            }
            if (targtedFixedModificationScore.containsKey(modPattern)) {
                continue;
            }
            if (!selectedVariableModificationOption.isEmpty() && modPattern.equalsIgnoreCase("M")) {
                continue;
            }
            potintialMods.add(modId);
        }
        prefix = "f_";
        Map<String, Map<String, RawScoreModel>> filterPotintialPtmMap = new LinkedHashMap<>();
        fullFixedModificationScore.putAll(this.checkModificationsScores(selectedFixedModificationOption, selectedVariableModificationOption, potintialMods, true, msFileName, tempIdParam, optProtDataset, identificationParametersFile, searchInputSetting, prefix, true));

        for (String modId : fullFixedModificationScore.keySet()) {
            RawScoreModel scoreModel = fullFixedModificationScore.get(modId);
            System.out.println("Fixed mod test modId: " + modId + " score: " + scoreModel);
            if (scoreModel.isSensitiveChange()) {
                resultsMap.put(modId, scoreModel);

            }
            Modification mod = ptmFactory.getModification(modId);
            //get modified intersection with avctive spectra 
            String modPattern = mod.getPattern().toString();
            if (modPattern.equals("") || ptmFactory.getModification(modId).getModificationType().isNTerm() || ptmFactory.getModification(modId).getModificationType().isCTerm()) {
                modPattern = modPattern + "-" + ptmFactory.getModification(modId).getModificationType().isNTerm() + "-" + ptmFactory.getModification(modId).getModificationType().isCTerm();
            }
            if (!filterPotintialPtmMap.containsKey(modPattern)) {
                filterPotintialPtmMap.put(modPattern, new TreeMap<>(Collections.reverseOrder()));
            }
            filterPotintialPtmMap.get(modPattern).put(modId, scoreModel);
        }

        potintialMods.clear();
        potintialMods.addAll(resultsMap.keySet());
        prefix = "FAV_";
        Map<String, RawScoreModel> clearResults = this.checkModificationsScores(selectedFixedModificationOption, selectedVariableModificationOption, potintialMods, false, msFileName, tempIdParam, optProtDataset, identificationParametersFile, searchInputSetting, prefix, false);
        for (String modId : clearResults.keySet()) {
            RawScoreModel vScore = clearResults.get(modId);
            if (vScore.isSensitiveChange()) {
                if (vScore.getFinalScore() > 1.1 * resultsMap.get(modId).getFinalScore() && vScore.getIdPSMNumber() > 1.02 * resultsMap.get(modId).getIdPSMNumber()) {// (SpectraUtilities.isBetterScore(resultsMap.get(modId).getSpectrumMatchResult(), vScore.getSpectrumMatchResult(), optProtDataset.getTotalSpectraNumber()) > 0) {
                    System.out.println("----------------------remove from fixed " + modId + "-------------------------------------------");
                    resultsMap.remove(modId);
                }
            }
        }
//        System.exit(0);
        //2nd fixed modifications

        if (!resultsMap.isEmpty()) {

            String bestMod = SpectraUtilities.compareScoresSet(resultsMap, true);
            String modPattern = ptmFactory.getModification(bestMod).getPattern().toString();
            if (modPattern.equals("") || ptmFactory.getModification(bestMod).getModificationType().isNTerm() || ptmFactory.getModification(bestMod).getModificationType().isCTerm()) {
                modPattern = modPattern + "-" + ptmFactory.getModification(bestMod).getModificationType().isNTerm() + "-" + ptmFactory.getModification(bestMod).getModificationType().isCTerm();
            }
            if (resultsMap.containsKey("Acetylation of protein N-term") && resultsMap.containsKey("Carbamilation of protein N-term")) {
                //we favor Acetylation of protein N-term in the case of multi n-teminal fixed mod
                if (resultsMap.get("Carbamilation of protein N-term").getFinalScore() <= resultsMap.get("Acetylation of protein N-term").getFinalScore() * 1.05 || resultsMap.get("Carbamilation of protein N-term").getIdPSMNumber() <= resultsMap.get("Acetylation of protein N-term").getIdPSMNumber()) {
                    resultsMap.remove("Carbamilation of protein N-term");
                    if (bestMod.equalsIgnoreCase("Carbamilation of protein N-term")) {
                        bestMod = "Acetylation of protein N-term";
                    }
                }
            }
            selectedFixedModificationOption.add(bestMod);
            optProtDataset.setActiveScoreModel(resultsMap.get(bestMod));
            potintialMods.clear();
            targtedFixedModificationScore.put(modPattern, resultsMap.get(bestMod));
            MainUtilities.cleanOutputFolder();
            resultsMap.remove(bestMod);
        }

        //3rd and 4th  fixed modification processing 
        boolean test = !resultsMap.isEmpty();
        int counter = 3;
        while (test) {
            for (String modId : resultsMap.keySet()) {
                String modPattern = ptmFactory.getModification(modId).getPattern().toString();
                if (modPattern.equals("") || ptmFactory.getModification(modId).getModificationType().isNTerm() || ptmFactory.getModification(modId).getModificationType().isCTerm()) {
                    modPattern = modPattern + "-" + ptmFactory.getModification(modId).getModificationType().isNTerm() + "-" + ptmFactory.getModification(modId).getModificationType().isCTerm();
                }

                if (targtedFixedModificationScore.containsKey(modPattern)) {
                    continue;
                }
                potintialMods.add(modId);
            }
            resultsMap.clear();
            if (!potintialMods.isEmpty()) {
                prefix = "f_";
                resultsMap.putAll(this.checkModificationsScores(selectedFixedModificationOption, selectedVariableModificationOption, potintialMods, true, msFileName, tempIdParam, optProtDataset, identificationParametersFile, searchInputSetting, prefix, false));
                if (!resultsMap.isEmpty()) {
                    String bestMod = SpectraUtilities.compareScoresSet(resultsMap, true);
                    selectedFixedModificationOption.add(bestMod);
                    optProtDataset.setActiveScoreModel(resultsMap.get(bestMod));
                    potintialMods.clear();
                    String modPattern = ptmFactory.getModification(bestMod).getPattern().toString();
                    if (modPattern.equals("") || ptmFactory.getModification(bestMod).getModificationType().isNTerm() || ptmFactory.getModification(bestMod).getModificationType().isCTerm()) {
                        modPattern = modPattern + "-" + ptmFactory.getModification(bestMod).getModificationType().isNTerm() + "-" + ptmFactory.getModification(bestMod).getModificationType().isCTerm();
                    }
                    targtedFixedModificationScore.put(modPattern, resultsMap.get(bestMod));
                    MainUtilities.cleanOutputFolder();
                    System.out.println("targtedFixedModificationScore " + targtedFixedModificationScore.keySet());
                    resultsMap.remove(bestMod);
                    test = !resultsMap.isEmpty() || counter <= 4;
                    counter++;
                } else {
                    test = false;
                }
            } else {
                break;
            }

        }
        modificationsResults.put("fixedModifications", new HashSet<>(selectedFixedModificationOption));
        modificationsResults.put("refinmentFixedModifications", new HashSet<>(selectedFixedModificationOption));
        preservedMods.removeAll(selectedFixedModificationOption);

        MainUtilities.cleanOutputFolder();
        fixedModParamScore.setScore(optProtDataset.getActiveIdentificationNum());
        fixedModParamScore.setParamValue(selectedFixedModificationOption.toString());
        parameterScoreSet.add(fixedModParamScore);
        resultsMap.clear();
        potintialMods.clear();

        for (String fixedMod : selectedFixedModificationOption) {
            String modPattern = ptmFactory.getModification(fixedMod).getPattern().toString();
            if (modPattern.equals("") || ptmFactory.getModification(fixedMod).getModificationType().isNTerm() || ptmFactory.getModification(fixedMod).getModificationType().isCTerm()) {
                modPattern = modPattern + "-" + ptmFactory.getModification(fixedMod).getModificationType().isNTerm() + "-" + ptmFactory.getModification(fixedMod).getModificationType().isCTerm();
            }
            if (filterPotintialPtmMap.containsKey(modPattern)) {
                filterPotintialPtmMap.remove(modPattern);
            }

        }

        for (String modPatternKey : filterPotintialPtmMap.keySet()) {
            if (modPatternKey.contains("-")) {
                for (String ptm : filterPotintialPtmMap.get(modPatternKey).keySet()) {
                    potintialMods.add(ptm);
                }
            } else {
                String bestTargeted = SpectraUtilities.getTopScoresSet(filterPotintialPtmMap.get(modPatternKey), optProtDataset.getCurrentScoreModel().getSpecTitles());
                potintialMods.add(bestTargeted.split("_-_")[0]);
                if (bestTargeted.contains("_-_")) {
                    potintialMods.add(bestTargeted.split("_-_")[1]);
                }
            }
            for (String str : preservedMods) {
                if (filterPotintialPtmMap.get(modPatternKey).containsKey(str)) {
                    potintialMods.add(str);
                }
            }

        }
        resultsMap.clear();
        prefix = "v_";
        counter = 0;
        double thre = 0;
        Map<String, TreeMap<Double, String>> filterVMMap = new LinkedHashMap<>();
        while (counter < 4) {
            resultsMap.putAll(this.checkModificationsScores(selectedFixedModificationOption, selectedVariableModificationOption, potintialMods, false, msFileName, tempIdParam, optProtDataset, identificationParametersFile, searchInputSetting, counter + "" + prefix, false));
            filterVMMap.clear();
            if (!resultsMap.isEmpty()) {
                for (String modId : resultsMap.keySet()) {
                    Modification mod = ptmFactory.getModification(modId);
//                    get modified intersection with avctive spectra 
                    String modPattern = mod.getPattern().toString();
                    if (modPattern.equals("") || mod.getModificationType().isNTerm() || mod.getModificationType().isCTerm()) {
                        modPattern = modPattern + "-" + mod.getModificationType().isNTerm() + "-" + mod.getModificationType().isCTerm();
                    }
                    if (modPattern.equalsIgnoreCase("true-false") && resultsMap.containsKey("Acetylation of protein N-term") && resultsMap.containsKey("Carbamilation of protein N-term")) {
                        if (resultsMap.get(modId).getFinalScore() <= resultsMap.get("Acetylation of protein N-term").getFinalScore() * 1.05 || resultsMap.get(modId).getIdPSMNumber() <= resultsMap.get("Acetylation of protein N-term").getIdPSMNumber()) {
                            modId = "Acetylation of protein N-term";
                        }
                    }

                    if (!filterVMMap.containsKey(modPattern)) {
                        filterVMMap.put(modPattern, new TreeMap<>(Collections.reverseOrder()));
                    }
                    filterVMMap.get(modPattern).put(resultsMap.get(modId).getFinalScore(), modId);

                }
//                System.out.println("filtered --------------->> " + counter + "  vm are " + filterVMMap);
                Set<String> toRemove = new HashSet();
                for (String patteren : filterVMMap.keySet()) {
                    if (patteren.contains("-")) {
                        continue;
                    }
                    TreeMap<Double, String> modPatMap = filterVMMap.get(patteren);
                    Map<String, RawScoreModel> subresultsMap = new LinkedHashMap<>();
                    for (String mod : modPatMap.values()) {
                        subresultsMap.put(mod, resultsMap.get(mod));
                        toRemove.add(mod);
                    }
                    String tokeepMod = SpectraUtilities.compareScoresSet(subresultsMap, true);
                    toRemove.remove(tokeepMod);
                    System.out.println("subResult--->> " + patteren + "--to keep " + toRemove + "  " + tokeepMod + "   " + subresultsMap.keySet());
                }

                for (String remove : toRemove) {
                    resultsMap.remove(remove);
                }
                String bestMod = SpectraUtilities.compareScoresSet(resultsMap, true);

                if (resultsMap.get(bestMod).getFinalScore() > thre) {
                    selectedVariableModificationOption.add(bestMod);
                    optProtDataset.setActiveScoreModel(resultsMap.get(bestMod));
                    MainUtilities.cleanOutputFolder();

//                    if (!optProtDataset.isFullDataSpectaInput()) {
                    System.out.println("before thers is " + thre);
                    thre += 0.01;
                    System.out.println("updated thers is " + thre);
                    resultsMap.remove(bestMod);
//                    } else {
//                        thre += optProtDataset.getComparisonsThreshold();
//                    }
                    System.out.println("updated thre " + thre + "   " + resultsMap.size());
                } else {
                    System.out.println("resultsMap.get(bestMod).getFinalScore()  " + resultsMap.get(bestMod).getFinalScore() + "   " + thre);
                    resultsMap.clear();
                }
                if (resultsMap.isEmpty()) {
                    break;
                }

                TreeSet<SortedPTMs> sorePtms = new TreeSet<>(Collections.reverseOrder());
                for (String mod : resultsMap.keySet()) {
                    if (resultsMap.get(mod).isSensitiveChange()) {
                        sorePtms.add(new SortedPTMs(mod, resultsMap.get(mod).getFinalScore(), 0));

                    }
                }
                resultsMap.clear();
                potintialMods.clear();
                int subCounter = selectedVariableModificationOption.size();
                for (SortedPTMs mod : sorePtms) {
                    if (subCounter > 4) {
                        break;
                    }
                    potintialMods.add(mod.getName());
                    subCounter++;
                    System.out.println(counter + "-->> mod " + mod.getName() + "  " + mod.getScore());

                }
                counter++;

            } else {
                break;
            }
        }
        variableModParamScore.setScore(optProtDataset.getActiveIdentificationNum());
        variableModParamScore.setParamValue(selectedVariableModificationOption.toString());
        parameterScoreSet.add(variableModParamScore);
        modificationsResults.put("variableModifications", new HashSet<>(selectedVariableModificationOption));
        MainUtilities.cleanOutputFolder();
        return modificationsResults;

    }

    private Map<String, RawScoreModel> checkModificationsScores(ArrayList<String> selectedFixedModificationOption, ArrayList<String> selectedVariableModificationOption, Set<String> modifications, boolean fixed, String msFileName, IdentificationParameters tempIdParam, SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting searchInputSetting, String prefix, boolean addAll) {
        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        for (String modId : modifications) {
            final String option = modId;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + prefix + option + "_" + msFileName;
            tempIdParam.getSearchParameters().getModificationParameters().clearFixedModifications();
            tempIdParam.getSearchParameters().getModificationParameters().getRefinementFixedModifications().clear();
            tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
//            tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(modId));
            for (String fixedMod : selectedFixedModificationOption) {
                tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(fixedMod));
                tempIdParam.getSearchParameters().getModificationParameters().addRefinementFixedModification(ptmFactory.getModification(fixedMod));
            }
            for (String variableMod : selectedVariableModificationOption) {
                tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(variableMod));
            }
            if (fixed) {
                tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(modId));
                tempIdParam.getSearchParameters().getModificationParameters().addRefinementFixedModification(ptmFactory.getModification(modId));
            } else {
                tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(modId));
            }
            System.out.println("modification " + modId + "  " + selectedFixedModificationOption + "   " + selectedVariableModificationOption);
            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, true, searchInputSetting, identificationParametersFile, false);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.getFinalScore() > 0 || addAll) {
                    resultsMap.put(modId, scoreModel);//                 
                }

            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }

        }
        tempIdParam.getSearchParameters().getModificationParameters().clearFixedModifications();
        tempIdParam.getSearchParameters().getModificationParameters().getRefinementFixedModifications().clear();
        tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();

        return resultsMap;

    }

    public abstract RawScoreModel excuteSearch(SearchingSubDataset optProtDataset, String defaultOutputFileName, String paramOption, IdentificationParameters tempIdParam, boolean addPeptideMasses, SearchInputSetting optProtSearchSettings, File identificationParametersFile, boolean pairData);

}
