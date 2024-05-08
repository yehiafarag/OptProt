/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package no.uib.probe.optprot.util;

import com.compomics.util.experiment.identification.Advocate;
import com.compomics.util.parameters.identification.IdentificationParameters;
import com.compomics.util.parameters.identification.tool_specific.XtandemParameters;
import no.uib.probe.optprot.dataset.model.SearchingSubDataset;

/**
 *
 * @author yfa041
 */
public class ReportExporter {

    public static void addElementToReport(String datasetId, String paramId, String paramOption, double idRate, double timeInSecond) {
        System.out.println("Report --->  datasetId: " + datasetId + "\tparamId:" + paramId + "\tparamOption:" + paramOption + "\tid_rate:" + idRate + "%\ttime:" + timeInSecond);

    }

    public static void printFullReport(IdentificationParameters optimisedSearchParameter, SearchingSubDataset dataset, Advocate searchEngine) {
        System.out.println("-------------------------------------------------------------------------------------------");
        System.out.println("Refrence score      :\t" + dataset.getIdentificationNum());
        System.out.println("Digestion           :\t" + optimisedSearchParameter.getSearchParameters().getDigestionParameters().getCleavageParameter().name());
        System.out.println("Enzyme              :\t" + optimisedSearchParameter.getSearchParameters().getDigestionParameters().getEnzymes().get(0).getName());
        System.out.println("Specificity         :\t" + optimisedSearchParameter.getSearchParameters().getDigestionParameters().getSpecificity(optimisedSearchParameter.getSearchParameters().getDigestionParameters().getEnzymes().get(0).getName()));
        System.out.println("Max Missed Cleavages:\t" + optimisedSearchParameter.getSearchParameters().getDigestionParameters().getnMissedCleavages(optimisedSearchParameter.getSearchParameters().getDigestionParameters().getEnzymes().get(0).getName()));
        System.out.println("Fragment Ion Types  :\t" + optimisedSearchParameter.getSearchParameters().getForwardIons().get(0) + "-" + optimisedSearchParameter.getSearchParameters().getRewindIons().get(0));
        System.out.println("Precursor Accuracy  :\t" + optimisedSearchParameter.getSearchParameters().getPrecursorAccuracy() + " " + optimisedSearchParameter.getSearchParameters().getPrecursorAccuracyType().name());
        System.out.println("Fragment Accuracy   :\t" + optimisedSearchParameter.getSearchParameters().getFragmentIonAccuracy() + " " + optimisedSearchParameter.getSearchParameters().getFragmentAccuracyType().name());
        System.out.println("PrecursorCharge     :\t" + optimisedSearchParameter.getSearchParameters().getMinChargeSearched() + " - " + optimisedSearchParameter.getSearchParameters().getMaxChargeSearched());
        System.out.println("Isotops             :\t" + optimisedSearchParameter.getSearchParameters().getMinIsotopicCorrection() + " - " + optimisedSearchParameter.getSearchParameters().getMaxIsotopicCorrection());
//            System.out.println("default Variable mod:\t" + optimisedSearchParameter.getSearchParameters() + "  Factor " + referenceFactor);
        String fm = "";
        if (optimisedSearchParameter.getSearchParameters().getModificationParameters().getFixedModifications() != null) {
            for (String fixedMod : optimisedSearchParameter.getSearchParameters().getModificationParameters().getFixedModifications()) {
                fm += (fixedMod + "\n");
            }
        }
        System.out.println("Fixed Mod mod:\n" + fm);
        fm = "";
        if (optimisedSearchParameter.getSearchParameters().getModificationParameters().getVariableModifications() != null) {
            for (String v : optimisedSearchParameter.getSearchParameters().getModificationParameters().getVariableModifications()) {
                fm += (v + "% \t" + "\n");
            }
        }
        System.out.println("Variable mod:\n" + fm);
        if (searchEngine.getIndex() == Advocate.xtandem.getIndex()) {
            XtandemParameters xtandemParameters = (XtandemParameters) optimisedSearchParameter.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
            System.out.println("-----------------------------xtandem advanced-------------------------------------------");
            System.out.println("Spectrum Dynamic Range:\t" + xtandemParameters.getDynamicRange());
            System.out.println("Number of Peaks       :\t" + xtandemParameters.getnPeaks());
            System.out.println("MinimumFragmentMz     :\t" + xtandemParameters.getMinFragmentMz());
            System.out.println("Minimum Peaks         :\t" + xtandemParameters.getMinPeaksPerSpectrum());
            System.out.println("Use NoiseSuppression  :\t" + xtandemParameters.isUseNoiseSuppression() + "  (" + xtandemParameters.getMinPrecursorMass() + ")");
            System.out.println("Use Parent isotop exp :\t" + xtandemParameters.getParentMonoisotopicMassIsotopeError());
            System.out.println("Use QuickAcetyl       :\t" + xtandemParameters.isProteinQuickAcetyl());
            System.out.println("Use QuickPyrolidone   :\t" + xtandemParameters.isQuickPyrolidone());
            System.out.println("Use stP Bias          :\t" + xtandemParameters.isStpBias());

            System.out.println("Use Refinement        :\t" + xtandemParameters.isRefine());
            System.out.println("UnanticipatedCleavage :\t" + xtandemParameters.isRefineUnanticipatedCleavages());
            System.out.println("SimiEnzymaticCleavage :\t" + xtandemParameters.isRefineSemi());
            System.out.println("Potintial Modification:\t" + xtandemParameters.isPotentialModificationsForFullRefinment());
            System.out.println("Use PointMutations    :\t" + xtandemParameters.isRefinePointMutations());
            System.out.println("Use SnAPs             :\t" + xtandemParameters.isRefineSnaps());
            System.out.println("Spectrum Synthesis    :\t" + xtandemParameters.isRefineSpectrumSynthesis());

            System.out.println("-------------------------------------------------------------------------------------------");
            String rfm = "";
            if (optimisedSearchParameter.getSearchParameters().getModificationParameters().getRefinementFixedModifications() != null) {
                for (String fixedMod : optimisedSearchParameter.getSearchParameters().getModificationParameters().getRefinementFixedModifications()) {
                    rfm += (fixedMod + "\n");
                }
            }
            System.out.println("Refined Fixed Mod mod:\n" + rfm);
            rfm = "";
            if (optimisedSearchParameter.getSearchParameters().getModificationParameters().getRefinementVariableModifications() != null) {
                for (String v : optimisedSearchParameter.getSearchParameters().getModificationParameters().getRefinementVariableModifications()) {
                    rfm += (v + "\n");
                }
            }
            System.out.println("Refined Variable mod:\n" + rfm);
            System.out.println("-------------------------------------------------------------------------------------------");
        }
    }
}
