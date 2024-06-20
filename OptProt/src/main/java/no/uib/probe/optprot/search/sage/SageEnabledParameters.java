/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package no.uib.probe.optprot.search.sage;

import no.uib.probe.optprot.search.myrimatch.*;
import no.uib.probe.optprot.configurations.SearchEngineParameterConfigurations;

/**
 *
 * @author yfa041
 */
public class SageEnabledParameters {

    private final SearchEngineParameterConfigurations paramsToOptimize;

    public SageEnabledParameters() {
        this.paramsToOptimize = new SearchEngineParameterConfigurations();
//        this.paramsToOptimize.disableFragmentTypes();
        this.paramsToOptimize.disableSpecificCTermOnlyEnzyme();
//        this.paramsToOptimize.disableMinCharge();
//        this.paramsToOptimize.disableMinIsotop();
//        this.paramsToOptimize.disableMaxIsotop();
        this.paramsToOptimize.disableParam("Acetylation of peptide N-term");
        this.paramsToOptimize.disableParam("specificNTermOnly");
    }

    public SearchEngineParameterConfigurations getParamsToOptimize() {
        return paramsToOptimize;
    }

}
