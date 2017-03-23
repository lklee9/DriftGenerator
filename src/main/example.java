package main;

import main.generator.AbruptTreeDriftGenerator;
import main.generator.CategoricalDriftGenerator;

/**
 * Created by LoongKuan on 23/03/2017.
 */
public class example {

    public static void main(String[] args) {
        int nBurnIn = 10;
        // Create new generator object
        AbruptTreeDriftGenerator dataStream = new AbruptTreeDriftGenerator();
        // Configure the generator's parameters
        configureDataSet(dataStream, nBurnIn);
        // Initialise generator with configured parameters to produce data stream
        dataStream.prepareForUse();
        for (int i = 0; i < nBurnIn*2; i++) {
            System.out.println(dataStream.nextInstance().getData().toString());
        }
    }

    private static void configureDataSet(CategoricalDriftGenerator dataStream, int nBurnIn) {
        // Number of attributes
        dataStream.nAttributes.setValue(5);
        // Number of values for each attributes
        dataStream.nValuesPerAttribute.setValue(2);
        // Error threshold for drift generated
        dataStream.precisionDriftMagnitude.setValue(0.01);
        // Drift Coavariate?
        dataStream.driftPriors.setValue(true);
        // Drift Posterior?
        dataStream.driftConditional.setValue(true);
        // Covariate Drift Magnitude
        dataStream.driftMagnitudePrior.setValue(0.5);
        // Posterior Drift Magnitude
        dataStream.driftMagnitudeConditional.setValue(0.9);
        // Number of instances before drift
        dataStream.burnInNInstances.setValue(nBurnIn);
    }
}
