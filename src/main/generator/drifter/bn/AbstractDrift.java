package main.generator.drifter.bn;

import main.generator.componets.BayesianNetworkGenerator;
import main.models.prior.BayesianNetwork;

/**
 * Created by loongkuan on 26/04/16.
 */
public abstract class AbstractDrift {
    BayesianNetworkGenerator bn_base;
    BayesianNetworkGenerator bn_max;
    BayesianNetworkGenerator bn_mid;

    public abstract double driftNetwork(double targetDist, double epsilon);

    public BayesianNetworkGenerator getBn_mid() {
        return bn_mid;
    }

    double computeMagnitudePX(BayesianNetworkGenerator bnBD, BayesianNetworkGenerator bnAD, int[][] allCombinations) {
        double driftDist = 0.0f;

        for (int[] combination : allCombinations) {
            driftDist += Math.abs(bnBD.getJointProbabilityOfX(combination) - bnAD.getJointProbabilityOfX(combination));
        }
        driftDist = driftDist / 2.0f;

        return driftDist;
    }

    int[][] generateCombinations(int attributeIndex, int[] auxCombination) {
        int nValuesPerAttribute = bn_mid.nodes[0].getOutcomeCount();
        int nAttributes = bn_mid.nodes.length;
        if (attributeIndex >= nAttributes - 1){
            int[][] combinations = new int[nValuesPerAttribute][nAttributes];
            for (int i = 0; i < nValuesPerAttribute; i++) {
                auxCombination[attributeIndex] = i;
                combinations[i] = auxCombination.clone();
            }
            return combinations;
        }
        else {
            int[][] possibleCombinations = new int[(int)Math.pow(nValuesPerAttribute, (nAttributes-attributeIndex))][nAttributes];
            for (int i = 0; i < nValuesPerAttribute; i++) {
                auxCombination[attributeIndex] = i;
                int[][] returnedCombinations = generateCombinations(attributeIndex+1, auxCombination);
                for (int j = 0; j < returnedCombinations.length; j++) {
                    possibleCombinations[(i * returnedCombinations.length) + j] = returnedCombinations[j];
                }
            }
            return possibleCombinations;
        }
    }
}
