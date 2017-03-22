package main.generator.drifter.bn;

import main.generator.componets.BayesianNetworkGenerator;
import org.eclipse.recommenders.jayes.BayesNode;

/**
 * Created by Lee on 25/04/2016.
 **/
public class LinearInterpolation extends AbstractDrift{
    private int[][] allCombinations;

    public LinearInterpolation(BayesianNetworkGenerator bn) {
        this.bn_base = bn;
        this.bn_max = generateDistantBN(bn);
        this.bn_mid = new BayesianNetworkGenerator(bn_max);
        System.out.println("Generating all Combinations");
        this.allCombinations = generateCombinations(0, new int[bn_mid.nodes.length]);
    }

    public double driftNetwork(double targetDist, double epsilon) {
        double wMid = 1.0;
        double wMin = 0.0;
        double wMax = 1.0;
        double currentDistance = computeMagnitudePX(this.bn_base, this.bn_mid, allCombinations);
        do {
            if (currentDistance < targetDist) {
                wMin = wMid;
                wMid = wMin + (wMax - wMin)/2;
            }
            else if (currentDistance > targetDist){
                wMax = wMid;
                wMid = wMin + (wMax - wMin)/2;
            }
            else break;
            System.out.println("Target not reached, drifting with " + wMid + "...");
            this.interpolate(wMid);
            currentDistance = computeMagnitudePX(this.bn_base, this.bn_mid, allCombinations);
            System.out.println("Drift Created: " + currentDistance + " Drift Needed: " + targetDist);
        } while (Math.abs(currentDistance - targetDist) > epsilon && wMax > wMin);
        System.out.println("Drift Created: " + currentDistance + " Drift Needed: " + targetDist);
        return currentDistance;
    }

    private BayesianNetworkGenerator generateDistantBN(BayesianNetworkGenerator bn) {
        BayesianNetworkGenerator bn_far = new BayesianNetworkGenerator(bn);
        for (BayesNode node : bn_far.nodes) {
            double[] prob1D = node.getProbabilities();
            for (int i = 0; i < prob1D.length; i++) {
                prob1D[i] = (1 - prob1D[i])/(node.getOutcomeCount() - 1);
            }
            node.setProbabilities(prob1D);
        }
        return bn_far;
    }

    private void interpolate(double w) {
        for (int i = 0; i < this.bn_mid.nodes.length; i++) {
            double[] MidProbs = this.bn_mid.nodes[i].getProbabilities();
            double[] BaseProbs = this.bn_base.nodes[i].getProbabilities();
            double[] MaxProbs = this.bn_max.nodes[i].getProbabilities();
            for (int j = 0; j < MidProbs.length; j++) {
                MidProbs[j] = BaseProbs[j] + (w * (MaxProbs[j] - BaseProbs[j]));
            }
            bn_mid.nodes[i].setProbabilities(MidProbs);
        }
    }
}
