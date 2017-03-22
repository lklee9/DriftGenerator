package main.generator;

import main.generator.componets.BayesianNetworkGenerator;
import main.generator.componets.TreeClassGenerator;
import com.yahoo.labs.samoa.instances.*;
import main.generator.drifter.bn.LinearInterpolation;
import main.models.distance.Distance;
import main.models.distance.TotalVariation;
import moa.core.FastVector;
import moa.core.InstanceExample;
import moa.core.ObjectRepository;
import moa.streams.InstanceStream;
import moa.tasks.TaskMonitor;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.optim.linear.LinearConstraint;
import org.apache.commons.math3.random.JDKRandomGenerator;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.eclipse.recommenders.jayes.BayesNode;

import java.math.BigInteger;
import java.util.*;


public class AbruptTreeDriftGenerator extends CategoricalDriftGenerator{
    private InstancesHeader streamHeader;
    private long nInstancesGeneratedSoFar;
    private Distance distanceMetric = new TotalVariation();

    // p(x) before drift
    private BayesianNetworkGenerator bnBD;
    // p(y|x) before drift
    private TreeClassGenerator<String> rootbd;
    // p(x) after drift
    private BayesianNetworkGenerator bnAD;
    // p(y|x) after drift
    private TreeClassGenerator<String> rootad;

    @Override
    public long estimatedRemainingInstances() {
        return -1;
    }

    @Override
    public boolean hasMoreInstances() {
        return true;
    }

    @Override
    public boolean isRestartable() {
        return true;
    }

    @Override
    public void restart() {
	nInstancesGeneratedSoFar = 0L;
    }

    @Override
    public void getDescription(StringBuilder sb, int indent) {}

    @Override
    public String getPurposeString() {
        return "Generates a stream with an abrupt drift of given magnitude.";
    }

    @Override
    public InstancesHeader getHeader() {
        return streamHeader;
    }

    private void generateHeader() {
        FastVector<Attribute> attributes = getHeaderAttributes(nAttributes
                .getValue(), nValuesPerAttribute.getValue());

        this.streamHeader = new InstancesHeader(new Instances(
                getCLICreationString(InstanceStream.class), attributes, 0));
        this.streamHeader.setClassIndex(this.streamHeader.numAttributes() - 1);
    }

    @Override
    public InstanceExample nextInstance() {
        BayesianNetworkGenerator bn = (nInstancesGeneratedSoFar <= burnInNInstances
                .getValue()) ? bnBD : bnAD;
        TreeClassGenerator<String> decisionTee = (nInstancesGeneratedSoFar <= burnInNInstances
                .getValue()) ? rootbd : rootad;

        Instance inst = new DenseInstance(streamHeader.numAttributes());
        inst.setDataset(streamHeader);

        List<String> instanceList = new ArrayList<>();

        // setting values for x_1,...,x_n
        int[] chosenValList = bn.getSingleDataPoint();
        for (int a = 0; a < chosenValList.length; a++) {
            int chosenVal = chosenValList[a];
            inst.setValue(a, chosenVal);
            instanceList.add(Double.toString((double)chosenVal));
        }

        int chosenClassValue = (int)Double.parseDouble(decisionTee.getClassLabel(instanceList));
        inst.setClassValue(chosenClassValue);

        nInstancesGeneratedSoFar++;
        // System.out.println("generated "+inst);
        return new InstanceExample(inst);
    }

    private TreeClassGenerator<String> buildTree() {
        // Get list of values for each attribute and class labels
        List<List<String>> allAttributesValues = new ArrayList<>();
        // Iterate through all possible values for each attributes
        for (int attributeIdx = 0; attributeIdx < nAttributes.getValue(); attributeIdx++) {
            List<String> attributeValues = new ArrayList<>();
            for (int attributeValue = 0; attributeValue < nValuesPerAttribute.getValue(); attributeValue++) {
                attributeValues.add(Double.toString((double)attributeValue));
            }
            allAttributesValues.add(new ArrayList<>(attributeValues));
        }
        // Get list of possible classes
        List<String> classLabels = new ArrayList<>();
        for (int classIdx = 0; classIdx < nValuesPerAttribute.getValue(); classIdx++) {
            classLabels.add(Double.toString((double)classIdx));
        }
        // Start building the tree and return it
        return new TreeClassGenerator<>(allAttributesValues, classLabels, null);
    }

    private TreeClassGenerator<String> driftTree(TreeClassGenerator<String> originalTree) {
        // Gte the number of possible combinations of values
        int nCombinationsValues = 1;
        for (int a = 0; a < nAttributes.getValue(); a++) {
            nCombinationsValues *= nValuesPerAttribute.getValue();
        }
        // Get the number of combinations needed to change in order to reach drift magnitude
        int nCombinationsToChange = (int) (nCombinationsValues * driftMagnitudeConditional.getValue());

        originalTree.changeClassLabel(nCombinationsToChange);
        return originalTree;
    }

    @Override
    protected void prepareForUseImpl(TaskMonitor monitor,
                                     ObjectRepository repository) {
        System.out.println("burnIn=" + burnInNInstances.getValue());
        generateHeader();

        BigInteger nCombinationsValuesForPX = new BigInteger("1");
        for (int a = 0; a < nAttributes.getValue(); a++) {
            nCombinationsValuesForPX = nCombinationsValuesForPX.multiply(new BigInteger(String.valueOf(nValuesPerAttribute.getValue())));
        }

        // generating distribution before drift
        // p(x)
        // Default values
        int nVariables = nAttributes.getValue();
        int maxNParents = nAttributes.getValue() - 1;
        int maxNValues= nValuesPerAttribute.getValue();
        int nDataPoints = 5;
        double alphaDirichlet = 1.0;
        long seed = 3071980L;
        bnBD = new BayesianNetworkGenerator(nVariables, maxNParents,maxNValues,nDataPoints,alphaDirichlet, seed);
        bnBD.prepareForUse();
        // p(y|X)
        System.out.println("Creating Class Tree");
        rootbd = buildTree();
        System.out.println("Done Creating Class Tree");


        // generating distribution after drift
        if (driftPriors.isSet()) {
            System.out.println("Creating p(x) drifter");
            LinearInterpolation drifter = new LinearInterpolation(bnBD);
            bnAD = drifter.getBn_mid();
            double drift  = 0.0f;
            do {
                bnBD.prepareForUse();
                System.out.println("Sampling p(x) for required magnitude...");
                drift = drifter.driftNetwork(driftMagnitudePrior.getValue(), precisionDriftMagnitude.getValue());
                System.out.println("Measuring magnitude...");
            } while (Math.abs(drift - driftMagnitudePrior.getValue()) > precisionDriftMagnitude.getValue());

            System.out.println("exact magnitude for p(x)="
                    + computeMagnitudePX() + "\tasked="
                    + driftMagnitudePrior.getValue());
        } else {
            bnAD = bnBD;
        }

        rootad = driftConditional.isSet() ?
                driftTree(new TreeClassGenerator<>(rootbd, null)) : new TreeClassGenerator<>(rootbd, null);
        double exactDrift = (double)rootbd.nDifferentClasses(rootad) / nCombinationsValuesForPX.doubleValue();
        if (driftConditional.isSet()) System.out.println("exact magnitude for p(y|x)=" + exactDrift
                + "\tasked=" + driftMagnitudeConditional.getValue());

        nInstancesGeneratedSoFar = 0L;
    }

    private double computeMagnitudePX() {
        double driftDist = 0.0f;

        int[][] allCombinations = generateCombinations(0, new int[nAttributes.getValue()]);
        for (int[] combination : allCombinations) {
            driftDist += Math.abs(bnBD.getJointProbabilityOfX(combination) - bnAD.getJointProbabilityOfX(combination));
        }
        driftDist /= 2.0f;

        return driftDist;
    }

    private int[][] generateCombinations(int attributeIndex, int[] auxCombination) {
        if (attributeIndex >= nAttributes.getValue() - 1){
            int[][] combinations = new int[nValuesPerAttribute.getValue()][nAttributes.getValue()];
            for (int i = 0; i < nValuesPerAttribute.getValue(); i++) {
                auxCombination[attributeIndex] = i;
                combinations[i] = auxCombination.clone();
            }
            return combinations;
        }
        else {
            int[][] possibleCombinations = new int[(int)Math.pow(nValuesPerAttribute.getValue(), (nAttributes.getValue()-attributeIndex))][nAttributes.getValue()];
            for (int i = 0; i < nValuesPerAttribute.getValue(); i++) {
                auxCombination[attributeIndex] = i;
                int[][] returnedCombinations = generateCombinations(attributeIndex+1, auxCombination);
                for (int j = 0; j < returnedCombinations.length; j++) {
                    possibleCombinations[(i * returnedCombinations.length) + j] = returnedCombinations[j];
                }
            }
            return possibleCombinations;
        }
    }

    // TODO: Systematically drift Network
    private void driftNetwork(BayesianNetworkGenerator bn) {
        // Total drift so far
        double driftNeeded = driftMagnitudePrior.getValue();
        BayesNode[] nodesSorted = sortNodesDecreasingParents(bn.nodes);
		// for every node in network
		for (int i = 0; i < nodesSorted.length; i++) {
			// Get node to drift in decreasing order of possible number of parents
			BayesNode node = nodesSorted[i];
            // Get parent nodes, 1 dimensional CPT, number of values of current node
            List<BayesNode> parentNodes = node.getParents();

            driftNeeded = driftMagnitudePrior.getValue() - computeMagnitudePX();
            System.out.println("Drifting node " + i + "...");
            if (driftNeeded > precisionDriftMagnitude.getValue()) {
                //System.out.println("Target not reached, drifting by: " + driftNeeded);
                double driftMade = driftCPT(node, new ArrayList<>(parentNodes), new HashMap<>(), driftNeeded);
                //System.out.println("Drift I think I did: " + driftMade);
            }
		}
    }

    private BayesNode[] sortNodesDecreasingParents(BayesNode[] nodes) {
        ArrayList<BayesNode> sortedNodes = new ArrayList<>();
        for (BayesNode node : nodes) {
            boolean inserted = false;
            for (int j = 0; j < sortedNodes.size(); j++) {
                BayesNode sortedNode = sortedNodes.get(j);
                if (sortedNode.getParents().size() > node.getParents().size()) {
                    sortedNodes.add(j, node);
                    inserted = true;
                }
                if (inserted) break;
            }
            if (!inserted) sortedNodes.add(node);
        }
        return sortedNodes.toArray(new BayesNode[sortedNodes.size()]);
    }

    private double driftCPT(BayesNode NodeToDrift, ArrayList<BayesNode> parentNodes,
                            HashMap<Integer, Integer> chosenValuesAux, double targetDrift) {
        // All information needed to describe a row in the CPT
        if (parentNodes.size() == 0) {
            //System.out.println("Drifting Row...");
            double drift = (targetDrift > precisionDriftMagnitude.getValue()) ?
                    driftCPTRow(NodeToDrift, chosenValuesAux, targetDrift) : 0.0f;
        }
        else {
            BayesNode currentNode = parentNodes.remove(0);
            double totalDrift = 0.0f;
            for (int i = 0; i < currentNode.getOutcomeCount(); i++) {
                chosenValuesAux.put(currentNode.getId(), i);
                targetDrift = driftMagnitudePrior.getValue() - computeMagnitudePX();
                totalDrift += (targetDrift > 0.0f) ?
                        driftCPT(NodeToDrift, parentNodes, chosenValuesAux, driftMagnitudePrior.getValue()-computeMagnitudePX()):0;
            }
        }
        return computeMagnitudePX() - (driftMagnitudePrior.getValue() - targetDrift);
    }

    private double driftCPTRow(BayesNode currentNode, HashMap<Integer, Integer> chosenValues, double targetDrift) {
        double transDistLeft = Math.pow(targetDrift * Math.sqrt(2), 2);
        double[] probas1D = currentNode.getProbabilities();
        // Iterate over all possible values 2 at a time
        for (int i = 0; i < currentNode.getOutcomeCount(); i+=2) {
            if (Math.sqrt(transDistLeft) / Math.sqrt(2) > precisionDriftMagnitude.getValue()) {
                // Get the base index in the 1 dimensional CPT
                int baseIndex = 0;
                List<Integer> nodeIndexes = new ArrayList<>(chosenValues.keySet());
                Collections.sort(nodeIndexes);
                for (int nodeIndex : nodeIndexes) {
                    int baseIndexTmp = chosenValues.get(nodeIndex);
                    for (int nodeIndex2 : nodeIndexes) {
                        baseIndexTmp *= (nodeIndex <= nodeIndex2) ? 1 : bnAD.nodes[nodeIndex2].getOutcomeCount();
                    }
                    baseIndex += baseIndexTmp;
                }
                // We got the row index now get the 1D baseIndex
                baseIndex *= currentNode.getOutcomeCount();

                // Add current node with value in current column to the Map
                chosenValues.put(currentNode.getId(), (i) % currentNode.getOutcomeCount());
                int currentColumnIndex = baseIndex + (i % currentNode.getOutcomeCount());
                // Get all the combinations that satisfies all the chosen values
                int[][] relevantCombinations1 = filterCombinations(chosenValues);
                chosenValues.remove(currentNode.getId());   // Remove the current node from the list of nodes to check

                // Add current node with value in next column to the Map
                chosenValues.put(currentNode.getId(), (i + 1) % currentNode.getOutcomeCount());
                int nextColumnIndex = baseIndex + ((i + 1) % currentNode.getOutcomeCount());
                // Get all the combinations that satisfies all the chosen values
                int[][] relevantCombinations2 = filterCombinations(chosenValues);
                chosenValues.remove(currentNode.getId());   // Remove the current node from the list of nodes to check

                double maxEpsilon = getMaxEpsilon(probas1D[currentColumnIndex], probas1D[nextColumnIndex]);
                maxEpsilon = (maxEpsilon < probas1D[currentColumnIndex]) ? maxEpsilon * -1 : maxEpsilon;
                /*
                double epsilon = findAppropriateEpsilon(relevantCombinations1, relevantCombinations2,
                        transDistLeft, maxEpsilon, probas1D[currentColumnIndex]);
                        */
                double epsilon = bruteFindEpsilon(Math.sqrt(transDistLeft) / Math.sqrt(2), maxEpsilon, probas1D,
                        currentColumnIndex, nextColumnIndex);
                /*
                double transDistMade = findTransDistanceMade(relevantCombinations1, epsilon, probas1D[currentColumnIndex])
                        + findTransDistanceMade(relevantCombinations2, epsilon*-1, probas1D[nextColumnIndex]);
                        */
                double transDistMade = Math.pow(computeMagnitudePX() * Math.sqrt(2), 2);
                transDistLeft -= transDistMade;
                transDistLeft = (transDistLeft <= 0.0f) ? 0.0f : transDistLeft;
                /*
                probas1D[currentColumnIndex] += epsilon;
                probas1D[nextColumnIndex] -= epsilon;
                */
            }
        }
        return targetDrift - (Math.sqrt(transDistLeft) / Math.sqrt(2));
    }

    private int[][] filterCombinations(Map<Integer, Integer> attribute_values) {
        int[][] allCombinations = this.generateCombinations(0, new int[nAttributes.getValue()]);
        for (int[] combination : allCombinations) {
            boolean isValid = true;
            for (int j : attribute_values.keySet()) {
                isValid = isValid && (attribute_values.get(j) == combination[j]);
            }
            if (!isValid) ArrayUtils.removeElement(allCombinations, combination);
        }
        return allCombinations;
    }

    private double bruteFindEpsilon(double targetDist, double maxEpsilon, double[] probas1D, int currIdx, int nextIdx) {
        double currOriginalProb = probas1D[currIdx];
        double nextOriginalProb = probas1D[nextIdx];

        double epsilon = (maxEpsilon >= 0.0f) ? maxEpsilon + 0.0001f : maxEpsilon - 0.0001f;
        do {
            epsilon -= (maxEpsilon >= 0.0f) ? 0.0001f : -0.0001f;
            // Change probability of current and next value by epsilon
            probas1D[currIdx] = currOriginalProb + epsilon;
            probas1D[nextIdx] = nextOriginalProb - epsilon;
        } while (computeMagnitudePX() - driftMagnitudePrior.getValue()> precisionDriftMagnitude.getValue() &&
                ((epsilon > 0.0f && Math.abs(maxEpsilon) > currOriginalProb) || (epsilon < 0.0f && Math.abs(maxEpsilon) < currOriginalProb)));
        return (targetDist > 0.0f) ? epsilon : 0.0f;
    }

    private double findAppropriateEpsilon(int[][] relevantComb1, int[][] relevantComb2, double targetTransDrift,
                                          double maxEpsilon, double currentProb) {
        double epsilon = (maxEpsilon >= 0.0f) ? maxEpsilon + 0.001f : maxEpsilon - 0.001f;
        double calculatedTransDist = 0.0f;
        do {
            epsilon -= (maxEpsilon >= 0.0f) ? 0.001f : -0.001f;
            calculatedTransDist = 0.0f;
            for (int[] combination : relevantComb1) {
                double jointProb = bnAD.getJointProbabilityOfX(combination);
                calculatedTransDist += Math.pow(Math.sqrt(jointProb) -
                        Math.sqrt(jointProb + (epsilon*jointProb/currentProb)), 2);
            }
            for (int[] combination : relevantComb2) {
                double jointProb = bnAD.getJointProbabilityOfX(combination);
                calculatedTransDist += Math.pow(Math.sqrt(jointProb) -
                        Math.sqrt(jointProb - (epsilon * jointProb / currentProb)), 2);
            }
            System.out.print("");
        } while (calculatedTransDist - targetTransDrift > 0.0f);
        return (targetTransDrift > 0.0f) ? epsilon : 0.0f;
    }

    private double findTransDistanceMade(int[][] someCombinations, double epsilon, double currentProb) {
        double tranDistance = 0.0f;
        for (int[] combination : someCombinations) {
            double jointProb = bnAD.getJointProbabilityOfX(combination);
            tranDistance += Math.pow(Math.sqrt(jointProb) -
                    Math.sqrt(jointProb + (epsilon*jointProb/currentProb)), 2);
        }
        return tranDistance;
    }

    private double getMaxEpsilon(double current_prob, double next_prob) {
        double epsilon = 1.0f;
        // Reduce epsilon so that it "fits"
        while (((epsilon >= current_prob || epsilon >= 1.0f - next_prob) &&
            (epsilon >= 1.0f - current_prob || epsilon >= next_prob)) && epsilon > 0.0f) {
            epsilon -= 0.001f;
        }
        return epsilon;
    }
}