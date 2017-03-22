/*
Modified code from: https://github.com/fpetitjean/DataGeneratorBN/blob/master/src/generator/RandomStructureGenerator.java
*/
package main.generator.componets;


import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.*;

import explorer.ChordalysisModelling;
import lattice.Lattice;
import lattice.LatticeNode;
import model.DecomposableModel;
import org.apache.commons.math3.random.JDKRandomGenerator;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.util.FastMath;
import org.eclipse.recommenders.jayes.BayesNet;
import org.eclipse.recommenders.jayes.BayesNode;
import org.eclipse.recommenders.jayes.util.BayesNodeUtil;
import org.jgrapht.experimental.dag.DirectedAcyclicGraph;
import org.jgrapht.graph.DefaultEdge;

public class BayesianNetworkGenerator {

	int nVariables;
	int nDataPoints;
	double alphaDirichlet;
	int maxNParents;
	int maxNValuesPerNode;
    int minNValuesPerNode;

	RandomDataGenerator r;
	BayesNet bn;
	public BayesNode[] nodes;
	//eliminationOrdering
	int[]eo;

	/**
	 * Creates the generator
	 * @param nVariables number of nodes in the BN
	 * @param maxNParents Maximum number of parents per nodes - sampled from U(0,maxNParents) 
	 * @param maxNValuesPerNode Maximum number of outcome per nodes - sampled from U(2,maxNValuesPerNode)
	 * @param nDataPoints Number of samples to generate
	 * @param alphaDirichlet will sample each multinomial in the CPT (ie line) from Dir(alphaDirichlet). The highest value of alpha, the harder it will be to detect the correlations. 
	 * @param seed Seed to pass to one of the random generator
	 */
	public BayesianNetworkGenerator(int nVariables, int maxNParents, int maxNValuesPerNode, int nDataPoints, double alphaDirichlet, long seed) {
		this.nVariables = nVariables;
		this.nDataPoints = nDataPoints;
		this.maxNParents = maxNParents;
		this.maxNValuesPerNode = maxNValuesPerNode;
		this.alphaDirichlet = alphaDirichlet;

		RandomGenerator rg = new JDKRandomGenerator();
		rg.setSeed(seed);
		this.r = new RandomDataGenerator(rg);
		this.bn = null;

        // TODO: Proper way to input min nodes
        minNValuesPerNode = maxNValuesPerNode - 1;
	}

    public BayesianNetworkGenerator(BayesianNetworkGenerator originStructure) {
        this.nVariables = originStructure.nVariables;
        nDataPoints = originStructure.nDataPoints;
        maxNParents = originStructure.maxNParents;
        maxNValuesPerNode = originStructure.maxNValuesPerNode;
        minNValuesPerNode = originStructure.minNValuesPerNode;
        alphaDirichlet = originStructure.alphaDirichlet;

        this.r = originStructure.r;
        this.eo = originStructure.eo;

        // Copy Bayesian Network
        bn = new BayesNet();
		nodes = new BayesNode[originStructure.nodes.length];
		for (int i = 0; i < originStructure.nodes.length; i++) {
			nodes[i] = bn.createNode("n"+(i+1));
			//int nValues = r.nextInt(minNValuesPerNode, maxNValuesPerNode);
            int nValues = maxNValuesPerNode;
			for (int j = 0; j < originStructure.nodes[i].getOutcomeCount(); j++) {
				nodes[i].addOutcome("s"+j);
			}
		}

        for (int i = 1; i < eo.length; i++) {

            BayesNode child = nodes[eo[i]];

            // looking for the parents of node eo[i]
            List<BayesNode> originParents = originStructure.nodes[eo[i]].getParents();
            int[] parentsIDsInEO = new int[originParents.size()];
            for (int j = 0; j < originParents.size(); j++){
                parentsIDsInEO[j] = originParents.get(j).getId();
            }

            ArrayList<BayesNode>parents = new ArrayList<>();
            for(int eoID:parentsIDsInEO){
                BayesNode parent = nodes[eoID];
                parents.add(parent);
            }

            child.setParents(parents);
        }

        // Copy CPTs
        for (int i = 0; i < nodes.length; i++) {
            nodes[i].setProbabilities(originStructure.nodes[i].getProbabilities().clone());
        }
    }

	// Creating bn from a chordalysis model
	public BayesianNetworkGenerator(ChordalysisModelling modeller, String[] variableNames, ArrayList<ArrayList<String>> outcomes) {
        DecomposableModel model = modeller.getModel();

        DirectedAcyclicGraph<Integer, DefaultEdge> bnGraph = new DirectedAcyclicGraph<>(DefaultEdge.class);
		try {
			bnGraph = model.getBayesianNetwork();
		} catch (DirectedAcyclicGraph.CycleFoundException e) {
			e.printStackTrace();
		}

        nodes = new BayesNode[bnGraph.vertexSet().size()];
		bn = new BayesNet();
        Integer[] nodeIDs = bnGraph.vertexSet().toArray(new Integer[bnGraph.vertexSet().size()]);
        Arrays.sort(nodeIDs);
        eo = new int[nodeIDs.length];
        int nodeSetIndex = 0;
        for (Integer i : nodeIDs) {
            eo[nodeSetIndex] = i;
            nodeSetIndex++;
        }
		for (Integer nodeID : nodeIDs) {
			String name = variableNames[nodeID];
			BayesNode node = bn.createNode(name);
			node.addOutcomes(outcomes.get(nodeID).toArray(new String[outcomes.get(nodeID).size()]));
            nodes[nodeID] = node;
		}

		for (Integer nodeID : nodeIDs) {
			BayesNode node = nodes[nodeID];
			ArrayList<BayesNode> parents = new ArrayList<BayesNode>();
			for (DefaultEdge e : bnGraph.edgesOf(nodeID)) {
				if (bnGraph.getEdgeTarget(e).equals(nodeID)) {
					BayesNode oneParent = nodes[bnGraph.getEdgeSource(e)];
					parents.add(oneParent);
				}
			}
			if (!parents.isEmpty()) {
				node.setParents(parents);
			}
		}

        Lattice lattice = modeller.getLattice();
        // Set Probabilities
		for (Integer nodeID : nodeIDs) {
			BayesNode n = nodes[nodeID];
			// System.out.println("setting CPT for "+n.getName());
			List<BayesNode> parents = n.getParents();
			List<BayesNode> parentsAndChild = new ArrayList<>(parents);
			parentsAndChild.add(n);

			int nbParents = parents.size();

			BitSet numbers = new BitSet();
			numbers.set(nodeID);

			int[] sizes = new int[nbParents];
			int nbRowsInCPT = 1;
			for (int i = 0; i < parents.size(); i++) {
				BayesNode parent = parents.get(i);
				numbers.set(parent.getId());
				sizes[i] = parents.get(i).getOutcomeCount();
				nbRowsInCPT *= sizes[i];
			}

			LatticeNode latticeNode = lattice.getNode(numbers);
			Map<Integer, Integer> fromNodeIDToPositionInSortedTable = new HashMap<>();

			Integer[] variablesNumbers = new Integer[numbers.cardinality()];
			int current = 0;
			for (int i = numbers.nextSetBit(0); i >= 0; i = numbers.nextSetBit(i + 1)) {
				variablesNumbers[current] = i;
				current++;
			}
			for (int i = 0; i < variablesNumbers.length; i++) {
				fromNodeIDToPositionInSortedTable.put(variablesNumbers[i], i);
			}

			int[] counts = new int[nbRowsInCPT * n.getOutcomeCount()];
			int[] indexes4lattice = new int[parentsAndChild.size()];
			int[] indexes4Jayes = new int[parentsAndChild.size()];
			// System.out.println(counts.length +" cases");
			// System.out.println("numbers for lattice "+Arrays.toString(variablesNumbers));

			for (int c = 0; c < counts.length; c++) {
				// System.out.println("case "+c);
				int index = c;
				// find indexes
				for (int i = indexes4Jayes.length - 1; i > 0; i--) {
					BayesNode associatedNode = parentsAndChild.get(i);
					int dim = associatedNode.getOutcomeCount();
					indexes4Jayes[i] = index % dim;
					index /= dim;
				}
				indexes4Jayes[0] = index;

				// System.out.println("indexes jayes = "+Arrays.toString(indexes4Jayes));

				for (int i = 0; i < indexes4Jayes.length; i++) {
					BayesNode nodeInPositionI = parentsAndChild.get(i);
					// System.out.println(nodeInPositionI);
					// System.out.println(fromNodeIDToPositionInSortedTable);
					int nodeInPositionIID = nodeInPositionI.getId();
					int indexInSortedTable = fromNodeIDToPositionInSortedTable
							.get(nodeInPositionIID);
					indexes4lattice[indexInSortedTable] = indexes4Jayes[i];
				}

				// System.out.println("indexes lattice = "+Arrays.toString(indexes4lattice));

				int count = latticeNode.getMatrixCell(indexes4lattice);
				counts[c] = count;
			}
			// System.out.println(Arrays.toString(counts));
			// System.out.println("total="+sumAllCounts);

			double mTerm = 0.5;
			double[] probas1D = new double[n.getOutcomeCount() * nbRowsInCPT];
			for (int s = 0; s < probas1D.length; s += n.getOutcomeCount()) {

				double sumOfCounts = 0.0;
				for (int j = 0; j < n.getOutcomeCount(); j++) {
					sumOfCounts += counts[s + j] + mTerm;
				}

				for (int j = 0; j < n.getOutcomeCount(); j++) {
					probas1D[s + j] = (counts[s + j] + mTerm) / sumOfCounts;
				}
			}
			// System.out.println(Arrays.toString(probas1D));
			n.setProbabilities(probas1D);
		}
	}

	public void generateRandomStructure() {
		RandomDataGenerator rdg = new RandomDataGenerator();

		bn = new BayesNet();
		nodes = new BayesNode[nVariables];
		for (int i = 0; i < nVariables; i++) {
			nodes[i] = bn.createNode("n"+(i+1));
			//int nValues = r.nextInt(minNValuesPerNode, maxNValuesPerNode);
			int nValues = maxNValuesPerNode;
			for (int j = 0; j < nValues; j++) {
				nodes[i].addOutcome("s"+j);
			}
		}
		
		eo = rdg.nextPermutation(nVariables, nVariables);

		for (int i = 1; i < eo.length; i++) {

			BayesNode child = nodes[eo[i]];
			// looking for the parents of node eo[i]
			//int nParents = r.nextInt(0, FastMath.min(i, maxNParents));
			int nParents = r.nextInt(0, FastMath.min(i, 10));
			if(nParents>0){
				int[] parentsIDsInEO = r.nextPermutation(i, nParents);

				ArrayList<BayesNode>parents = new ArrayList<>();
				for(int eoID:parentsIDsInEO){
					BayesNode parent = nodes[eo[eoID]];
					parents.add(parent);
				}
				child.setParents(parents);
			}
		}
		System.out.println("Generated BN structure as follows:");
		System.out.println(toString());
	}

	public String toString() {
		String str = "BN with "+nVariables+" nodes\nnode_name:\n"
				+ "\toutcome_1,...,outcome_v\n"
				+ "\tparent_1,...,parent_k\n";
		
		for (int eoIndex:eo) {
			BayesNode child = nodes[eo[eoIndex]];
			str+=child.getName()+":\n";
			str+="\toutcomes: "+child.getOutcomes().toString()+"\n";
			str+="\tparents: ";
			List<BayesNode> parents = child.getParents();
			if(parents.size()>0){
				str+=parents.toString()+"\n";
			}else{
				str+="no parents\n";
			}
		}
		return str;
	}

	public void generateRandomCPTs(BayesNet net, double hardness){
		List<BayesNode> allNodes = net.getNodes();
		// ~ learning the CPT of node i
		for (int i = 0; i < net.getNodes().size(); i++) {
			BayesNode n = allNodes.get(i);
			List<BayesNode> parents = n.getParents();
			int nbParents = parents.size();

			// ~ CPT = cartesian product (sizes multiplied)
			int[] sizes = new int[nbParents];
			int nbRowsInCPT = 1;
			for (int parent = 0; parent < parents.size(); parent++) {
				sizes[parent] = parents.get(parent).getOutcomeCount();
				nbRowsInCPT *= sizes[parent];
			}

			double[][] probas = new double[nbRowsInCPT][n.getOutcomeCount()];

			// random CPTs
			for (int row = 0; row < nbRowsInCPT; row++) {
				
				//sample from Dirichlet(alphaDirichlet,...,alphaDirichlet)
				//to assign probabilities of multinomial
				double sum = 0.0;
				for (int s = 0; s < n.getOutcomeCount(); s++) {
					probas[row][s]=this.r.nextGamma(hardness, 1.0);
					sum+=probas[row][s];
				}
				for (int s = 0; s < n.getOutcomeCount(); s++) {
					probas[row][s]/=sum;
				}
			}
			double[] probas1D = new double[nbRowsInCPT * n.getOutcomeCount()];
			for (int row = 0; row < nbRowsInCPT; row++) {
				int offset = n.getOutcomeCount() * row;
				for (int s = 0; s < n.getOutcomeCount(); s++) {
					int index = offset+s;
					probas1D[index] = probas[row][s];
				}
			}
//			System.out.println(Arrays.toString(probas1D));
			n.setProbabilities(probas1D);
		}
	}

    public void prepareForUse() {
        if(bn==null){
            generateRandomStructure();
        }
        generateRandomCPTs(bn, alphaDirichlet);
    }

    public int[] getSingleDataPoint() {
        int[] dataPoint = new int[eo.length];
        BayesNode n;
        HashMap<BayesNode, String> evidence = new HashMap<>();
		try {
			SecureRandom secureRandomGenerator = SecureRandom.getInstance("SHA1PRNG");
            for (int i = 0; i < eo.length; i++) {
                n = nodes[eo[i]];
                double[]probas = BayesNodeUtil.getSubCpt(n, evidence);
                double rand = secureRandomGenerator.nextDouble();
                int chosenValue = 0;
                double sumProba = probas[chosenValue];
                while (rand > sumProba) {
                    chosenValue++;
                    assert(chosenValue<probas.length);
                    sumProba += probas[chosenValue];
                }
                dataPoint[eo[i]] = chosenValue;
                String outcome = n.getOutcomeName(chosenValue);
                evidence.put(n, outcome);
            }
		}
		catch (NoSuchAlgorithmException ex) {
			System.out.println("error");
		}
        return dataPoint;
    }

    public double getJointProbabilityOfX(int[] dataVector) {
		assert dataVector.length == nodes.length;

		double jointProb = 1.0f;
		BayesNode n;
		HashMap<BayesNode, String> evidence = new HashMap<>();

        // Add all evidence
		for (int i = 0; i < eo.length; i++) {
			n = nodes[eo[i]];
			int chosenValue = dataVector[eo[i]];
			String outcome = n.getOutcomeName(chosenValue);
			evidence.put(n, outcome);
		}
        // Get joint probability
        for (int i = 0; i < dataVector.length; i++) {
            n = nodes[eo[i]];
            double[] probas = BayesNodeUtil.getSubCpt(n, evidence);
            int chosenValue = dataVector[eo[i]];
            if (chosenValue >= probas.length) return 0.0f;
            jointProb *= probas[chosenValue];
        }
		return jointProb;
	}
}
