package main.generator.componets;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Created by Lee on 4/02/2016.
 **/
public class TreeClassGenerator<T> {
    protected List<T> possibleValues = new ArrayList<>();
    private List<TreeClassGenerator<T>> childNodes = new ArrayList<>();
    private List<T> classLabels = new ArrayList<>();
    private TreeClassGenerator<T> parentNode = null;
    private int classLabelIndex = -1;
    protected int attributeIndex = 0;
    private int weight = 1;

    public TreeClassGenerator(List<List<T>> remainingAttributesValues, List<T> possibleClasses, TreeClassGenerator<T> parent) {
        // Store parent of this node and depth
        parentNode = parent;
        classLabels = possibleClasses;
        attributeIndex = (parentNode == null)? 0 : parentNode.attributeIndex + 1;

        // If there is still attributes left, create new branch
        if (remainingAttributesValues.size() > 0) {
            possibleValues = remainingAttributesValues.remove(0);

            // Set the first attribute value to a terminating node
            childNodes.add(new TreeClassGenerator<>(new ArrayList<>(), classLabels, this));
            for (int i = 1; i < possibleValues.size(); i++) {
                childNodes.add(new TreeClassGenerator<>(new ArrayList<>(remainingAttributesValues), classLabels, this));
            }
        }

        // Choose a fixed class label using a RNG
        Random rg = new Random();
        classLabelIndex = (rg.nextInt(classLabels.size())) % classLabels.size();

        // Assign weight based on how big the impact would be if the class for this attribute were to change
        for (List<T> possibleValues : remainingAttributesValues) {
            weight = weight * possibleValues.size();
        }
    }

    // Copy method
    public TreeClassGenerator(TreeClassGenerator<T> originalTree, TreeClassGenerator<T> parent) {
        possibleValues = new ArrayList<>(originalTree.possibleValues);
        classLabels = new ArrayList<>(originalTree.classLabels);
        parentNode = parent;
        classLabelIndex = originalTree.classLabelIndex;
        attributeIndex = originalTree.attributeIndex;
        weight = originalTree.weight;
        childNodes = new ArrayList<>();
        for (TreeClassGenerator<T> child : originalTree.childNodes) {
            childNodes.add(new TreeClassGenerator<>(child, this));
        }
    }

    /**
     * Check if the current node is a leaf
     * <p>
     *     This method checks if the current node is a leaf.
     *     This is done by checking if the node does not have any remaining possible values.
     * </p>
     * @return
     */
    private boolean isLeaf() {
        if (this.possibleValues.size() == 0) {
            return true;
        }
        else {
            return false;
        }
    }

    /**
     * Get the class label of the instance given by recursively traversing the tree
     * @param instance List of attribute values that represents the remaining values of a instance
     * @return The class label of the instance given
     */
    public T getClassLabel(List<T> instance) {
        // Get the child node based on value passed for this attribute in instance
        T curAttributeValue = instance.remove(0);
        int curAttributeValueIdx = possibleValues.indexOf(curAttributeValue);
        TreeClassGenerator<T> child = childNodes.get(curAttributeValueIdx);
        // If child node is a leaf return the chosen classLabel
        if (child.isLeaf()) {
            return classLabels.get(child.classLabelIndex);
        }
        else {
            // Ask for the class label of the child node
            return child.getClassLabel(instance);
        }
    }

    /**
     * Change the class labels for a number of combinations by changing the class index of leaves in the tree
     * @param totalCombinationsToChange The number of vector combinations' classes to change in the tree
     * @return The total number of combinations whose class labels have been changed
     */
    public int changeClassLabel(int totalCombinationsToChange) {
        int combinationsLeftToChange = totalCombinationsToChange;
        // Check if the current weight is smaller than the weight needed and we are at a terminating node
        if (this.isLeaf() && parentNode.weight <= totalCombinationsToChange) {
            classLabelIndex = (classLabelIndex + 1) % classLabels.size();
            combinationsLeftToChange = combinationsLeftToChange - parentNode.weight;
        }
        else {
            // Distribute weight left to child nodes
            for (TreeClassGenerator<T> node : childNodes) {
                combinationsLeftToChange = combinationsLeftToChange - node.changeClassLabel(combinationsLeftToChange);
            }
        }
        return totalCombinationsToChange - combinationsLeftToChange;
    }

    /**
     * Get the number of identical vector combinations between
     * this tree and the treeToCompare that have different class labels
     * @param treeToCompare The tree to compare with
     * @return The number of vector combinations with different class labels between the 2 trees
     */
    public int nDifferentClasses(TreeClassGenerator<T> treeToCompare) {
        assert this.possibleValues.size() == treeToCompare.possibleValues.size();
        if (this.possibleValues.size() == 0 && treeToCompare.possibleValues.size() == 0) {
            int difference = (this.classLabelIndex == treeToCompare.classLabelIndex) ? 0 : 1;
            return difference * parentNode.weight;
        }
        else {
            int nDifferences = 0;
            for (int i = 0; i < this.possibleValues.size(); i++) {
                TreeClassGenerator<T> child1 = this.childNodes.get(i);
                TreeClassGenerator<T> child2 = treeToCompare.childNodes.get(i);
                nDifferences += child1.nDifferentClasses(child2);
            }
            return nDifferences;
        }
    }

    /**
     * Get the string representation of this node which is the node's possible values
     * @return The node's list of possible values
     */
    public String strRep(){
        return possibleValues.toString();
    }

    public void printNodes() {
        System.out.println(attributeIndex + " : " + strRep());
        for (TreeClassGenerator<T> node : childNodes) {
            node.printNodes();
        }
    }
}
