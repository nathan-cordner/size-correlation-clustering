package com.google.ortools;

import com.google.ortools.Loader;
import com.google.ortools.linearsolver.MPConstraint;
import com.google.ortools.linearsolver.MPObjective;
import com.google.ortools.linearsolver.MPSolver;
import com.google.ortools.linearsolver.MPVariable;

import java.util.*;


public class RunNonUniformAlg1 {
    
    public static void main(String args[]) {
        String data_set = args[0]; // "cor_landmarks";
        String delimiter = "\\s"; 
        ArrayList<ArrayList<Integer>> prob_matrix = Helper.read_large_network_relabel("Data/"+ data_set + "/graph.txt", delimiter);
        int num_nodes = prob_matrix.size();
        System.out.println("Data set with " + num_nodes + " nodes");

        Loader.loadNativeLibraries();

        // STANDARD TEST

        int ROUNDS = 10;

        int[] k_vals = {0, 1, 2}; //, 15, 20, 25, 30, 35, 40, 45, 50};

        for (int q = 0; q < k_vals.length; q++) {

        ArrayList<Integer> size_bounds = Helper.read_size_bounds("Data/"+ data_set + "/" + data_set +"_bounds" + k_vals[q] + ".txt");

        int K = k_vals[q];
        System.out.println("Start: k = " + K);

        long lpStart = System.currentTimeMillis();

        MPSolver solver = MPSolver.createSolver("CLP");

        MPVariable[][] var_array = new MPVariable[num_nodes][num_nodes];

        for (int i = 0; i < num_nodes; i++) {
            for (int j = i+1; j < num_nodes; j++) {
                // String name = i + "_" + j;
                var_array[i][j] = solver.makeNumVar(0.0, 1.0, i + "_" + j); // + if i < j
                var_array[j][i] = solver.makeNumVar(0.0, 1.0, j + "_" + i); // - if i > j
            }
        }
    
        System.out.println("Number of variables = " + solver.numVariables());

        // probability constraints
        for (int i = 0; i < num_nodes; i++) {
            for (int j = i + 1; j < num_nodes; j++) {
            MPConstraint ct = solver.makeConstraint(1.0, 1.0, i + "_" + j);
            ct.setCoefficient(var_array[i][j], 1); // assume all other variables are 0?
            ct.setCoefficient(var_array[j][i], 1);
            }
        }

        // size constraints
        for (int i = 0; i < num_nodes; i++) {
            MPConstraint ct = solver.makeConstraint(0.0, size_bounds.get(i), i + ""); // reduce K by 1 for PM notation
            for (int j = 0; j < num_nodes; j++) {
                // count number of nodes in same cluster with i
                if (i < j)
                    ct.setCoefficient(var_array[i][j], 1); // assume all other variables are 0?
                else if (i > j)
                    ct.setCoefficient(var_array[j][i], 1);
            }
        }

        System.out.println("Constraints so far = " + solver.numConstraints());


        // triangle inequality constraints
        for (int i = 0; i < num_nodes; i++) {
            for (int j = i + 1; j < num_nodes; j++) {
                for (int k = j + 1; k < num_nodes; k++) {
                    MPConstraint ct = solver.makeConstraint(0.0, 2.0, i + "_" + j + "_" + k);
                    ct.setCoefficient(var_array[j][i], 1); // assume all other variables are 0?
                    ct.setCoefficient(var_array[k][j], 1);     
                    ct.setCoefficient(var_array[k][i], -1); 
                }
    
            }
        }

        System.out.println("Number of constraints = " + solver.numConstraints());

        MPObjective objective = solver.objective();
        for (int i = 0; i < num_nodes; i++) {
            for (int j = i + 1; j < num_nodes; j++) {
                int prob = 0;
                if (prob_matrix.get(i).contains(j))
                    prob = 1; // edge exists!
                objective.setCoefficient(var_array[i][j], 1 - prob);
                objective.setCoefficient(var_array[j][i], prob);
            }
        }

        objective.setMinimization();
        solver.setTimeLimit(600000); // 10 minute limit
        solver.solve();
        long lpTime = System.currentTimeMillis() - lpStart;

        System.out.println("LP Objective value = " + objective.value());
        System.out.println("Total time: " + (lpTime / 1000.0));
        System.out.println();


        // Collect numbers here
        double[] pivotTimes = new double[ROUNDS];
        long[] pivotScores = new long[ROUNDS];

        double[] blendTimes = new double[ROUNDS];
        long[] blendScores = new long[ROUNDS];

        double[] hybridTimes = new double[ROUNDS];
        long[] hybridScores = new long[ROUNDS];

        double[] lsTimes = new double[ROUNDS];
        long[] lsScores = new long[ROUNDS];


        double pivotTimeTotal = 0;
        double blendTimeTotal = 0; // using "blend" for constrained RNode
        double hybridTimeTotal = 0;
        double lsTimeTotal = 0;

        int pivotNumClusters = 0;
        int blendNumClusters = 0;
        int hybridNumClusters = 0;
        int lsNumClusters = 0;

        int largestPivotCluster = 0;
        int largestBlendCluster = 0;
        int largestHybridCluster = 0;
        int largestLsCluster = 0;

        long pivotScoreTotal = 0;
    	long blendScoreTotal = 0;
	    long hybridScoreTotal = 0;
        long lsScoreTotal = 0;

        int k = k_vals[q]; 

        System.out.println("Start: k = " + k);
        for (int j = 0; j < ROUNDS; j++) {

        long pivotStart = System.currentTimeMillis();
        ArrayList<ArrayList<Integer>> pivot_result = RoundingAlg1(var_array, size_bounds, 0.5);
        long pivotTime = System.currentTimeMillis() - pivotStart;
        pivotTimeTotal += (pivotTime / 1000.0);
        pivotTimes[j] = pivotTime;
        // System.out.println("finish pivot " + j);

        long blendStart = System.currentTimeMillis();
        ArrayList<ArrayList<Integer>> blend_result = RoundingAlg1(var_array, size_bounds, 0.1);
        long blendTime = System.currentTimeMillis() - blendStart;
        blendTimeTotal += (blendTime / 1000.0); 
        blendTimes[j] = blendTime;
        // System.out.println("finish rnode " + j);

        long hybridStart = System.currentTimeMillis();
        ArrayList<ArrayList<Integer>> fix = RoundingAlg1(var_array, size_bounds, 0.3);
        long hybridTime = System.currentTimeMillis() - hybridStart;
        hybridTimeTotal += (hybridTime / 1000.0);
        hybridTimes[j] = hybridTime;
        // System.out.println("finish hybrid " + j);

        long lsStart = System.currentTimeMillis();
        ArrayList<ArrayList<Integer>> ls_result = RoundingAlg2(var_array, size_bounds);
        long lsTime = System.currentTimeMillis() - lsStart;
        lsTimeTotal += (lsTime / 1000.0);
        lsTimes[j] = lsTime;

        pivotNumClusters += pivot_result.size();
        blendNumClusters += blend_result.size();
        hybridNumClusters += fix.size();
        lsNumClusters += ls_result.size();

        int max_pivot = 0;
        for (int i = 0; i < pivot_result.size(); i++) {
            if (pivot_result.get(i).size() > max_pivot)
                max_pivot = pivot_result.get(i).size();
        }
        largestPivotCluster += max_pivot;

        int max_blend = 0;
        for (int i = 0; i < blend_result.size(); i++) {
            if (blend_result.get(i).size() > max_blend)
                max_blend = blend_result.get(i).size();
        }
        largestBlendCluster += max_blend;
        
        int max_hybrid = 0;
        for (int i = 0; i < fix.size(); i++) {
            if (fix.get(i).size() > max_hybrid)
                max_hybrid = fix.get(i).size();
        }
        largestHybridCluster += max_hybrid;
        
        int max_ls = 0;
        for (int i = 0; i < ls_result.size(); i++) {
            if (ls_result.get(i).size() > max_ls)
                max_ls = ls_result.get(i).size();
        }
        largestLsCluster += max_ls;

        long pivot_score = Helper.quick_edit_dist(pivot_result, prob_matrix);
        pivotScoreTotal += pivot_score;
        pivotScores[j] = pivot_score;

        long blend_score = Helper.quick_edit_dist(blend_result, prob_matrix);
        blendScoreTotal += blend_score;
        blendScores[j] = blend_score;

        long hybrid_score = Helper.quick_edit_dist(fix, prob_matrix); 
        hybridScoreTotal += hybrid_score;
        hybridScores[j] = hybrid_score;

        long ls_score = Helper.quick_edit_dist(ls_result, prob_matrix); 
        lsScoreTotal += ls_score;
        lsScores[j] = ls_score;

        }

        System.out.println("Finish");
        System.out.println();

        System.out.println("Alg1 0.5 times: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(pivotTimes[i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("Alg1 0.5 scores: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(pivotScores[i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("Alg1 0.3 times: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(hybridTimes[i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("Alg1 0.3 scores: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(hybridScores[i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("Alg1 0.1 times: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(blendTimes[i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("Alg1 0.1 scores: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(blendScores[i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("Alg2 times: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(lsTimes[i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("Alg2 scores: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(lsScores[i] + " ");
        System.out.println();
        System.out.println();

        System.out.println("Average Alg1 0.5 time: " + pivotTimeTotal / ((double) ROUNDS));
        System.out.println("Average Alg1 0.3 time: " + hybridTimeTotal / ((double) ROUNDS));
        System.out.println("Average Alg1 0.1 time: " + blendTimeTotal / ((double) ROUNDS));
        System.out.println("Average Alg2 time: " + lsTimeTotal / ((double) ROUNDS));
        System.out.println();
        System.out.println("Average Alg1 0.5 score: " + pivotScoreTotal / ((double) ROUNDS));
	    System.out.println("Average Alg1 0.3 score: " + hybridScoreTotal / ((double) ROUNDS));
	    System.out.println("Average Alg1 0.1 score: " + blendScoreTotal / ((double) ROUNDS));
	    System.out.println("Average Alg2 score: " + lsScoreTotal / ((double) ROUNDS));

        System.out.println();
        System.out.println("Average Alg1 0.5 num clusters: " + pivotNumClusters / ((double) ROUNDS));
        System.out.println("Average Alg1 0.3 num clusters: " + hybridNumClusters / ((double) ROUNDS));
        System.out.println("Average Alg1 0.1 num clusters: " + blendNumClusters / ((double) ROUNDS));
        System.out.println("Average Alg2 num clusters: " + lsNumClusters / ((double) ROUNDS));

        System.out.println();
        System.out.println("Average Alg1 0.5 max cluster size: " + largestPivotCluster / ((double) ROUNDS));
        System.out.println("Average Alg1 0.3 max cluster size: " + largestHybridCluster / ((double) ROUNDS));
        System.out.println("Average Alg1 0.1 max cluster size: " + largestBlendCluster / ((double) ROUNDS));
        System.out.println("Average Alg2 max cluster size: " + largestLsCluster / ((double) ROUNDS));

        System.out.println();


        }



    }


    private static class BoundComparator implements Comparator<Integer[]> {  
        public int compare (Integer[] a1, Integer[] a2) {  

            if (a1[0] == a2[0]) {   
                return 0;
            } else if (a1[0] > a2[0]) {  
                return 1;
            } else if (a1[0] < a2[0]) {
                return -1;
            }
            return 0;  
        }  
    }     

    public static ArrayList<ArrayList<Integer>> RoundingAlg1(MPVariable[][] var_array, ArrayList<Integer> size_bounds, double alpha) {

        int num_nodes = var_array.length;

        ArrayList<Integer[]> boundSort = new ArrayList<Integer[]>();
        for (int i = 0; i < num_nodes; i++) {
            Integer[] entry = {size_bounds.get(i), i};
            boundSort.add(entry);
        }
        boundSort.sort(new BoundComparator()); // sort by increasing bound sizes

        ArrayList<Integer> permutation = new ArrayList<Integer>();
        for (int i = 0; i < num_nodes; i++) 
            permutation.add(boundSort.get(i)[1]);

        HashSet<Integer> available_nodes = new HashSet<Integer>();
        available_nodes.addAll(permutation);
        Collections.shuffle(permutation); // create random order for choosing pivot
    
        ArrayList<ArrayList<Integer>> clusters = new ArrayList<ArrayList<Integer>>(); 
    
        int cur_index = 0;
    
        while (!available_nodes.isEmpty()) {
            // pick random node
            while (!available_nodes.contains(permutation.get(cur_index)))
                cur_index++;
    
            int cur_node = permutation.get(cur_index);
    
            ArrayList<Integer> cur_cluster = new ArrayList<Integer>();
            double sum = 0.0;
    
            for (int i = 0; i < num_nodes; i++) {
                if (available_nodes.contains(i) && i != cur_node) {
                    double prob = var_array[cur_node][i].solutionValue();
                    if (i > cur_node)
                        prob = var_array[i][cur_node].solutionValue(); // get x- value
                    if (prob <= alpha) {
                        sum += prob;
                        cur_cluster.add(i);
                    }
                }
            }

            if (sum >= alpha * (cur_cluster.size() + 1) / 2) {
                cur_cluster = new ArrayList<Integer>();
            }

            cur_cluster.add(cur_node);

    
            available_nodes.removeAll(cur_cluster);
            clusters.add(cur_cluster);            
        }
            
        return clusters;

    }

    public static <K, V extends Comparable<V> > TreeMap<K, V>
    valueSort(final Map<K, V> map) // https://www.geeksforgeeks.org/how-to-sort-a-treemap-by-value-in-java/
    {
        // Static Method with return type Map and
        // extending comparator class which compares values
        // associated with two keys
        Comparator<K> valueComparator = new Comparator<K>() {
            
                  // return comparison results of values of
                  // two keys
                  public int compare(K k1, K k2)
                  {
                      int comp = map.get(k1).compareTo(
                          map.get(k2));
                      if (comp == 0)
                          return 1;
                      else
                          return comp;
                  }
            
              };
        
        // SortedMap created using the comparator
        TreeMap<K, V> sorted = new TreeMap<K, V>(valueComparator);
        
        sorted.putAll(map);
        
        return sorted;
    }



    public static ArrayList<ArrayList<Integer>> RoundingAlg2(MPVariable[][] var_array, ArrayList<Integer> size_bounds) {

        int num_nodes = var_array.length;
        ArrayList<Integer[]> boundSort = new ArrayList<Integer[]>();
        for (int i = 0; i < num_nodes; i++) {
            Integer[] entry = {size_bounds.get(i), i};
            boundSort.add(entry);
        }
        boundSort.sort(new BoundComparator()); // sort by increasing bound sizes

        ArrayList<Integer> permutation = new ArrayList<Integer>();
        for (int i = 0; i < num_nodes; i++) 
            permutation.add(boundSort.get(i)[1]);
        HashSet<Integer> available_nodes = new HashSet<Integer>();
        available_nodes.addAll(permutation);
        Collections.shuffle(permutation); // create random order for choosing pivot
    
        ArrayList<ArrayList<Integer>> clusters = new ArrayList<ArrayList<Integer>>(); 
    
        int cur_index = 0;
        // double alpha = 0.37228; // (-5 + sqrt(33)) / 2

    
        while (!available_nodes.isEmpty()) {
            // pick random node
            while (!available_nodes.contains(permutation.get(cur_index)))
                cur_index++;
    
            int cur_node = permutation.get(cur_index);
            int cur_size_limit = size_bounds.get(cur_node);
    
            TreeMap<Integer, Double> cur_values = new TreeMap<Integer, Double>();
    
            for (int i = 0; i < num_nodes; i++) {
                if (available_nodes.contains(i) && i != cur_node) {
                    double prob = var_array[cur_node][i].solutionValue();
                    if (i > cur_node)
                        prob = var_array[i][cur_node].solutionValue(); // get x- value
                    if (prob <= 0.5)
                        cur_values.put(i, prob);
                }
            }

            if (cur_values.size() >= cur_size_limit) { // pick K-1 smallest values (not including cur_node)
                cur_values = valueSort(cur_values);
            }

            double sum = 0.0;
            ArrayList<Integer> cur_cluster = new ArrayList<Integer>();

            for (Map.Entry<Integer, Double> item : cur_values.entrySet()) {
                cur_cluster.add(item.getKey());
                sum += item.getValue();
                if (cur_cluster.size() == cur_size_limit - 1) // assume K >= 2 ...
                    break;
            }

            if (sum > (cur_cluster.size() + 1) / (2 * cur_size_limit)) {
                cur_cluster = new ArrayList<Integer>();
            }

            cur_cluster.add(cur_node);
    
            available_nodes.removeAll(cur_cluster);
            clusters.add(cur_cluster);            
        }
            
        return clusters;

    }

}
