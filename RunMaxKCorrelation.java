import java.util.*;

public class RunMaxKCorrelation {

    public static void main(String args[]){

        String data_set = args[0];
        String delimiter = "\\s"; 

        ArrayList<ArrayList<Integer>> prob_matrix = Helper.read_large_network_relabel("Data/"+ data_set + "/graph.txt", delimiter);

        int ROUNDS = 10;

        int[] k_vals = {5, 10}; //, 15, 20, 25, 30, 35, 40, 45, 50};

        for (int q = 0; q < k_vals.length; q++) {

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
        ArrayList<ArrayList<Integer>> pivot_result = PKwik.maxKPivot(prob_matrix, k);
        long pivotTime = System.currentTimeMillis() - pivotStart;
        pivotTimeTotal += (pivotTime / 1000.0);
        pivotTimes[j] = pivotTime;
        // System.out.println("finish pivot " + j);

        long blendStart = System.currentTimeMillis();
        ArrayList<ArrayList<Integer>> blend_result = DNode.max_k_random_node_network(prob_matrix, k);
        long blendTime = System.currentTimeMillis() - blendStart;
        blendTimeTotal += (blendTime / 1000.0); 
        blendTimes[j] = blendTime;
        // System.out.println("finish rnode " + j);

        ArrayList<ArrayList<Integer>> copy1 = new ArrayList<ArrayList<Integer>>();
        for (int i = 0; i < pivot_result.size(); i++) {
            ArrayList<Integer> cluster_copy = new ArrayList<Integer>();
            cluster_copy.addAll(pivot_result.get(i));
            copy1.add(cluster_copy);
        }

        long hybridStart = System.currentTimeMillis();
        ArrayList<ArrayList<Integer>> fix = DNode.local_search_network_size_limit(copy1, prob_matrix, k, -1); // Hybrid.large_graph_fix_clusters(pivot_result, prob_matrix);
        long hybridTime = System.currentTimeMillis() - hybridStart;
        hybridTimeTotal += (hybridTime / 1000.0);
        hybridTimes[j] = hybridTime;
        // System.out.println("finish hybrid " + j);

        ArrayList<ArrayList<Integer>> copy2 = new ArrayList<ArrayList<Integer>>();
        for (int i = 0; i < blend_result.size(); i++) {
            ArrayList<Integer> cluster_copy = new ArrayList<Integer>();
            cluster_copy.addAll(blend_result.get(i));
            copy2.add(cluster_copy);
        }

        long lsStart = System.currentTimeMillis();
        ArrayList<ArrayList<Integer>> ls_result = DNode.local_search_network_size_limit(copy2, prob_matrix, k, -1);
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

        System.out.println("Pivot times: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(pivotTimes[i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("Pivot scores: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(pivotScores[i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("PLS times: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(hybridTimes[i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("PLS scores: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(hybridScores[i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("Vote times: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(blendTimes[i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("Vote scores: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(blendScores[i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("VLS times: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(lsTimes[i] + " ");
        System.out.println();
        System.out.println();
        System.out.println("VLS scores: ");
        for (int i = 0; i < ROUNDS; i++)
            System.out.print(lsScores[i] + " ");
        System.out.println();
        System.out.println();

        System.out.println("Average Pivot time: " + pivotTimeTotal / ((double) ROUNDS));
        System.out.println("Average PLS time: " + hybridTimeTotal / ((double) ROUNDS));
        System.out.println("Average Vote time: " + blendTimeTotal / ((double) ROUNDS));
        System.out.println("Average VLS time: " + lsTimeTotal / ((double) ROUNDS));
        System.out.println();
        System.out.println("Average Pivot score: " + pivotScoreTotal / ((double) ROUNDS));
	System.out.println("Average PLS score: " + hybridScoreTotal / ((double) ROUNDS));
	System.out.println("Average Vote score: " + blendScoreTotal / ((double) ROUNDS));
	System.out.println("Average VLS score: " + lsScoreTotal / ((double) ROUNDS));

        System.out.println();
        System.out.println("Average Pivot num clusters: " + pivotNumClusters / ((double) ROUNDS));
        System.out.println("Average PLS num clusters: " + hybridNumClusters / ((double) ROUNDS));
        System.out.println("Average Vote num clusters: " + blendNumClusters / ((double) ROUNDS));
        System.out.println("Average VLS num clusters: " + lsNumClusters / ((double) ROUNDS));

        System.out.println();
        System.out.println("Average Pivot max cluster size: " + largestPivotCluster / ((double) ROUNDS));
        System.out.println("Average PLS max cluster size: " + largestHybridCluster / ((double) ROUNDS));
        System.out.println("Average Vote max cluster size: " + largestBlendCluster / ((double) ROUNDS));
        System.out.println("Average VLS max cluster size: " + largestLsCluster / ((double) ROUNDS));

        System.out.println();


        }


    }
    
}
