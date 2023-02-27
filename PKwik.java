import java.util.*;

public class PKwik {

    // --- UNCONSTRAINED PIVOT ---

    // LARGE NETWORK: READ AS ADJACENCY LIST
    public static ArrayList<ArrayList<Integer>> pKwikClustering(ArrayList<ArrayList<Integer>> prob_matrix) {
        /*
            prob_matrix is now a dictionary of lists 
            index of graph is a node, list contains all neighbors of node
        */
        
        int num_nodes = prob_matrix.size();
        ArrayList<Integer> permutation = new ArrayList<Integer>();
        for (int i = 0; i < num_nodes; i++)
            permutation.add(i);
        HashSet<Integer> available_nodes = new HashSet<Integer>();
        available_nodes.addAll(permutation);
        Collections.shuffle(permutation); // create random order for choosing pivot

        ArrayList<ArrayList<Integer>> clusters = new ArrayList<ArrayList<Integer>>(); 
    
        // int loop1 = 0;
        // int loop2 = 0;
        // long start = 0;

        int cur_index = 0;

        while (!available_nodes.isEmpty()) {
            // pick random node
            while (!available_nodes.contains(permutation.get(cur_index)))
                cur_index++;

            int cur_node = permutation.get(cur_index);
            // available_nodes.remove(cur_index);

            ArrayList<Integer> cur_cluster = new ArrayList<Integer>();
            cur_cluster.add(cur_node);
      

            // start = System.currentTimeMillis();
            for (int neighbor : prob_matrix.get(cur_node)) {
                if (available_nodes.contains(neighbor)) {
                    cur_cluster.add(neighbor);
                }
            }
            // loop1 += System.currentTimeMillis() - start;

            // start = System.currentTimeMillis();
            available_nodes.removeAll(cur_cluster);
            // for (int i = 1; i < cur_cluster.size(); i++) 
            //     available_nodes.remove(cur_cluster.get(i));
            // loop2 += System.currentTimeMillis() - start;
            clusters.add(cur_cluster);            
        }
        // System.out.println("Loop 1 finished in: " + loop1 / 1000.0 + " s");
        // System.out.println("Loop 2 finished in: " + loop2 / 1000.0 + " s");
            
        return clusters;
    }

    // INPUT: probability matrix
    public static ArrayList<ArrayList<Integer>> pKwikClusteringProb(ArrayList<ArrayList<Double[]>> prob_matrix) {
        /*
            prob_matrix is now a dictionary of lists 
            index of graph is a node, list contains all neighbors of node
        */
        
        int num_nodes = prob_matrix.size();
        ArrayList<Integer> permutation = new ArrayList<Integer>();
        for (int i = 0; i < num_nodes; i++)
            permutation.add(i);
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
            cur_cluster.add(cur_node);

            for (Double[] values : prob_matrix.get(cur_node)) {
                double temp = values[0];
                int neighbor = (int) temp;
                double prob = values[1];

                if (prob >= 0.5 && available_nodes.contains(neighbor))
                    cur_cluster.add(neighbor);
            }

            available_nodes.removeAll(cur_cluster);
            clusters.add(cur_cluster);            
        }
            
        return clusters;
    }

    // --- MAX CLUSTER SIZE K ---

    // LARGE NETWORK: READ AS ADJACENCY LIST
    public static ArrayList<ArrayList<Integer>> maxKPivot(ArrayList<ArrayList<Integer>> prob_matrix, int max_size) {
        /*
            prob_matrix is now a dictionary of lists 
            index of graph is a node, list contains all neighbors of node
        */
        
        int num_nodes = prob_matrix.size();
        ArrayList<Integer> permutation = new ArrayList<Integer>();
        for (int i = 0; i < num_nodes; i++)
            permutation.add(i);
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
            cur_cluster.add(cur_node);

            ArrayList<Integer> neighbors = prob_matrix.get(cur_node);
            Collections.shuffle(neighbors); // add neighbors in random order up to max size

            for (int neighbor : neighbors) {
                if (available_nodes.contains(neighbor) && cur_cluster.size() < max_size) {
                    cur_cluster.add(neighbor);
                }
            }

            available_nodes.removeAll(cur_cluster);

            clusters.add(cur_cluster);            
        }
            
        return clusters;
    }


    // Non-Uniform Max Cluster Size

    public static ArrayList<ArrayList<Integer>> nonUniformMaxPivot(ArrayList<ArrayList<Integer>> prob_matrix, ArrayList<Integer> size_bounds) {
        /*
            prob_matrix is now a dictionary of lists 
            index of graph is a node, list contains all neighbors of node
        */
        
        int num_nodes = prob_matrix.size();
        ArrayList<Integer> permutation = new ArrayList<Integer>();
        for (int i = 0; i < num_nodes; i++)
            permutation.add(i);
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
            int max_size = size_bounds.get(cur_node);
            // System.out.println("cur_node " + cur_node + " has size limit " + max_size);

            ArrayList<Integer> cur_cluster = new ArrayList<Integer>();
            cur_cluster.add(cur_node);

            ArrayList<Integer> neighbors = prob_matrix.get(cur_node);
            Collections.shuffle(neighbors); // add neighbors in random order up to max size

            for (int neighbor : neighbors) {
                if (available_nodes.contains(neighbor) && cur_cluster.size() < max_size) {
                    // check size constraint of neighbor
                    int new_size = size_bounds.get(neighbor);
                    // System.out.println("Neighbor " + neighbor + " has size limit " + new_size);
                    if (new_size > cur_cluster.size()) {
                        cur_cluster.add(neighbor);
                        if (new_size < max_size)
                            max_size = new_size; // re-adjust cluster size limit 
                    }
                }
            }

            available_nodes.removeAll(cur_cluster);

            clusters.add(cur_cluster);            
        }
            
        return clusters;
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

    public static ArrayList<ArrayList<Integer>> nonUniformMaxPivotOrder(ArrayList<ArrayList<Integer>> prob_matrix, ArrayList<Integer> size_bounds) {
        /*
            prob_matrix is now a dictionary of lists 
            index of graph is a node, list contains all neighbors of node

            ORDER NODES BY INCREASING SIZE BOUNDS
        */
        
        int num_nodes = prob_matrix.size();

        ArrayList<Integer[]> boundSort = new ArrayList<Integer[]>();
        for (int i = 0; i < num_nodes; i++) {
            Integer[] entry = {size_bounds.get(i), i};
            boundSort.add(entry);
        }
        boundSort.sort(new BoundComparator()); // sort by increasing bound sizes

        ArrayList<Integer> permutation = new ArrayList<Integer>();
        for (int i = 0; i < num_nodes; i++) {
            permutation.add(boundSort.get(i)[1]); // use this order for Pivot
            // System.out.println(boundSort.get(i)[0]);
        }
        HashSet<Integer> available_nodes = new HashSet<Integer>();
        available_nodes.addAll(permutation);
        // Collections.shuffle(permutation); // create random order for choosing pivot

        ArrayList<ArrayList<Integer>> clusters = new ArrayList<ArrayList<Integer>>(); 
    
        int cur_index = 0;

        while (!available_nodes.isEmpty()) {
            // pick random node
            while (!available_nodes.contains(permutation.get(cur_index)))
                cur_index++;

            int cur_node = permutation.get(cur_index);
            int max_size = size_bounds.get(cur_node);

            ArrayList<Integer> cur_cluster = new ArrayList<Integer>();
            cur_cluster.add(cur_node);

            ArrayList<Integer> neighbors = prob_matrix.get(cur_node);
            Collections.shuffle(neighbors); // add neighbors in random order up to max size

            for (int neighbor : neighbors) {
                if (available_nodes.contains(neighbor) && cur_cluster.size() < max_size) {
                    // check size constraint of neighbor
                    int new_size = size_bounds.get(neighbor);
                    if (new_size > cur_cluster.size()) {
                        cur_cluster.add(neighbor);
                        if (new_size < max_size)
                            max_size = new_size; // re-adjust cluster size limit 
                    }
                }
            }

            available_nodes.removeAll(cur_cluster);

            clusters.add(cur_cluster);            
        }
            
        return clusters;
    }
    
}
