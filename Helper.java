import java.io.*;
import java.util.*;

public class Helper {

    ArrayList<ArrayList<Double[]>> matches;
    ArrayList<ArrayList<Integer>> index_mapping;
    ArrayList<HashSet<Integer>> my_labels;
    long num_matches;

    public Helper(ArrayList<ArrayList<Double[]>> matches, ArrayList<ArrayList<Integer>> index_mapping) {
        /*
         * Workaround for multple returns from Consensus Clustering pre-processing
         */
        this.matches = matches;
        this.index_mapping = index_mapping;
    }

    public Helper(ArrayList<HashSet<Integer>> my_labels, long num_matches) {
        // WORKAROUND FOR MULTIPLE RETURNS FROM COMMUNITY DISAGREEMENT
        this.my_labels = my_labels;
        this.num_matches = num_matches;
    }

    public ArrayList<ArrayList<Double[]>> getMatches() {
        return matches;
    }

    public ArrayList<ArrayList<Integer>> getMapping () {
        return index_mapping;
    }

    public ArrayList<HashSet<Integer>> getLabels () {
        return my_labels;
    }

    public long getNumMatches () {
        return num_matches;
    }

    public static ArrayList<ArrayList<Integer>> read_large_network_relabel(String file_name, String delimiter) {
        /* Read in 0/1 graph from social network, input as adjacency list 
         * Relabel nodes as we encounter them, in case of non-contiguous numbering in data file
         */

        HashMap<Integer, Integer> my_mapping = new HashMap<Integer, Integer>();

        try {
            BufferedReader myReader = new BufferedReader(new FileReader(file_name));
            
            String data = myReader.readLine();
            String[] split = data.split(" ");
            int num_pts = Integer.parseInt(split[0]);

            ArrayList<ArrayList<Integer>> my_graph = new ArrayList<ArrayList<Integer>>();
            for (int i = 0; i < num_pts; i++)
                my_graph.add(new ArrayList<Integer>());

            int cur_index = 0;
            while ((data = myReader.readLine()) != null) {

                split = data.split(delimiter);
                int x = Integer.parseInt(split[0]);
                int y = Integer.parseInt(split[1]);

                // RELABEL
                if (!my_mapping.containsKey(x)) {
                    my_mapping.put(x, cur_index);
                    cur_index += 1;
                }    
                if (!my_mapping.containsKey(y)) {
                    my_mapping.put(y, cur_index);
                    cur_index += 1; 
                }

                // ADD TO GRAPH
                int map_x = my_mapping.get(x);
                int map_y = my_mapping.get(y);
                if (map_x != map_y) { // avoid self-edges
                    my_graph.get(map_x).add(map_y);
                    my_graph.get(map_y).add(map_x);
                }
            }
            myReader.close();
            return my_graph;
          } catch (Exception e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
            return null;
          }
    
    }

    public static ArrayList<Integer> read_size_bounds(String file_name) {

        ArrayList<Integer> size_bounds = new ArrayList<Integer>();

        try {
            BufferedReader myReader = new BufferedReader(new FileReader(file_name));
            
            String data;
            while ((data = myReader.readLine()) != null) {

                int x = Integer.parseInt(data);
                size_bounds.add(x);
            }
            myReader.close();
            return size_bounds;

          } catch (Exception e) {
            System.out.println("An error occurred. Flip!");
            e.printStackTrace();
            return null;
          }
    
    }

    public static HashMap<Integer, Integer> large_network_get_node_labels(String file_name, String delimiter) {
        /* Read in 0/1 graph from social network, input as adjacency list 
         * Relabel nodes as we encounter them, in case of non-contiguous numbering in data file
         */

        HashMap<Integer, Integer> my_mapping = new HashMap<Integer, Integer>();

        try {
            BufferedReader myReader = new BufferedReader(new FileReader(file_name));
            
            String data = myReader.readLine();
            String[] split = data.split(" ");
            int num_pts = Integer.parseInt(split[0]);

            int cur_index = 0;
            while ((data = myReader.readLine()) != null) {

                split = data.split(delimiter);
                int x = Integer.parseInt(split[0]);
                int y = Integer.parseInt(split[1]);

                // RELABEL
                if (!my_mapping.containsKey(x)) {
                    my_mapping.put(x, cur_index);
                    cur_index += 1;
                }    
                if (!my_mapping.containsKey(y)) {
                    my_mapping.put(y, cur_index);
                    cur_index += 1; 
                }

            }
            myReader.close();
            return my_mapping;
          } catch (Exception e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
            return null;
          }
    
    }

    public static Helper get_community_labels(HashMap<Integer, Integer> my_mapping, String communities, String delimiter) {

        // set of original node labels
        ArrayList<HashSet<Integer>> community_labels = new ArrayList<HashSet<Integer>>();
        int num_nodes = my_mapping.keySet().size();
        for (int i = 0; i < num_nodes; i++)
            community_labels.add(new HashSet<Integer>());
        long num_edges = 0;

        try {            
            System.out.println("start reading community edges");
            BufferedReader myReader2 = new BufferedReader(new FileReader(communities));
            String[] split = {};
            String data = "";
            int community_num = 0;
            while ((data = myReader2.readLine()) != null) {
                split = data.split(delimiter);
                for (int i = 0; i < split.length; i++) {
                    int first_val = my_mapping.get(Integer.parseInt(split[i]));
                    for (int j = i+1; j < split.length; j++) {
                        int second_val = my_mapping.get(Integer.parseInt(split[j]));
                        if (Collections.disjoint(community_labels.get(first_val), community_labels.get(second_val)))
                            num_edges++;
                    }
                    community_labels.get(first_val).add(community_num);
                }
                community_num++;
            }

            myReader2.close();
            System.out.println("finished reading community edges");

            return new Helper(community_labels, num_edges);            


        } catch (Exception e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
            return null;
        }
    }


    // TAKES TOO LONG
    public static HashSet<String> get_community_edges(HashMap<Integer, Integer> my_mapping, String communities, String delimiter) {

        // set of original node labels
        HashSet<String> community_labels = new HashSet<String>();


        try {            

            BufferedReader myReader2 = new BufferedReader(new FileReader(communities));
            String[] split = {};
            String data = "";
            while ((data = myReader2.readLine()) != null) {
                split = data.split(delimiter);
                ArrayList<String> repeat_labels = new ArrayList<String>();
                for (int i = 0; i < split.length; i++) {
                    for (int j = i+1; j < split.length; j++) {

                        int first_val = my_mapping.get(Integer.parseInt(split[i]));
                        int second_val = my_mapping.get(Integer.parseInt(split[j]));
                        if (first_val < second_val)
                            repeat_labels.add(first_val + "," + second_val);
                        else
                            repeat_labels.add(second_val + "," + first_val);

                    }

                }
                community_labels.addAll(repeat_labels);

            }

            myReader2.close();
            System.out.println("finished reading community edges");


            return community_labels;            


        } catch (Exception e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
            return null;
        }
    }

    public static long communityDisagreement(ArrayList<ArrayList<Integer>> clustering, ArrayList<HashSet<Integer>> community_labels, long shared) {

        long score = 0;
        long num_used = 0;

        for (int k = 0; k < clustering.size(); k++) {
            ArrayList<Integer> cur_cluster = clustering.get(k);
            for (int i = 0; i < cur_cluster.size(); i++) {
                int first_val = cur_cluster.get(i);
                for (int j = i+1; j < cur_cluster.size(); j++) {                  
                    int second_val = cur_cluster.get(j);

                    if (Collections.disjoint(community_labels.get(first_val), community_labels.get(second_val)))
                        score++;
                    else
                        num_used++;
                }

            }
        }
        long remainder = shared - num_used;

        return score + remainder;
    
    }


    // NOT PRACTICAL 
    public static long communityDisagreementString(ArrayList<ArrayList<Integer>> clustering, HashSet<String> community_labels) {

        long shared = community_labels.size();
        long score = 0;
        long num_used = 0;

        for (int k = 0; k < clustering.size(); k++) {
            ArrayList<Integer> cur_cluster = clustering.get(k);
            for (int i = 0; i < cur_cluster.size(); i++) {
                int first_val = cur_cluster.get(i);
                for (int j = i+1; j < cur_cluster.size(); j++) {                  
                    int second_val = cur_cluster.get(j);
                    String comparison = "";
                    if (first_val < second_val)
                        comparison = first_val +"," + second_val;
                    else
                        comparison = second_val + "," + first_val;

                    if (!community_labels.contains(comparison))
                        score++;
                    else
                        num_used++;
                }

            }
        }
        long remainder = shared - num_used;

        return score + remainder;
    
    }



    public static ArrayList<ArrayList<ArrayList<Integer>>> cluster_categorical_data(String file_name, String delimiter) {
        /* Read in 0/1 graph from social network, input as adjacency list 
         * Relabel nodes as we encounter them, in case of non-contiguous numbering in data file
         */

        try {
            BufferedReader myReader = new BufferedReader(new FileReader(file_name));
            
            String data = myReader.readLine();
            String[] split = data.split(delimiter);
            int num_values = split.length;

            ArrayList<ArrayList<ArrayList<Integer>>> my_clusterings = new ArrayList<ArrayList<ArrayList<Integer>>>();
            ArrayList<HashMap<String, Integer>> my_mapping = new ArrayList<HashMap<String, Integer>>();

            // Initialize with first line of data
            for (int i = 0; i < num_values; i++) {
                ArrayList<ArrayList<Integer>> cur_clustering = new ArrayList<ArrayList<Integer>>();
                ArrayList<Integer> cur_cluster = new ArrayList<Integer>();
                cur_cluster.add(0); // node 0
                cur_clustering.add(cur_cluster);
                my_clusterings.add(cur_clustering);
                
                // add mapping to current value of attribute
                HashMap<String, Integer> cur_map = new HashMap<String, Integer>();
                cur_map.put(split[i], 0); // split[i] maps to cluster 0
                my_mapping.add(cur_map);
            }

            int cur_index = 1;
            while ((data = myReader.readLine()) != null) {

                split = data.split(delimiter);
                /*
                if (split.length != num_values) {
                    System.out.println("Index: " + cur_index);
                    System.out.println(split.toString());
                }
                */


                for (int i = 0; i < num_values; i++) {
                    ArrayList<ArrayList<Integer>> cur_clustering = my_clusterings.get(i);
                    HashMap<String, Integer> cur_map = my_mapping.get(i);
                    if (cur_map.containsKey(split[i])) {
                        int cluster_index = cur_map.get(split[i]);
                        cur_clustering.get(cluster_index).add(cur_index);
                    } else {
                        int num_clusters = cur_clustering.size();
                        cur_map.put(split[i], num_clusters);
                        ArrayList<Integer> cur_cluster = new ArrayList<Integer>();
                        cur_cluster.add(cur_index);
                        cur_clustering.add(cur_cluster);

                    }
                    
                }

                if (cur_index % 1000 == 0)
                    System.out.println("Finished reading " + cur_index + " lines");

                cur_index += 1;
            }
            myReader.close();

            // FILTER CLUSTERINGS
            /*
            ArrayList<ArrayList<ArrayList<Integer>>> filtered_clusterings = new ArrayList<ArrayList<ArrayList<Integer>>>();
            for (int i = 0; i < my_clusterings.size(); i++) {
                ArrayList<ArrayList<Integer>> candidate_clustering = my_clusterings.get(i);
                if (candidate_clustering.size() <= cur_index / 2)
                    filtered_clusterings.add(candidate_clustering);
            }
            */

            return my_clusterings;
          } catch (Exception e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
            return null;
          }
    
    }

    public static HashMap<Integer, Integer> get_clustering_map (ArrayList<ArrayList<Integer>> clustering) {

        HashMap<Integer, Integer> cluster_map = new HashMap<Integer, Integer>();
        for (int i = 0; i < clustering.size(); i++) {
            ArrayList<Integer> cur_cluster = clustering.get(i);
            for (int j = 0; j < cur_cluster.size(); j++) {
                cluster_map.put(cur_cluster.get(j), i);               
            }
        }

        return cluster_map;

    }

    public static long quick_edit_dist(ArrayList<ArrayList<Integer>> clustering, ArrayList<ArrayList<Integer>> prob_matrix) {
        /*
           Quick edit distance for large 0/1 graphs
        */

        long dist = 0;
        HashMap<Integer, Integer> cluster_map = get_clustering_map(clustering);
        int[] cluster_sizes = new int[clustering.size()];
        for (int i = 0; i < cluster_sizes.length; i++)
            cluster_sizes[i] = clustering.get(i).size();

        // loop over nodes in prob_graph
        for (int i = 0; i < prob_matrix.size(); i++) {
            int cur_label = cluster_map.get(i);
            ArrayList<Integer> neighbors = prob_matrix.get(i);

            long x1 = 0; // neighbors of i with same cluster label
            long x2 = 0; // neighbors of i with different cluster label

            for (int j = 0; j < neighbors.size(); j++) {
                int cur_node = neighbors.get(j);
                if (i < cur_node) {
                    if (cluster_map.get(cur_node) == cur_label) {
                        x1++;
                    } else {
                        x2++;
                    }
                }
            }
            long update = x2 + cluster_sizes[cur_label] - x1 - 1;
            // if (cluster_sizes[cur_label] - x1 - 1 < 0) {
            //     System.out.println(cluster_sizes[cur_label] - x1 - 1);
            // }
            dist += update;
            cluster_sizes[cur_label]--; 

        }

        /*
        for (int i = 0; i < cluster_sizes.length; i++) {
            if(cluster_sizes[i] != 0)
                System.out.print(i + ": " + cluster_sizes[i]);
        }
        */


        return dist;
    }

    public static double[] get_precision_recall(ArrayList<ArrayList<Integer>> clustering, ArrayList<ArrayList<Integer>> prob_matrix) {
        /*
           Quick edit distance for large 0/1 graphs
        */

        int dist = 0;
        int total_pos = 0;
        HashMap<Integer, Integer> cluster_map = get_clustering_map(clustering);

        // loop over nodes in prob_graph
        for (int i = 0; i < prob_matrix.size(); i++) {
            int cur_label = cluster_map.get(i);
            ArrayList<Integer> neighbors = prob_matrix.get(i);

            for (int j = 0; j < neighbors.size(); j++) {
                int cur_node = neighbors.get(j);
                if (i < cur_node) {
                    if (cluster_map.get(cur_node) == cur_label) {
                        dist++;
                    } 
                    total_pos++;
                }
            }
           

        }
        int clustered_pairs = 0;
        for (int i = 0; i < clustering.size(); i++) {
            int cur_cluster_size = clustering.get(i).size();
            int cur_pairs = cur_cluster_size * (cur_cluster_size - 1) / 2;
            clustered_pairs += cur_pairs;
        }

        double precision = ((double) dist) / clustered_pairs;
        double recall = ((double) dist) / total_pos;

        return new double[] {precision, recall};


    }


    // --- OLD CODE --- 

    public static ArrayList<ArrayList<Integer>> read_large_network_relabel2 (String file_name, String delimiter) {
        // 0/1 graph from social network, input as adjacency list
        // replacing prob_matrix with a dictionary...
    
        // read dataset from file
        HashMap<Integer, Integer> my_mapping = new HashMap<Integer, Integer>();

        try {
            // File myObj = new File(file_name);
            // Scanner myReader = new Scanner(myObj);
            BufferedReader myReader = new BufferedReader(new FileReader(file_name));
            
            String data = myReader.readLine();
            String[] split = data.split(" ");
            // int num_pts = Integer.parseInt(split[0]);

            int cur_index = 0;
            while ((data = myReader.readLine()) != null) {
                // data = myReader.nextLine();
                split = data.split(delimiter);
                int x = Integer.parseInt(split[0]);
                int y = Integer.parseInt(split[1]);
                // System.out.println(x + " and " + y);
                if (!my_mapping.containsKey(x)) {
                    my_mapping.put(x, cur_index);
                    cur_index += 1;
                }    
                if (!my_mapping.containsKey(y)) {
                    my_mapping.put(y, cur_index);
                    cur_index += 1; 
                }             
            }
            myReader.close();
            // return my_graph;
          } catch (Exception e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
            // return null;
          }

        try {
            // File myObj = new File(file_name);
            // Scanner myReader = new Scanner(myObj);
            BufferedReader myReader = new BufferedReader(new FileReader(file_name));
            
            String data = myReader.readLine();
            String[] split = data.split(" ");
            int num_pts = Integer.parseInt(split[0]);

            ArrayList<ArrayList<Integer>> my_graph = new ArrayList<ArrayList<Integer>>();
            for (int i = 0; i < num_pts; i++)
                my_graph.add(new ArrayList<Integer>());

            while ((data = myReader.readLine()) != null) {
                // data = myReader.nextLine();
                split = data.split(delimiter);
                int x = Integer.parseInt(split[0]);
                int y = Integer.parseInt(split[1]);

                int map_x = my_mapping.get(x);
                int map_y = my_mapping.get(y);

                my_graph.get(map_x).add(map_y);
                my_graph.get(map_y).add(map_x);

            }
            myReader.close();
            return my_graph;
          } catch (Exception e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
            return null;
          }
    
    }


    public static ArrayList<ArrayList<Integer>> read_large_network_graph3(String file_name, String delimiter) {
        // 0/1 graph from social network, input as adjacency list
        // replacing prob_matrix with a dictionary...
    
        // read dataset from file
        try {
            File myObj = new File(file_name);
            Scanner myReader = new Scanner(myObj);
            
            String data = myReader.nextLine();
            String[] split = data.split(" ");
            int num_pts = Integer.parseInt(split[0]);

            ArrayList<ArrayList<Integer>> my_graph = new ArrayList<ArrayList<Integer>>();
            for (int i = 0; i < num_pts; i++)
                my_graph.add(new ArrayList<Integer>());

            while (myReader.hasNextLine()) {
                data = myReader.nextLine();
                split = data.split(delimiter);
                int x = Integer.parseInt(split[0]);
                int y = Integer.parseInt(split[1]);
                my_graph.get(x).add(y);
                my_graph.get(y).add(x);


            }
            myReader.close();
            return my_graph;
          } catch (FileNotFoundException e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
            return null;
          }
    
    }

    public static ArrayList<HashSet<Integer>> read_large_network_graph_hash (String file_name, String delimiter) {
        // 0/1 graph from social network, input as adjacency list
        // replacing prob_matrix with a dictionary...
    
        // read dataset from file
        try {
            File myObj = new File(file_name);
            Scanner myReader = new Scanner(myObj);
            
            String data = myReader.nextLine();
            String[] split = data.split(" ");
            int num_pts = Integer.parseInt(split[0]);

            ArrayList<HashSet<Integer>> my_graph = new ArrayList<HashSet<Integer>>();
            for (int i = 0; i < num_pts; i++)
                my_graph.add(new HashSet<Integer>());

            while (myReader.hasNextLine()) {
                data = myReader.nextLine();
                split = data.split(delimiter);
                int x = Integer.parseInt(split[0]);
                int y = Integer.parseInt(split[1]);
                my_graph.get(x).add(y);
                my_graph.get(y).add(x);


            }
            myReader.close();
            return my_graph;
          } catch (FileNotFoundException e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
            return null;
          }
    
    }

    public static ArrayList<HashSet<Integer>> read_large_network_relabel (String file_name) {
        // 0/1 graph from social network, input as adjacency list
        // replacing prob_matrix with a dictionary...
    
        // read dataset from file
        HashMap<Integer, Integer> my_mapping = new HashMap<Integer, Integer>();

        try {
            // File myObj = new File(file_name);
            // Scanner myReader = new Scanner(myObj);
            BufferedReader myReader = new BufferedReader(new FileReader(file_name));
            
            String data = myReader.readLine();
            // String[] split = data.split(" ");
            // int num_pts = Integer.parseInt(split[0]);
            int cur_index = 0;
            while ((data = myReader.readLine()) != null) {
                // data = myReader.nextLine();
                String [] split = data.split(",");
                int x = Integer.parseInt(split[0]);
                int y = Integer.parseInt(split[1]);
                if (!my_mapping.containsKey(x)) {
                    my_mapping.put(x, cur_index);
                    cur_index += 1;
                }    
                if (!my_mapping.containsKey(y)) {
                    my_mapping.put(y, cur_index);
                    cur_index += 1; 
                }             
            }
            myReader.close();
            // return my_graph;
          } catch (Exception e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
            // return null;
          }

        try {
            // File myObj = new File(file_name);
            // Scanner myReader = new Scanner(myObj);
            BufferedReader myReader = new BufferedReader(new FileReader(file_name));
            
            String data = myReader.readLine();
            String[] split = data.split(" ");
            int num_pts = Integer.parseInt(split[0]);

            ArrayList<HashSet<Integer>> my_graph = new ArrayList<HashSet<Integer>>();
            for (int i = 0; i < num_pts; i++)
                my_graph.add(new HashSet<Integer>());

            while ((data = myReader.readLine()) != null) {
                // data = myReader.nextLine();
                split = data.split(","); // \\s
                int x = Integer.parseInt(split[0]);
                int y = Integer.parseInt(split[1]);

                int map_x = my_mapping.get(x);
                int map_y = my_mapping.get(y);

                my_graph.get(map_x).add(map_y);
                my_graph.get(map_y).add(map_x);

            }
            myReader.close();
            return my_graph;
          } catch (Exception e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
            return null;
          }
    
    }




    public static boolean read_very_large_network_relabel (String file_name, String path, int step) {
        // 0/1 graph from social network, input as adjacency list
        // replacing prob_matrix with a dictionary...
    
        // read dataset from file
        HashMap<Integer, Integer> my_mapping = new HashMap<Integer, Integer>();

        int nodesAtATime = 10000;

        try {
            // File myObj = new File(file_name);
            // Scanner myReader = new Scanner(myObj);
            
            BufferedReader myReader = new BufferedReader(new FileReader(file_name));

            String data = myReader.readLine();
            // String[] split = data.split(" ");
            // int num_pts = Integer.parseInt(split[0]);
            int cur_index = 0;
            while ((data = myReader.readLine()) != null) {
                // data = myReader.nextLine();
                String [] split = data.split(","); // \\s
                int x = Integer.parseInt(split[0]);
                int y = Integer.parseInt(split[1]);
                if (!my_mapping.containsKey(x)) {
                    my_mapping.put(x, cur_index);
                    cur_index += 1;
                }    
                if (!my_mapping.containsKey(y)) {
                    my_mapping.put(y, cur_index);
                    cur_index += 1; 
                }             
            }
            myReader.close();
            // return my_graph;
          } catch (Exception e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
            // return null;
          } 

        try {
            // File myObj = new File(file_name);

            // iterate through 10,000 nodes at a time
            int num_nodes = my_mapping.size();
            int lower_limit = 0;          
            int upper_limit = nodesAtATime;
            int round = 0;
            int span = upper_limit;
            if (span > num_nodes)
                span = num_nodes;

 
            while (lower_limit < num_nodes){
                BufferedReader myReader = new BufferedReader(new FileReader(file_name));

                String data = myReader.readLine();
                String[] split = data.split(" ");

                ArrayList<ArrayList<Integer>> my_graph = new ArrayList<ArrayList<Integer>>();
                for (int i = 0; i < span; i++)
                    my_graph.add(new ArrayList<Integer>());

                while ((data = myReader.readLine()) != null) {
                    // data = myReader.nextLine();
                    split = data.split("\\s");
                    int x = Integer.parseInt(split[0]);
                    int y = Integer.parseInt(split[1]);

                    int map_x = my_mapping.get(x);
                    int map_y = my_mapping.get(y);

                    if (lower_limit <= map_x && map_x < upper_limit)
                        my_graph.get(map_x % nodesAtATime).add(map_y);
                    if (lower_limit <= map_y && map_y < upper_limit)
                        my_graph.get(map_y % nodesAtATime).add(map_x);


                }
                myReader.close();

                // write results to file

                int iterations = nodesAtATime / step;
                int lowerStep = 0;
                int upperStep = step;
                if (upperStep > my_graph.size())
                    upperStep = my_graph.size();
                for (int k = 0; k < iterations; k++) {

                    String x_file = path + "/" + round + ".txt";
                    File myFile = new File(x_file); // create file
                    BufferedWriter bw = new BufferedWriter(new FileWriter(x_file));
                    for (int i = lowerStep; i < upperStep; i++) {
                        ArrayList<Integer> neighbors = my_graph.get(i);
                        if (!neighbors.isEmpty()) {
                            bw.write("" + neighbors.get(0));
                        }
                        for (int j = 1; j < neighbors.size(); j++) {
                            bw.write("," + neighbors.get(j));
                        }
                        bw.newLine();
                    }
                    bw.close();

                    round += 1;
                    lowerStep += step;
                    upperStep += step;
                    if (upperStep > my_graph.size())
                        upperStep = my_graph.size();
                }

                lower_limit += span;          
                upper_limit += span;
                if (span > num_nodes)
                    span = num_nodes;

            }
            
            return true;
          } catch (Exception e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
            return false;
          }
    
    }


    public static int large_graph_prob_edit_dist(ArrayList<ArrayList<Integer>> clustering, ArrayList<HashSet<Integer>> prob_matrix) {
    /*
        Compute the edit distance between a given clustering and the 
        original probablisitic graph (in expectation)
        
        input graph is adjacency list
    */

    int dist = 0;

    for (int i = 0; i < clustering.size(); i++) {
        ArrayList<Integer> cur_cluster = clustering.get(i);
        for (int j = 0; j < cur_cluster.size(); j++) {

            int cur_node = cur_cluster.get(j);
            for (int k = j + 1; k < cur_cluster.size(); k++) {
                int other_node = cur_cluster.get(k);
                // add 1 if nodes aren't actually neighbors
                if (!prob_matrix.get(cur_node).contains(other_node))
                    dist += 1; 
            }
            for (int k = i + 1; k < clustering.size(); k++) {
                ArrayList<Integer> other_cluster = clustering.get(k);
                for (int l = 0; l < other_cluster.size(); l++) {
                    int other_node = other_cluster.get(l);
                    // add 1 if nodes are actually neighbors
                    if (prob_matrix.get(cur_node).contains(other_node))
                        dist += 1;
                }

            }

        }  
        // System.out.println("Finished cluster: " + i);              
    }
    
    return dist;
    }





    public static int large_graph_prob_edit_dist_hash (ArrayList<ArrayList<Integer>> clustering, ArrayList<HashSet<Integer>> prob_matrix) {
        /*
            Compute the edit distance between a given clustering and the 
            original probablisitic graph (in expectation)
            
            input graph is adjacency list
        */
    
        int dist = 0;
        HashMap<Integer, Integer> cluster_map = get_clustering_map(clustering);
        int neighbor = 0;
    
        for (int i = 0; i < prob_matrix.size(); i++) {            
            for (int j = i + 1; j < prob_matrix.size(); j++) {
                if (prob_matrix.get(i).contains(j)) {
                    neighbor = 1;
                } else {
                    neighbor = 0;
                }
                if (cluster_map.get(i).equals(cluster_map.get(j))) { // .equals() important for comparing Integer objects!
                    dist += (1 - neighbor);
                } else {
                    dist += neighbor;
                }
                
            }  
            // if ((i + 1) % 100 == 0)
            //     System.out.println("finished " + (i + 1) + " nodes");
        }
        
        return dist;
    }


    
}
