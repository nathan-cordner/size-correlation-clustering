import java.util.*;
import java.io.*;

public class WriteNonUniformSizes {

    public static int getRandomNumber(int min, int max) {
        // https://www.baeldung.com/java-generating-random-numbers-in-range
        // too lazy to write it myself D:
        return (int) ((Math.random() * (max - min)) + min);
    }

    public static void main(String args[]){

        String data_set = args[0];
        String delimiter = "\\s"; 

        ArrayList<ArrayList<Integer>> prob_matrix = Helper.read_large_network_relabel("Data/"+ data_set + "/graph.txt", delimiter);

        int ROUNDS = 3;

        for (int j = 0; j < ROUNDS; j++) {

            ArrayList<Integer> size_bounds = new ArrayList<Integer>();
            for (int i = 0; i < prob_matrix.size(); i++){
                // random uniform between 1 and |N(i)| + 1 inclusive
                int cur_bound = getRandomNumber(1, prob_matrix.get(i).size() + 2); 
                // System.out.println(i + ": " + cur_bound + " out of " + (prob_matrix.get(i).size() + 1));
                size_bounds.add(cur_bound);
            }

            // output to file 
            try {
                String output_file = data_set + "_bounds" + j + ".txt";

                BufferedWriter writer = new BufferedWriter(new FileWriter(output_file));            

                for (int i = 0; i < size_bounds.size(); i++) {
                    writer.write(size_bounds.get(i).toString());
                    writer.newLine();
                }            

                writer.close();
            
            } catch (Exception e) {
                System.out.println("Ding nabit, an error occurred!");
                e.printStackTrace();
                return;
            }       
        }

    }
    
}
