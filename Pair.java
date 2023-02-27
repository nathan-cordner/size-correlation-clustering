// --- HELPER CLASS ---

public class Pair implements Comparable<Pair> {

    int key;
    int value;

    public Pair(int key, int value) {
        this.key = key;
        this.value = value;
    }

    public int getKey(){
       return key;
    }

    public void updateKey(int newKey) {
        this.key = newKey;
    }

    public void incrementKey() {
        this.key++;
    }

    public void decrementKey() {
        this.key--;
    }

    public int getValue() {
        return value;
    }

    public boolean equals(Object other) {

        if (!(other instanceof Pair)) {
            return false;
        }

        Pair o = (Pair) other;

        if (this.value == o.getValue())
            return true;
        return false;
    }

    public int compareTo(Pair other) {
        if (this.key < other.getKey())
            return -1;
        else if (this.key == other.getKey())
            return 0;
        else
            return 1;

    }

}
