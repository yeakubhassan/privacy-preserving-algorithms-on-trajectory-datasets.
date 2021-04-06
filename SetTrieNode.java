import java.util.ArrayList;

public class SetTrieNode {
    String node;
    ArrayList<SetTrieNode> child = new ArrayList<SetTrieNode>();
    int flag;

    public String getNode() {
        return node;
    }

    public void setNode(String node) {
        this.node = node;
    }

    public ArrayList<SetTrieNode> getChild() {
        return child;
    }

    public void setChild(ArrayList<SetTrieNode> child) {
        this.child = child;
    }

    public int getFlag() {
        return flag;
    }

    public void setFlag(int flag) {
        this.flag = flag;
    }
}
