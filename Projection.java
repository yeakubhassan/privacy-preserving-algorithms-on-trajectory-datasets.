import java.io.Serializable;
import java.util.ArrayList;

public class Projection{
    ArrayList<String> proj = new ArrayList<String>();
    int ST;
    ArrayList<Projection> children = new ArrayList<Projection>();
    ArrayList<Integer> trajectories = new ArrayList<Integer>();
    boolean isProblem = false;
    public ArrayList<String> getProj() {
        return proj;
    }

    public ArrayList<Integer> getTrajectories() {
        return trajectories;
    }

    public void setTrajectories(ArrayList<Integer> trajectories) {
        this.trajectories = trajectories;
    }

    public void setProj(ArrayList<String> proj) {
        this.proj = proj;
    }

    public int getST() {
        return ST;
    }

    public void setST(int ST) {
        this.ST = ST;
    }

    public ArrayList<Projection> getChildren() {
        return children;
    }

    public void setChildren(ArrayList<Projection> children) {
        this.children = children;
    }

    public boolean isProblem() {
        return isProblem;
    }

    public void setProblem(boolean problem) {
        isProblem = problem;
    }
    public void add_trajectory(int i)
    {
        this.trajectories.add(i);
    }
    public void remove_trajectories(int t)
    {
        int index = this.trajectories.indexOf(t);
        this.trajectories.remove(index);
    }
    public void add_children(Projection c)
    {
        children.add(c);
    }
    public void remove_child(int index)
    {
        this.children.remove(index);
    }

    public void add_projection_data(String st)
    {
        this.proj.add(st);
    }
}
