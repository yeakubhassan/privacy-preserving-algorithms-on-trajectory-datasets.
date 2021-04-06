import java.util.ArrayList;

public class ptt {
    ArrayList<String> projection = new ArrayList<String>();
    ArrayList<Integer> trajectories = new ArrayList<Integer>();
    ArrayList<ptt> children = new ArrayList<ptt>();

    public ArrayList<String> getProjection() {
        return projection;
    }

    public void setProjection(ArrayList<String> projection) {
        this.projection = projection;
    }

    public ArrayList<Integer> getTrajectories() {
        return trajectories;
    }

    public void setTrajectories(ArrayList<Integer> trajectories) {
        this.trajectories = trajectories;
    }

    public ArrayList<ptt> getChildren() {
        return children;
    }

    public void setChildren(ArrayList<ptt> children) {
        this.children = children;
    }

    public void add_trajectory(int i)
    {
        this.trajectories.add(i);
    }
    public void add_children(ptt c)
    {
        children.add(c);
    }

    public void remove_trajectories(int t)
    {
        int index = this.trajectories.indexOf(t);
        this.trajectories.remove(index);
    }

    public void remove_child(int index)
    {
        this.children.remove(index);
    }
}

