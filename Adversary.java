import java.io.Serializable;
import java.util.ArrayList;

public class Adversary {
        int id;
        ArrayList<Projection> projections = new ArrayList<Projection>();

    public int getId() {
        return id;
    }

    public void setId(int id) {
        this.id = id;
    }

    public ArrayList<Projection> getProjections() {
        return projections;
    }

    public void setProjections(ArrayList<Projection> projections) {
        this.projections = projections;
    }

    public void remove_projection(int index)
    {
        this.projections.remove(index);
    }
}
