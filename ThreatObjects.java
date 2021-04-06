import java.util.ArrayList;

public class ThreatObjects {
    ArrayList<ProjectionList> projectionLists = new ArrayList<ProjectionList>();
    ArrayList<ProjectionSupport>allProjection = new ArrayList<ProjectionSupport>();
    ArrayList<Integer> problematicPair = new ArrayList<Integer>();
    int number_of_problems;

    public ArrayList<ProjectionSupport> getAllProjection() {
        return allProjection;
    }

    public void setAllProjection(ArrayList<ProjectionSupport> allProjection) {
        this.allProjection = allProjection;
    }

    public ArrayList<Integer> getProblematicPair() {
        return problematicPair;
    }

    public void setProblematicPair(ArrayList<Integer> problematicPair) {
        this.problematicPair = problematicPair;
    }

    public int getNumber_of_problems() {
        return number_of_problems;
    }

    public void setNumber_of_problems(int number_of_problems) {
        this.number_of_problems = number_of_problems;
    }

    public ArrayList<ProjectionList> getProjectionLists() {
        return projectionLists;
    }

    public void setProjectionLists(ArrayList<ProjectionList> projectionLists) {
        this.projectionLists = projectionLists;
    }
}
