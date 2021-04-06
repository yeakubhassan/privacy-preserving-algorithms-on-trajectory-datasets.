import java.util.ArrayList;

public class ProjectionSupport {
    private ArrayList<String> projection;
    private int ST;
    private int adversary;

    public int getAdversary() {
        return adversary;
    }

    public void setAdversary(int adversary) {
        this.adversary = adversary;
    }



    public ArrayList<String> getProjection() {
        return projection;
    }

    public ProjectionSupport() {
        this.ST = 0;
    }

    public void setProjection(ArrayList<String> projection) {
        this.projection = projection;
    }

    public int getST() {
        return ST;
    }

    public void setST(int ST) {
        this.ST = ST;
    }

    public void inc_st()
    {
        this.ST= this.ST +1;
    }

    @Override
    public String toString() {
//        return "ProjectionSupport{" +
//                "projection=" + projection +
//                ", ST=" + ST +
//                '}';

        return "projection=" + projection +
                ", ST=" + ST + ", adversary=" + adversary ;
    }
}
