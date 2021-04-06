import java.util.ArrayList;

public class ProjectionList {
    private ArrayList<String> projection;
    private String lembda;
    private int nS;
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

    public String getLembda() {
        return lembda;
    }

    public int getnS() {
        return nS;
    }

    public int getST() {
        return ST;
    }

    public ProjectionList() {
        this.nS = 0;
        this.ST = 0;
    }

    public void inc_nS()
    {
        this.nS = this.nS +1;
    }

    public void inc_st()
    {
        this.ST = this.ST+1;
    }

    public void setProjection(ArrayList<String> projection) {
        this.projection = projection;
    }

    public void setLembda(String lembda) {
        this.lembda = lembda;
    }

    public void setnS(int nS) {
        this.nS = nS;
    }

    public void setST(int ST) {
        this.ST = ST;
    }

    @Override
    public String toString() {
//        return "ProjectionList{" +
//                "projection=" + projection +
//                ", lembda='" + lembda + '\'' +
//                ", nS=" + nS +
//                '}';

        return "(" + lembda + ", " + projection +"),  " + nS;

    }

}


