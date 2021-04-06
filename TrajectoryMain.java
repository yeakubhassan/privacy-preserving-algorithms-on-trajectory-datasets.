import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Objects;
import java.util.zip.Adler32;

public class TrajectoryMain {
    private static double pbr = 0.5;
    private static String csvFile = "C:/Study/Security and Privacy/Project/Code/data/data_c_100_7_200.csv";
    private static int adversary_number = 10;
    private static int adversary_location = 10;

    public static void main(String args[]) {
        ArrayList<ArrayList<String>> trajectories = new ArrayList<ArrayList<String>>();
        trajectories = getData();
        ArrayList<ArrayList<String>> adversaries = new ArrayList<ArrayList<String>>();
        adversaries = defineAdversaries();

        //newAlgorithm(trajectories,adversaries);
        //gsup(trajectories,adversaries);
        //ArrayList<Adversary> projection_tree = construct_tree(trajectories,adversaries);
        long startTime = System.nanoTime();
        //ThreatObjects threatObjects = new ThreatObjects();
        //threatObjects = threatID(trajectories,adversaries);
        newAlgorithm(trajectories,adversaries);
        long endTime = System.nanoTime();
        long total_time = endTime - startTime;
        System.out.println("New Algorithm total time: "+ (double)total_time/1000000000.00);
        //System.out.println("Tree threat time: "+ (double)total_time/1000000000.00 );

        trajectories.clear();
        trajectories = getData();
        startTime = System.nanoTime();
        gsup(trajectories,adversaries);
        endTime = System.nanoTime();
        total_time = endTime - startTime;
        System.out.println("Gsup total time: "+ (double)total_time/1000000000.00);

        /*
        startTime = System.nanoTime();
        //gsup_ptt(trajectories,adversaries);
        gsup_ptt_update(trajectories,adversaries);
        endTime = System.nanoTime();
        total_time = endTime - startTime;
        System.out.println("Gsup_PTT update total time: "+ (double)total_time/1000000000.00);

        //gsup(trajectories,adversaries);*/
        //gsup_ptt_update(trajectories,adversaries);
        //newAlgorithm(trajectories,adversaries);

    }

    //*********************************************New algorithm************************************************************

    public static ArrayList<Adversary> construct_tree(ArrayList<ArrayList<String>> trajectories, ArrayList<ArrayList<String>> adversaries)
    {
        ArrayList<Adversary> projection_tree = new ArrayList<Adversary>();
        for(int i=0;i<adversary_number;i++)
        {
            Adversary adversary = new Adversary();
            adversary.setId(i);
            projection_tree.add(adversary);
        }
        for(int i=0;i<trajectories.size();i++)
        {
            for(int j=0;j<adversary_number;j++) {
                ArrayList<String> projection = new ArrayList<String>();
                ArrayList<String> not_in_projection = new ArrayList<String>();

                for (int k = 0; k < trajectories.get(i).size(); k++) {
                    if (adversaries.get(j).contains(trajectories.get(i).get(k))) {
                        projection.add(trajectories.get(i).get(k));
                    } else {
                        not_in_projection.add(trajectories.get(i).get(k));
                    }
                }
                int pos = checkProjection(projection_tree.get(j), projection); ///check already a child or not
                Projection p;
                if (pos == -1) //not a child
                {
                    p = new Projection();
                    p.setProj(projection);
                    p.add_trajectory(i);
                } else {
                    p = projection_tree.get(j).getProjections().get(pos);
                    p.add_trajectory(i);
                }
                //add childen from not in trajectory
                for (int k = 0; k < not_in_projection.size(); k++)
                {
                    Projection pc = new Projection();
                    ArrayList<String> st = new ArrayList<String>();
                    st.add(not_in_projection.get(k));
                    int index = checkProjectionChild(p, st);//check if the child exists
                    if (index == -1)//a new child
                    {
                        pc.setProj(st);
                        pc.add_trajectory(i);
                        p.add_children(pc);
                    }
                    else
                    {
                        if(!p.getChildren().get(index).getTrajectories().contains(i))
                        {
                            p.getChildren().get(index).add_trajectory(i);
                        }
                    }
                }
                if (pos == -1) //not a child
                {
                  projection_tree.get(j).getProjections().add(p);
                } else {
                    projection_tree.get(j).getProjections().set(pos,p);
                }
            }
        }
        return projection_tree;
    }

    public static void printProjectionTree(ArrayList<Adversary> projection_tree)
    {
        System.out.println("Printing Projection Tree");
        for(int i=0;i<projection_tree.size();i++)
        {
            System.out.println("Adversary = "+i);
            for(int j =0;j<projection_tree.get(i).getProjections().size();j++)
            {
                System.out.print(projection_tree.get(i).getProjections().get(j).getProj() + " " + projection_tree.get(i).getProjections().get(j).isProblem() + " " + projection_tree.get(i).getProjections().get(j).getTrajectories() +" --> ");
                for(int k=0;k<projection_tree.get(i).getProjections().get(j).getChildren().size();k++)
                {
                    System.out.print(projection_tree.get(i).getProjections().get(j).getChildren().get(k).getProj() + " " + projection_tree.get(i).getProjections().get(j).getChildren().get(k).getTrajectories() +" , ");
                }
                System.out.println();
            }
            System.out.println();
            System.out.println();
        }
        System.out.println("End Printing Projection Tree");
    }

    public static int initial_threat(ArrayList<Adversary> projection_tree)
    {
        int N=0;
        for(int i=0;i<projection_tree.size();i++)
        {
            for(int j=0;j<projection_tree.get(i).getProjections().size();j++)
            {
                int f=0;
                for(int k=0;k<projection_tree.get(i).getProjections().get(j).getChildren().size();k++)
                {
                    double pr = (double) projection_tree.get(i).getProjections().get(j).getChildren().get(k).getTrajectories().size() / (double)projection_tree.get(i).getProjections().get(j).getTrajectories().size();
                    if(pr > pbr)
                    {
                        f=1;
                        projection_tree.get(i).getProjections().get(j).setProblem(true);
                        N = N + projection_tree.get(i).getProjections().get(j).getChildren().get(k).getTrajectories().size();
                    }
                }
                if(f==0)
                {
                    projection_tree.get(i).getProjections().get(j).setProblem(false);
                }
            }
        }
        return N;
    }

    public static void newAlgorithm(ArrayList<ArrayList<String>> trajectories, ArrayList<ArrayList<String>> adversaries)
    {
        ArrayList<Adversary> projection_tree = construct_tree(trajectories,adversaries);
        //printProjectionTree(projection_tree);
        int N = initial_threat(projection_tree);
        //System.out.println(N);
        while(N>0)
        {
            double max_ugain = -9999999;
            int mtr_index = -1;
            int mtR_index = -1;
            int ad = -1;
            ArrayList<String>suppressed = new ArrayList<String>();
            ArrayList<String>mtR = new ArrayList<String>();
            ArrayList<String>mtr = new ArrayList<String>();
            for(int i=0;i<projection_tree.size();i++)
            {
                for(int j=0;j<projection_tree.get(i).getProjections().size();j++)
                {
                    suppressed.clear();
                    if(projection_tree.get(i).getProjections().get(j).isProblem())
                    {
                        mtR = projection_tree.get(i).getProjections().get(j).getProj();
                        int mtR_size = projection_tree.get(i).getProjections().get(j).getProj().size();
                        for(int k=0;k<projection_tree.get(i).getProjections().size();k++)
                        {
                            if(k!=j && mtR.size() > projection_tree.get(i).getProjections().get(k).getProj().size())
                            {
                                int flag = 1;
                                for(int a=0;a<projection_tree.get(i).getProjections().get(k).getProj().size();a++)
                                {
                                    if(!mtR.contains(projection_tree.get(i).getProjections().get(k).getProj().get(a)))
                                    {
                                        flag = 0;
                                        break;
                                    }
                                }
                                if(flag == 1)
                                {
                                    for(int l =0;l<projection_tree.get(i).getProjections().get(j).getProj().size();l++)
                                    {
                                        if(!projection_tree.get(i).getProjections().get(k).getProj().contains(projection_tree.get(i).getProjections().get(j).getProj().get(l)))
                                        {
                                          suppressed.add(projection_tree.get(i).getProjections().get(j).getProj().get(l));
                                        }
                                    }
                                    int mtr_size = projection_tree.get(i).getProjections().get(k).getProj().size();
                                    mtr = projection_tree.get(i).getProjections().get(k).getProj();
                                   // double g = getUgain(N,projection_tree,trajectories,mtR_index,mtr_index,ad,mtr_size,mtR_size);
                                    double g = ugian(trajectories,trajectories,adversaries,mtR, mtr,N, i);
                                    //System.out.println(mtR + " " + mtr + " " +g);
                                    if(g >= max_ugain )
                                    {
                                        if(g == max_ugain && i < ad)
                                        {
                                            max_ugain = g;
                                            mtr_index = k;
                                            mtR_index = j;
                                            ad = i;
                                        }
                                        else if(g > max_ugain)
                                        {
                                            max_ugain = g;
                                            mtr_index = k;
                                            mtR_index = j;
                                            ad = i;
                                        }

                                    }
                                }
                            }
                        }
                        //double g = getUgain(N,projection_tree,trajectories,j,-1,i,0,mtR_size);
                        mtr = new ArrayList<String>();
                        double g = ugian(trajectories,trajectories,adversaries,mtR, mtr,N, i);
                        //System.out.println(mtR + " " + mtr + " " +g);
                        if(g >= max_ugain)
                        {
                            if(g == max_ugain && i < ad)
                            {
                                max_ugain = g;
                                mtr_index = -1;
                                mtR_index = j;
                                ad = i;
                                for(int l =0;l<projection_tree.get(i).getProjections().get(j).getProj().size();l++)
                                {
                                    suppressed.add(projection_tree.get(i).getProjections().get(j).getProj().get(l));
                                }
                            }
                            else if(g > max_ugain)
                            {
                                max_ugain = g;
                                mtr_index = -1;
                                mtR_index = j;
                                ad = i;
                                for(int l =0;l<projection_tree.get(i).getProjections().get(j).getProj().size();l++)
                                {
                                    suppressed.add(projection_tree.get(i).getProjections().get(j).getProj().get(l));
                                }
                            }

                        }
                    }
                }
            }
            if(max_ugain != -9999999)
            {
                //System.out.println(projection_tree.get(ad).getProjections().get(mtR_index).getProj());
                for(int i=0;i<projection_tree.get(ad).getProjections().get(mtR_index).getTrajectories().size();i++)
                {
                    ArrayList<String>t = new ArrayList<String>();
                    for(int j=0;j<trajectories.get(projection_tree.get(ad).getProjections().get(mtR_index).getTrajectories().get(i)).size();j++)
                    {
                        if(!suppressed.contains(trajectories.get(projection_tree.get(ad).getProjections().get(mtR_index).getTrajectories().get(i)).get(j)))
                        {
                            t.add(trajectories.get(projection_tree.get(ad).getProjections().get(mtR_index).getTrajectories().get(i)).get(j));
                        }
                    }
                    trajectories.set(projection_tree.get(ad).getProjections().get(mtR_index).getTrajectories().get(i),t);
                }
                if(mtr_index != -1)
                {
                    updateTree(projection_tree,mtR_index,mtr_index,ad);
                }
                else
                {
                    updateTreeVoid(projection_tree,mtR_index,ad);
                }
                int temp = initial_threat(projection_tree);
                //System.out.println(N + " " +temp );
                N = temp;
                //System.out.println("Initial Threat : " + N);
                //printData(trajectories);
                //break;
            }
            else
            {
                break;
            }
        }
        //printProjectionTree(projection_tree);
    }

    public static double getUgain(int N, ArrayList<Adversary> projection_tree, ArrayList<ArrayList<String>> trajectories, int mtR, int mtr, int ad,  int mtr_size, int mtR_size)
    {
        double gain=0;
        double ploss = 0;
        ArrayList<Adversary> new_tree= new  ArrayList<Adversary>();

        for(int i=0;i<projection_tree.size();i++)
        {
            Adversary add = new Adversary();
            ArrayList<Projection> projections = new ArrayList<Projection>();
            for(int j=0;j<projection_tree.get(i).getProjections().size();j++)
            {
                Projection p = new Projection();
                p.setProblem(projection_tree.get(i).getProjections().get(j).isProblem());
                ArrayList<String> ppp = new ArrayList<String>();
                for(int k=0;k<projection_tree.get(i).getProjections().get(j).getProj().size();k++)
                {
                    ppp.add(projection_tree.get(i).getProjections().get(j).getProj().get(k));
                }
                p.setProj(ppp);
                for(int k=0;k<projection_tree.get(i).getProjections().get(j).getTrajectories().size();k++)
                {
                    p.add_trajectory(projection_tree.get(i).getProjections().get(j).getTrajectories().get(k));

                }
                for(int k=0;k<projection_tree.get(i).getProjections().get(j).getChildren().size();k++)
                {
                    Projection pc = new Projection();
                    ppp = new ArrayList<String>();
                    for(int l=0;l<projection_tree.get(i).getProjections().get(j).getChildren().get(k).getProj().size();l++)
                    {
                        ppp.add(projection_tree.get(i).getProjections().get(j).getChildren().get(k).getProj().get(l));
                    }
                    pc.setProj(ppp);
                    for(int l=0;l<projection_tree.get(i).getProjections().get(j).getChildren().get(k).getTrajectories().size();l++)
                    {
                        pc.add_trajectory(projection_tree.get(i).getProjections().get(j).getChildren().get(k).getTrajectories().get(l));
                    }
                    p.add_children(pc);
                }
                projections.add(p);
            }
            add.setProjections(projections);
            new_tree.add(add);
        }
        //printProjectionTree(projection_tree);
        //printProjectionTree(new_tree);
        //new_tree = projection_tree;
        for(int i=0;i<projection_tree.get(ad).getProjections().get(mtR).getTrajectories().size();i++)
        {
            double t = trajectories.get(projection_tree.get(ad).getProjections().get(mtR).getTrajectories().get(i)).size();
            double t_prime = t - mtR_size +  mtr_size;
            ploss = ploss + (1 - ((t_prime*(t_prime-1))/(t*(t-1))));
        }
        if(mtr!=-1)
        {
            updateTree(new_tree,mtR,mtr,ad);
        }else
        {
            updateTreeVoid(new_tree,mtR,ad);
        }
        int N_prime = initial_threat(new_tree);
        gain = (((double)N - (double)N_prime)/(double)N)*(1/ploss);
        return gain;
    }

    public static  Adversary deepClone(Adversary object){
        try {
            ByteArrayOutputStream byteArrayOutputStream = new ByteArrayOutputStream();
            ObjectOutputStream objectOutputStream = new ObjectOutputStream(byteArrayOutputStream);
            objectOutputStream.writeObject(object);
            ByteArrayInputStream bais = new ByteArrayInputStream(byteArrayOutputStream.toByteArray());
            ObjectInputStream objectInputStream = new ObjectInputStream(bais);
            return (Adversary) objectInputStream.readObject();
        }
        catch (Exception e) {
            e.printStackTrace();
            return null;
        }
    }

    public static void updateTreeVoid(ArrayList<Adversary> projection_tree, int mtR, int ad)
    {
        //trajectories that contains mtR, update their children due to the suppression
        for(int i=0;i<projection_tree.size();i++)
        {
            if(i!=ad)
            {
                for(int j=0;j<projection_tree.get(i).getProjections().size();j++)
                {
                    for(int k = 0;k<projection_tree.get(ad).getProjections().get(mtR).getTrajectories().size();k++)
                    {
                        if(projection_tree.get(i).getProjections().get(j).getTrajectories().contains(projection_tree.get(ad).getProjections().get(mtR).getTrajectories().get(k)))
                        {
                            ArrayList<Integer> empty_node = new ArrayList<Integer>();
                            for(int l=0;l<projection_tree.get(i).getProjections().get(j).getChildren().size();l++)
                            {
                                if(projection_tree.get(ad).getProjections().get(mtR).getProj().contains(projection_tree.get(i).getProjections().get(j).getChildren().get(l).getProj().get(0)))
                                {
                                    if(projection_tree.get(i).getProjections().get(j).getChildren().get(l).getTrajectories().contains(projection_tree.get(ad).getProjections().get(mtR).getTrajectories().get(k)))
                                    {
                                        projection_tree.get(i).getProjections().get(j).getChildren().get(l).remove_trajectories(projection_tree.get(ad).getProjections().get(mtR).getTrajectories().get(k));
                                    }
                                }
                                if(projection_tree.get(i).getProjections().get(j).getChildren().get(l).getTrajectories().size()==0)
                                {
                                    empty_node.add(l);
                                }
                            }
                            for(int l=0;l<empty_node.size();l++)
                            {
                                projection_tree.get(i).getProjections().get(j).remove_child(empty_node.get(l));
                            }
                        }
                    }
                }
            }
        }
        projection_tree.get(ad).remove_projection(mtR);
    }

    public static void updateTree(ArrayList<Adversary> projection_tree, int mtR, int mtr,int ad)
    {
        ArrayList<String>suppressed = new ArrayList<String>();
        for(int i=0;i<projection_tree.get(ad).getProjections().get(mtR).getProj().size();i++)
        {
            if(!projection_tree.get(ad).getProjections().get(mtr).getProj().contains(projection_tree.get(ad).getProjections().get(mtR).getProj().get(i)))
            {
                suppressed.add(projection_tree.get(ad).getProjections().get(mtR).getProj().get(i));
            }
        }
        //add trajectories of mtR to mtr
        for(int i=0;i<projection_tree.get(ad).getProjections().get(mtR).getTrajectories().size();i++)
        {
            if(!projection_tree.get(ad).getProjections().get(mtr).getTrajectories().contains(projection_tree.get(ad).getProjections().get(mtR).getTrajectories().get(i)))
            {
                projection_tree.get(ad).getProjections().get(mtr).add_trajectory(projection_tree.get(ad).getProjections().get(mtR).getTrajectories().get(i));
            }
        }
        //update mtR children's to mtr children
        for(int i=0;i<projection_tree.get(ad).getProjections().get(mtR).getChildren().size();i++)
        {
            int pos = checkProjectionChild(projection_tree.get(ad).getProjections().get(mtr), projection_tree.get(ad).getProjections().get(mtR).getChildren().get(i).getProj());
            if(pos == -1)
            {
                Projection p = new Projection();
                p = projection_tree.get(ad).getProjections().get(mtR).getChildren().get(i);
                projection_tree.get(ad).getProjections().get(mtr).add_children(p);
            }
            else
            {
                for(int j=0;j<projection_tree.get(ad).getProjections().get(mtR).getChildren().get(i).getTrajectories().size();j++)
                {
                    if(!projection_tree.get(ad).getProjections().get(mtr).getChildren().get(pos).getTrajectories().contains(projection_tree.get(ad).getProjections().get(mtR).getChildren().get(i).getTrajectories().get(j)))
                    {
                        projection_tree.get(ad).getProjections().get(mtr).getChildren().get(pos).add_trajectory(projection_tree.get(ad).getProjections().get(mtR).getChildren().get(i).getTrajectories().get(j));
                    }
                }

            }
        }
        //trajectories that contains mtR, update their children due to the suppression
        for(int i=0;i<projection_tree.size();i++)
        {
            if(i!=ad)
            {
                for(int j=0;j<projection_tree.get(i).getProjections().size();j++)
                {
                    for(int k = 0;k<projection_tree.get(ad).getProjections().get(mtR).getTrajectories().size();k++)
                    {
                        if(projection_tree.get(i).getProjections().get(j).getTrajectories().contains(projection_tree.get(ad).getProjections().get(mtR).getTrajectories().get(k)))
                        {
                            ArrayList<Integer> empty_node = new ArrayList<Integer>();
                            for(int l=0;l<projection_tree.get(i).getProjections().get(j).getChildren().size();l++)
                            {
                                if(suppressed.contains(projection_tree.get(i).getProjections().get(j).getChildren().get(l).getProj().get(0)))
                                {
                                   if(projection_tree.get(i).getProjections().get(j).getChildren().get(l).getTrajectories().contains(projection_tree.get(ad).getProjections().get(mtR).getTrajectories().get(k)))
                                   {
                                       projection_tree.get(i).getProjections().get(j).getChildren().get(l).remove_trajectories(projection_tree.get(ad).getProjections().get(mtR).getTrajectories().get(k));
                                   }
                                }
                                if(projection_tree.get(i).getProjections().get(j).getChildren().get(l).getTrajectories().size()==0)
                                {
                                    empty_node.add(l);
                                }
                            }
                            for(int l=0;l<empty_node.size();l++)
                            {
                                //projection_tree.get(i).getProjections().get(j).remove_child(empty_node.get(l));
                            }
                        }
                    }
                }
            }
        }
        ///delete mtR
        projection_tree.get(ad).remove_projection(mtR);
    }

    public static int checkProjection(Adversary ad, ArrayList<String> projection)
    {
        int pos = -1;
        for(int i=0;i<ad.getProjections().size();i++)
        {
            if(ad.getProjections().get(i).getProj().equals(projection))
            {
                pos = i;
            }
        }
        return pos;
    }

    public static int checkProjectionChild(Projection p, ArrayList<String> projection)
    {
        int pos = -1;
        for(int i=0;i<p.getChildren().size();i++)
        {
           if(p.getChildren().get(i).getProj().equals(projection))
           {
               pos = i;
               break;
           }
        }
        return  pos;
    }

    public static int getAdversaryNumber(String st,ArrayList<ArrayList<String>> adversaries)
    {
        int f=-1;
        for(int i=0;i<adversaries.size();i++)
        {
            if(adversaries.get(i).contains(st))
            {
                f=i;
                break;
            }
        }
        return f;
    }


    //*******************************************************************************************************

    public static void gsup_ptt_update(ArrayList<ArrayList<String>> trajectories, ArrayList<ArrayList<String>> adversaries)
    {
        ArrayList<ptt> ptt_nodes = new ArrayList<ptt>();
        for(int i=0;i<trajectories.size();i++) {
            for (int j = 0; j < adversaries.size(); j++) {
                ptt node = new ptt();
                ArrayList<String> proj = new ArrayList<String>();
                ArrayList<String> not_in_proj = new ArrayList<String>();
                for (int k = 0; k < trajectories.get(i).size(); k++) {
                    if (adversaries.get(j).contains(trajectories.get(i).get(k))) {
                        proj.add(trajectories.get(i).get(k));
                    } else {
                        not_in_proj.add(trajectories.get(i).get(k));
                    }
                }
                if(proj.size()>0)
                {
                    int pos = check_ptt_first_level(ptt_nodes,proj);
                    if( pos == -1)
                    {
                        node.setProjection(proj);
                        node.add_trajectory(i);
                    }
                    else
                    {
                        ptt n = ptt_nodes.get(pos);
                        n.add_trajectory(i);
                        ptt_nodes.set(pos,n);
                        node = ptt_nodes.get(pos);
                    }
                    for(int x=0;x<not_in_proj.size();x++) {

                        ptt c = new ptt();
                        ArrayList<String> p = new ArrayList<String>();
                        p.add(not_in_proj.get(x));
                        int poss = check_ptt_children(node, p);
                        if (poss == -1) {
                            c.setProjection(p);
                            c.add_trajectory(i);
                            node.add_children(c);
                        } else
                        {
                            c = node.getChildren().get(poss);
                            c.add_trajectory(i);
                            node.getChildren().set(poss,c);
                        }
                    }
                    if( pos == -1)
                    {
                        ptt_nodes.add(node);
                    }
                    else
                    {
                        ptt_nodes.set(pos, node);
                    }
                }

            }
        }
        ThreatObjects threatObjects = pttThreatID(ptt_nodes,adversaries);
        int N = threatObjects.getNumber_of_problems();
        System.out.println(N);
        //print_ptt(ptt_nodes);
        while(N>0)
        {
            ArrayList<String> mtR = new ArrayList<String>();
            ArrayList<String> mtr = new ArrayList<String>();
            double max_ugain = -9999999;
            int adversary_ugain = -1;
            for (int i = 0; i < threatObjects.getProjectionLists().size(); i++) {
                ArrayList<String> tR = new ArrayList<String>();
                tR = threatObjects.getProjectionLists().get(i).getProjection();
                int tR_adversary = threatObjects.getProjectionLists().get(i).getAdversary();
                ArrayList<String> tr = new ArrayList<String>();
                int f=0;
                double Ugain;

                Ugain = ugian(trajectories, trajectories, adversaries, tR, tr, N, tR_adversary);
                if (Ugain > max_ugain) {
                    max_ugain = Ugain;
                    mtR = tR;
                    mtr = tr;
                    adversary_ugain = tR_adversary;
                }
                for (int j = 0; j < threatObjects.getProjectionLists().size(); j++) {
                    if (j != i) {
                        tr = new ArrayList<String>();
                        for (int k = 0; k < threatObjects.getProjectionLists().get(j).getProjection().size(); k++) {
                            if (tR.contains(threatObjects.getProjectionLists().get(j).getProjection().get(k))) {
                                tr.add(threatObjects.getProjectionLists().get(j).getProjection().get(k));
                            }
                        }
                        if (tr.size() < tR.size() && !tr.isEmpty()) {
                            Ugain = ugian(trajectories, trajectories, adversaries, tR, tr, N, tR_adversary);
                            if (Ugain > max_ugain) {
                                max_ugain = Ugain;
                                mtR = tR;
                                mtr = tr;
                                adversary_ugain = tR_adversary;
                            }
                        }
                    }
                }
            }
            if(max_ugain == -9999999)break;
            int mtR_index = -1;
            int mtr_index = -1;

            for(int i=0;i<ptt_nodes.size();i++)
            {
                if(check_projction(ptt_nodes.get(i).getProjection(),mtR))
                {
                    mtR_index = i;
                    break;
                }
            }
            System.out.println(mtR + " " + mtr);
            for(int i=0;i<ptt_nodes.size();i++)
            {
                if(check_projction(ptt_nodes.get(i).getProjection(),mtr))
                {
                    mtr_index = i;
                    break;
                }
            }

            ptt node = new ptt();
            if(mtr_index != -1)
            {
                node = ptt_nodes.get(mtr_index);
            }
            else
            {
                if(mtR_index != -1)
                {
                    for(int i=0;i<ptt_nodes.get(mtR_index).getTrajectories().size();i++)
                    {
                        for(int j=0;j<ptt_nodes.size();j++)
                        {
                            if(j!=mtR_index)
                            {
                                if(ptt_nodes.get(j).getTrajectories().contains(ptt_nodes.get(mtR_index).getTrajectories().get(i)))
                                {
                                    for(int k=0;k<ptt_nodes.get(j).getChildren().size();k++)
                                    {
                                        if(ptt_nodes.get(mtR_index).getProjection().contains(ptt_nodes.get(j).getChildren().get(k).getProjection().get(0)))
                                        {
                                            ptt_nodes.get(j).getChildren().get(k).remove_trajectories(ptt_nodes.get(mtR_index).getTrajectories().get(i));
                                            if(ptt_nodes.get(j).getChildren().get(k).getTrajectories().size() <= 0)
                                            {
                                                ptt_nodes.get(j).remove_child(k);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    ptt_nodes.remove(mtR_index);

                }
            }

            if(mtR_index != -1 && mtr_index != -1)
            {
                ArrayList<String> replace = new ArrayList<String>();
                for(int i=0;i<ptt_nodes.get(mtR_index).getProjection().size();i++)
                {
                    if(!ptt_nodes.get(mtr_index).getProjection().contains(ptt_nodes.get(mtR_index).getProjection().get(i)))
                    {
                        replace.add(ptt_nodes.get(mtR_index).getProjection().get(i));
                    }
                }
                for(int i=0;i<ptt_nodes.get(mtR_index).getTrajectories().size();i++)
                {
                    for(int j=0;j<ptt_nodes.size();j++)
                    {
                        if(j!=mtR_index)
                        {
                            if(ptt_nodes.get(j).getTrajectories().contains(ptt_nodes.get(mtR_index).getTrajectories().get(i)))
                            {
                                for(int k=0;k<ptt_nodes.get(j).getChildren().size();k++)
                                {
                                    if(ptt_nodes.get(mtR_index).getProjection().contains(ptt_nodes.get(j).getChildren().get(k).getProjection().get(0)) && replace.contains(ptt_nodes.get(j).getChildren().get(k).getProjection().get(0)))
                                    {
                                        ptt_nodes.get(j).getChildren().get(k).remove_trajectories(ptt_nodes.get(mtR_index).getTrajectories().get(i));
                                        if(ptt_nodes.get(j).getChildren().get(k).getTrajectories().size() <= 0)
                                        {
                                            ptt_nodes.get(j).remove_child(k);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                for(int i=0;i<ptt_nodes.get(mtR_index).getTrajectories().size();i++)
                {
                    if(!ptt_nodes.get(mtr_index).getTrajectories().contains(ptt_nodes.get(mtR_index).getTrajectories().get(i)))
                    {
                        node.add_trajectory(ptt_nodes.get(mtR_index).getTrajectories().get(i));
                    }
                }

                for(int i=0;i<ptt_nodes.get(mtR_index).getChildren().size();i++)
                {
                    int pos = check_ptt_children(ptt_nodes.get(mtr_index),ptt_nodes.get(mtR_index).getChildren().get(i).getProjection());
                    if(pos != -1)
                    {
                        for(int j=0;j<ptt_nodes.get(mtR_index).getChildren().get(i).getTrajectories().size();j++)
                        {
                            if(!ptt_nodes.get(mtr_index).getChildren().get(pos).getTrajectories().contains(ptt_nodes.get(mtR_index).getChildren().get(i).getTrajectories().get(j)))
                            {
                                node.getChildren().get(pos).add_trajectory(ptt_nodes.get(mtR_index).getChildren().get(i).getTrajectories().get(j));
                            }
                        }
                    }
                    else
                    {
                        ptt c_node = ptt_nodes.get(mtR_index).getChildren().get(i);
                        node.add_children(c_node);
                    }
                }
                ptt_nodes.set(mtr_index,node);
                ptt_nodes.remove(mtR_index);
            }

            for (int i = 0; i < trajectories.size(); i++) {
                ArrayList<String> proj = new ArrayList<String>();
                for (int j = 0; j < trajectories.get(i).size(); j++) {
                    if (adversaries.get(adversary_ugain).contains(trajectories.get(i).get(j))) {
                        proj.add(trajectories.get(i).get(j));
                    }
                }
                if (mtR.equals(proj)) {
                    ArrayList<String> sup = new ArrayList<String>();
                    for (int j = 0; j < mtR.size(); j++) {
                        if (!mtr.contains(mtR.get(j))) {
                            sup.add(mtR.get(j));
                        }

                    }
                    ArrayList<String> traj = new ArrayList<String>();
                    for (int j = 0; j < trajectories.get(i).size(); j++) {
                        if (!sup.contains(trajectories.get(i).get(j))) {
                            traj.add(trajectories.get(i).get(j));
                        }
                    }
                    trajectories.set(i, traj);
                }
            }

            threatObjects = pttThreatID(ptt_nodes,adversaries);
            N = threatObjects.getNumber_of_problems();

            System.out.println(N);
            print_ptt(ptt_nodes);
            //break;
        }
    }

    public static void gsup_ptt(ArrayList<ArrayList<String>> trajectories, ArrayList<ArrayList<String>> adversaries)
    {
        //PTT-creation
        ArrayList<ptt> ptt_nodes = new ArrayList<ptt>();
        for(int i=0;i<trajectories.size();i++) {
            for (int j = 0; j < adversaries.size(); j++) {
                ptt node = new ptt();
                ArrayList<String> proj = new ArrayList<String>();
                ArrayList<String> not_in_proj = new ArrayList<String>();
                for (int k = 0; k < trajectories.get(i).size(); k++) {
                    if (adversaries.get(j).contains(trajectories.get(i).get(k))) {
                        proj.add(trajectories.get(i).get(k));
                    } else {
                        not_in_proj.add(trajectories.get(i).get(k));
                    }
                }
                if(proj.size()>0)
                {
                    int pos = check_ptt_first_level(ptt_nodes,proj);
                    if( pos == -1)
                    {
                        node.setProjection(proj);
                        node.add_trajectory(i);
                    }
                    else
                    {
                        ptt n = ptt_nodes.get(pos);
                        n.add_trajectory(i);
                        ptt_nodes.set(pos,n);
                        node = ptt_nodes.get(pos);
                    }
                    for(int x=0;x<not_in_proj.size();x++) {

                        ptt c = new ptt();
                        ArrayList<String> p = new ArrayList<String>();
                        p.add(not_in_proj.get(x));
                        int poss = check_ptt_children(node, p);
                        if (poss == -1) {
                            c.setProjection(p);
                            c.add_trajectory(i);
                            node.add_children(c);
                        } else
                        {
                            c = node.getChildren().get(poss);
                            c.add_trajectory(i);
                            node.getChildren().set(poss,c);
                        }
                    }
                    if( pos == -1)
                    {
                        ptt_nodes.add(node);
                    }
                    else
                    {
                        ptt_nodes.set(pos, node);
                    }
                }

            }
        }
        //print_ptt(ptt_nodes);

        ThreatObjects threatObjects = pttThreatID(ptt_nodes,adversaries);
        int N = threatObjects.getNumber_of_problems();
        //System.out.println(N);
        //gsup_ptt
        while(N>0)
        {
            ArrayList<String> mtR = new ArrayList<String>();
            ArrayList<String> mtr = new ArrayList<String>();
            double max_ugain = -9999999;

            int adversary_ugain = -1;
            for (int i = 0; i < threatObjects.getProjectionLists().size(); i++) {
                ArrayList<String> tR = new ArrayList<String>();
                tR = threatObjects.getProjectionLists().get(i).getProjection();
                int tR_adversary = threatObjects.getProjectionLists().get(i).getAdversary();
                ArrayList<String> tr = new ArrayList<String>();
                int f=0;
                double Ugain;
                /*for (int a = 0; a < threatObjects.getProblematicPair().size(); a++) {
                    if (threatObjects.getProjectionLists().get(threatObjects.getProblematicPair().get(a)).getProjection().equals(tR)) {
                        f = 1;
                    }
                }*/

                Ugain = ugian(trajectories, trajectories, adversaries, tR, tr, N, tR_adversary);
                if (Ugain > max_ugain) {
                    max_ugain = Ugain;
                    mtR = tR;
                    mtr = tr;
                    adversary_ugain = tR_adversary;
                }
                for (int j = 0; j < threatObjects.getProjectionLists().size(); j++) {
                    if (j != i) {
                        tr = new ArrayList<String>();
                        for (int k = 0; k < threatObjects.getProjectionLists().get(j).getProjection().size(); k++) {
                            if (tR.contains(threatObjects.getProjectionLists().get(j).getProjection().get(k))) {
                                tr.add(threatObjects.getProjectionLists().get(j).getProjection().get(k));
                            }
                        }
                        if (tr.size() < tR.size() && !tr.isEmpty()) {
                            Ugain = ugian(trajectories, trajectories, adversaries, tR, tr, N, tR_adversary);
                            if (Ugain > max_ugain) {
                                max_ugain = Ugain;
                                mtR = tR;
                                mtr = tr;
                                adversary_ugain = tR_adversary;
                            }
                        }
                    }
                }
            }
            if(max_ugain == -9999999)break;
            int mtR_index = -1;
            int mtr_index = -1;

            for(int i=0;i<ptt_nodes.size();i++)
            {
                if(check_projction(ptt_nodes.get(i).getProjection(),mtR))
                {
                    mtR_index = i;
                    break;
                }
            }
            for(int i=0;i<ptt_nodes.size();i++)
            {
                if(check_projction(ptt_nodes.get(i).getProjection(),mtr))
                {
                    mtr_index = i;
                    break;
                }
            }
            ptt node = new ptt();
            if(mtr_index != -1)
            {
                node = ptt_nodes.get(mtr_index);

            }
            else
            {
                if(mtR_index != -1)
                {
                    ptt_nodes.remove(mtR_index);
                }
            }
            if(mtR_index != -1 && mtr_index != -1)
            {
                for(int i=0;i<ptt_nodes.get(mtR_index).getTrajectories().size();i++)
                {
                    if(!ptt_nodes.get(mtr_index).getTrajectories().contains(ptt_nodes.get(mtR_index).getTrajectories().get(i)))
                    {
                        node.add_trajectory(ptt_nodes.get(mtR_index).getTrajectories().get(i));
                    }
                }

                for(int i=0;i<ptt_nodes.get(mtR_index).getChildren().size();i++)
                {
                    int pos = check_ptt_children(ptt_nodes.get(mtr_index),ptt_nodes.get(mtR_index).getChildren().get(i).getProjection());
                    if(pos != -1)
                    {
                        for(int j=0;j<ptt_nodes.get(mtR_index).getChildren().get(i).getTrajectories().size();j++)
                        {
                                if(!ptt_nodes.get(mtr_index).getChildren().get(pos).getTrajectories().contains(ptt_nodes.get(mtR_index).getChildren().get(i).getTrajectories().get(j)))
                                {
                                        node.getChildren().get(pos).add_trajectory(ptt_nodes.get(mtR_index).getChildren().get(i).getTrajectories().get(j));
                                }
                        }
                    }
                    else
                    {
                        ptt c_node = ptt_nodes.get(mtR_index).getChildren().get(i);
                        node.add_children(c_node);
                    }
                }
                ptt_nodes.set(mtr_index,node);
                ptt_nodes.remove(mtR_index);
            }
            threatObjects = pttThreatID(ptt_nodes,adversaries);
            N = threatObjects.getNumber_of_problems();
            //System.out.println(N);
            //print_ptt(ptt_nodes);
            //break;
        }
    }

    private static ThreatObjects pttThreatID(ArrayList<ptt> ptt_nodes,ArrayList<ArrayList<String>> adversaries)
    {
        int number_of_problems =0;
        ArrayList<ProjectionList> projectionLists = new ArrayList<ProjectionList>();
        for(int i=0;i<ptt_nodes.size();i++)
        {
            ptt n = ptt_nodes.get(i);
            for(int j=0;j<ptt_nodes.get(i).getChildren().size();j++)
            {
                ptt c = ptt_nodes.get(i).getChildren().get(j);
                double pr = (double) c.getTrajectories().size() / (double)n.getTrajectories().size();
                if(pr > pbr)
                {
                    ProjectionList projectionList = new ProjectionList();
                    projectionList.setProjection(n.getProjection());
                    projectionList.setLembda(c.getProjection().get(0));
                    projectionList.setnS(c.getTrajectories().size());
                    projectionList.setST(n.getTrajectories().size());
                    String str = n.getProjection().get(0);
                    int ad = -1;
                    for(int x=0;x<adversaries.size();x++)
                    {
                        if(adversaries.get(x).contains(str))
                        {
                            ad = x;
                        }
                    }
                    projectionList.setAdversary(ad);
                    projectionLists.add(projectionList);
                    number_of_problems += c.getTrajectories().size();
                }
            }
        }

        ThreatObjects threatObjects = new ThreatObjects();
        threatObjects.setProjectionLists(projectionLists);
        threatObjects.setNumber_of_problems(number_of_problems);
        return threatObjects;
    }


    public static void print_ptt(ArrayList<ptt> ptt_nodes)
    {
        for(int i=0;i<ptt_nodes.size();i++)
        {
            System.out.print(ptt_nodes.get(i).getProjection() + "{" + ptt_nodes.get(i).getTrajectories() + "} -->> ");
            for(int j=0;j<ptt_nodes.get(i).getChildren().size();j++)
            {
                System.out.print(ptt_nodes.get(i).getChildren().get(j).getProjection() + "{" + ptt_nodes.get(i).getChildren().get(j).getTrajectories() + "},  ");
            }
            System.out.println();
        }
    }

    public static int check_ptt_children(ptt ptt_node, ArrayList<String>pro)
    {
        int pos =-1;
        for(int i=0; i<ptt_node.getChildren().size();i++)
        {
            ptt p = ptt_node.getChildren().get(i);
            if(check_projction(p.getProjection(),pro))
            {
                pos = i;
                break;
            }
        }
        return pos;
    }

    public static int check_ptt_first_level(ArrayList<ptt> ptt_nodes, ArrayList<String>pro)
    {
        int position = -1;
        for(int i=0;i<ptt_nodes.size();i++)
        {
            ArrayList<String>proj = ptt_nodes.get(i).getProjection();
            if(check_projction(proj,pro))
            {
                position = i;
                break;
            }
        }
        return position;
    }

    public static boolean check_projction(ArrayList<String>proj1, ArrayList<String>proj2)
    {
        if(proj1.size() != proj2.size())
        {
            return false;
        }else
        {
            int flag =0;
            for(int i=0;i<proj1.size();i++)
            {
                if(!proj1.get(i).equals(proj2.get(i)))
                {
                    flag = 1;
                    break;
                }
            }
            if(flag == 0) {return true;}
            else {return false;}
        }
    }

    public static void gsup(ArrayList<ArrayList<String>> trajectories, ArrayList<ArrayList<String>> adversaries)
    {
        ThreatObjects threatObjects = new ThreatObjects();
        threatObjects = threatID(trajectories,adversaries);
        int N = threatObjects.getNumber_of_problems();
        ArrayList<ArrayList<String>> old_trajectories = new ArrayList<ArrayList<String>>();
        for(int i=0;i<trajectories.size();i++)
        {
            old_trajectories.add(trajectories.get(i));
        }
        //System.out.println(N);
        while (N>0) {
            old_trajectories = new ArrayList<ArrayList<String>>();
            for(int i=0;i<trajectories.size();i++)
            {
                old_trajectories.add(trajectories.get(i));
            }
            //System.out.println("start");
            threatObjects = threatID(trajectories,adversaries);
            ArrayList<String> mtR = new ArrayList<String>();
            ArrayList<String> mtr = new ArrayList<String>();
            double max_ugain = -9999999;
            int adversary_ugain = -1;
            for (int i = 0; i < threatObjects.getAllProjection().size(); i++) {
                ArrayList<String> tR = new ArrayList<String>();
                tR = threatObjects.getAllProjection().get(i).getProjection();
                int tR_adversary = threatObjects.getAllProjection().get(i).getAdversary();
                ArrayList<String> tr = new ArrayList<String>();
                int f=0;
                double Ugain;
                for (int a = 0; a < threatObjects.getProblematicPair().size(); a++) {
                    if (threatObjects.getProjectionLists().get(threatObjects.getProblematicPair().get(a)).getProjection().equals(tR)) {
                        f = 1;
                    }
                }
                if(f==1)
                {
                    Ugain = ugian(trajectories, trajectories, adversaries, tR, tr, N, tR_adversary);
                    //System.out.println(tR + " " + tr + " " +Ugain);
                    if (Ugain >= max_ugain) {
                        if(Ugain == max_ugain & tR_adversary < adversary_ugain)
                        {
                            max_ugain = Ugain;
                            mtR = tR;
                            mtr = tr;
                            adversary_ugain = tR_adversary;
                        }
                        else if(Ugain > max_ugain)
                        {
                            max_ugain = Ugain;
                            mtR = tR;
                            mtr = tr;
                            adversary_ugain = tR_adversary;
                        }
                    }
                    for (int j = 0; j < threatObjects.getAllProjection().size(); j++) {
                        if (j != i) {
                            tr = new ArrayList<String>();
                            tr = threatObjects.getAllProjection().get(j).getProjection();
                            int fll = 0;
                            for(int l =0;l<tr.size();l++)
                            {
                                if(!tR.contains(tr.get(l)))
                                {
                                    fll=1;
                                }
                            }
                            if(fll==0)
                            {
                                Ugain = ugian(trajectories, trajectories, adversaries, tR, tr, N, tR_adversary);
                                //System.out.println(tR + " " + tr + " " +Ugain);
                                if (Ugain >= max_ugain) {
                                    if(Ugain == max_ugain & tR_adversary < adversary_ugain)
                                    {
                                        max_ugain = Ugain;
                                        mtR = tR;
                                        mtr = tr;
                                        adversary_ugain = tR_adversary;
                                    }
                                    else if(Ugain > max_ugain)
                                    {
                                        max_ugain = Ugain;
                                        mtR = tR;
                                        mtr = tr;
                                        adversary_ugain = tR_adversary;
                                    }

                                }
                            }
                        }
                    }
                }
                /*for (int j = 0; j < threatObjects.getAllProjection().size(); j++) {
                    if (j != i) {
                        tr = new ArrayList<String>();
                        for (int k = 0; k < threatObjects.getAllProjection().get(j).getProjection().size(); k++) {
                            if (tR.contains(threatObjects.getAllProjection().get(j).getProjection().get(k))) {
                                tr.add(threatObjects.getAllProjection().get(j).getProjection().get(k));
                            }
                        }
                        int flag = 0;
                        for (int a = 0; a < threatObjects.getProblematicPair().size(); a++) {
                            if (threatObjects.getProjectionLists().get(threatObjects.getProblematicPair().get(a)).getProjection().equals(tR) || threatObjects.getProjectionLists().get(threatObjects.getProblematicPair().get(a)).getProjection().equals(tr)) {
                                flag = 1;
                            }
                        }
                        if (tr.size() < tR.size() && !tr.isEmpty() && flag == 1) {
                            Ugain = ugian(trajectories, trajectories, adversaries, tR, tr, N, tR_adversary);
                            System.out.println(tR + " " + tr + " " +Ugain);
                            if (Ugain > max_ugain) {
                                max_ugain = Ugain;
                                mtR = tR;
                                mtr = tr;
                                adversary_ugain = tR_adversary;
                            }
                        }
                    }
                }*/
            }
            //System.out.println(max_ugain);
            if(max_ugain == -9999999)break;
            for (int i = 0; i < trajectories.size(); i++) {
                ArrayList<String> proj = new ArrayList<String>();
                for (int j = 0; j < trajectories.get(i).size(); j++) {
                    if (adversaries.get(adversary_ugain).contains(trajectories.get(i).get(j))) {
                        proj.add(trajectories.get(i).get(j));
                    }
                }
                if (mtR.equals(proj)) {
                    ArrayList<String> sup = new ArrayList<String>();
                    for (int j = 0; j < mtR.size(); j++) {
                        if (!mtr.contains(mtR.get(j))) {
                            sup.add(mtR.get(j));
                        }

                    }
                    ArrayList<String> traj = new ArrayList<String>();
                    for (int j = 0; j < trajectories.get(i).size(); j++) {
                        if (!sup.contains(trajectories.get(i).get(j))) {
                            traj.add(trajectories.get(i).get(j));
                        }
                    }
                    trajectories.set(i, traj);
                }
            }
            threatObjects = threatID(trajectories, adversaries);
            //System.out.println(mtR + " " + mtr);
            //System.out.println(N + " " + threatObjects.getNumber_of_problems());
            N = threatObjects.getNumber_of_problems();
            //printData(trajectories);
        }
    }


    private static ArrayList<ArrayList<String>> getData()
    {
        ArrayList<ArrayList<String>> data = new ArrayList<ArrayList<String>>();

        BufferedReader br = null;
        String line = "";
        String cvsSplitBy = ",";
        String newStr ="abcdefghijklmnopqrstuvwxyz";
        int ii =0;
        try {
            br = new BufferedReader(new FileReader(csvFile));
            while ((line = br.readLine()) != null) {
                String[] str = line.split(cvsSplitBy);
                ArrayList<String> t = new ArrayList<String>();
                //System.out.println(str.length);
                for (int i = 0; i < str.length; i++) {
                    t.add(str[i]);
                }
                data.add(t);
                ii++;
                if(ii>8000)break;
            }

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (br != null) {
                try {
                    br.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
        return data;
    }

    private static ArrayList<ArrayList<String>> defineAdversaries()
    {
        ArrayList<ArrayList<String>> ad = new ArrayList<ArrayList<String>>();
        String newStr ="abcdefghijklmnopqrstuvwxyz";

        for(int i=0;i<adversary_number;i++)
        {
            ArrayList<String> str = new ArrayList<String>();
            for(int j=0;j<adversary_location;j++)
            {
                String st = newStr.charAt(i) + "" + j + "";
                str.add(st);
            }
            ad.add(str);
        }

        return ad;
    }

    private static void printData(ArrayList<ArrayList<String>> trajectories)
    {
        for(int i=0;i<trajectories.size();i++)
        {
            for (int j=0;j<trajectories.get(i).size();j++)
            {
                System.out.print(trajectories.get(i).get(j)+" ");
            }
            System.out.println("");
        }
    }

    private static ThreatObjects threatID(ArrayList<ArrayList<String>> trajectories, ArrayList<ArrayList<String>> adversaries)
    {
        int number_of_problems =0;
        ArrayList<ArrayList<String>> nS = new ArrayList<ArrayList<String>>();
        ArrayList<ProjectionList> projectionLists = new ArrayList<ProjectionList>();
        ArrayList<ProjectionSupport>allProjection = new ArrayList<ProjectionSupport>();
        ArrayList<Integer> problematicPair = new ArrayList<Integer>();
        for(int i=0;i<trajectories.size();i++)
        {
            for(int j=0;j<adversaries.size();j++)
            {
                ArrayList<String> projection = new ArrayList<String>();
                ArrayList<String> not_in_projection = new ArrayList<String>();
                for(int k = 0;k<trajectories.get(i).size();k++)
                {
                    if(adversaries.get(j).contains(trajectories.get(i).get(k)))
                    {
                        projection.add(trajectories.get(i).get(k));
                    }else
                    {
                        not_in_projection.add(trajectories.get(i).get(k));
                    }
                }
                int position = findInAllProjection(allProjection, projection);
                if(position == -1)
                {
                    ProjectionSupport ps = new ProjectionSupport();
                    ps.setProjection(projection);
                    ps.inc_st();
                    ps.setAdversary(j);
                    allProjection.add(ps);
                }
                else
                {
                    ProjectionSupport ps  = new ProjectionSupport();
                    ps = allProjection.get(position);
                    ps.inc_st();
                    allProjection.set(position, ps);
                }
                for(int a=0;a<not_in_projection.size();a++) {
                    ProjectionList projectionList = new ProjectionList();
                    int pos = findProjection(projectionLists, projection, not_in_projection.get(a));
                    if (pos == -1) {
                        projectionList.setProjection(projection);
                        projectionList.setLembda(not_in_projection.get(a));
                        projectionList.inc_nS();
                        projectionList.setAdversary(j);
                        projectionLists.add(projectionList);

                    } else {
                        projectionList = projectionLists.get(pos);
                        projectionList.inc_nS();
                        projectionLists.set(pos, projectionList);
                    }
                }
            }
        }
        for(int i=0;i<projectionLists.size();i++)
        {
            for(int j=0;j<allProjection.size();j++)
            {
                if(projectionLists.get(i).getProjection().equals(allProjection.get(j).getProjection()))
                {
                    double prob = (double)projectionLists.get(i).getnS()/(double)allProjection.get(j).getST();
                    if(prob > pbr)
                    {
                        number_of_problems = number_of_problems + projectionLists.get(i).getnS();
                        problematicPair.add(i);

                    }
                }
            }
        }
        ThreatObjects threatObjects = new ThreatObjects();
        threatObjects.setAllProjection(allProjection);
        threatObjects.setProblematicPair(problematicPair);
        threatObjects.setProjectionLists(projectionLists);
        threatObjects.setNumber_of_problems(number_of_problems);

        return threatObjects;
    }

    public static double ugian (ArrayList<ArrayList<String>> original_trajectories, ArrayList<ArrayList<String>> trajectories,ArrayList<ArrayList<String>> adversaries,ArrayList<String> tR, ArrayList<String> tr, int N, int adversary)
    {
        double Ugain = 0;
        ArrayList<Integer> affectedTraj = new ArrayList<Integer>();
        ArrayList<ArrayList<String>> temp_trajectories = new ArrayList<ArrayList<String>>();
        for(int i=0;i<trajectories.size();i++)
        {
            temp_trajectories.add(trajectories.get(i));
        }
        for(int i=0;i<temp_trajectories.size();i++)
        {
            ArrayList<String> projection = new ArrayList<String>();
                for(int j=0;j<temp_trajectories.get(i).size();j++)
                {
                    if(adversaries.get(adversary).contains(temp_trajectories.get(i).get(j)))
                    {
                        projection.add(temp_trajectories.get(i).get(j));
                    }
                }
                if(tR.equals(projection))
                {
                    affectedTraj.add(i);

                    ArrayList<String> sup = new ArrayList<String>();
                    for(int j=0;j<tR.size();j++)
                    {
                        if(!tr.contains(tR.get(j)))
                        {
                            sup.add(tR.get(j));
                        }

                    }
                    ArrayList<String> traj = new ArrayList<String>();

                    for(int j=0;j<temp_trajectories.get(i).size();j++)
                    {
                        if(!sup.contains(temp_trajectories.get(i).get(j)))
                        {
                            traj.add(temp_trajectories.get(i).get(j));
                        }
                    }
                    temp_trajectories.set(i,traj);
                }
        }
        ThreatObjects threatObjects = new ThreatObjects();
        threatObjects = threatID(temp_trajectories,adversaries);
        double N_prime = threatObjects.getNumber_of_problems();
        double ploss =0;
        for(int i=0;i<affectedTraj.size();i++)
        {
            double p =0;
            double t = original_trajectories.get(affectedTraj.get(i)).size();
            double t_prime = temp_trajectories.get(affectedTraj.get(i)).size();
            p = 1 - ((t_prime * (t_prime - 1))/(t * (t-1)));
            ploss = ploss + p;
        }
        Ugain =((N - N_prime) / N) * (1/ploss);
        return Ugain;
    }

    private static int findProjection( ArrayList<ProjectionList> projectionLists,ArrayList<String> proj, String lembda)
    {
        int pos = -1;
        for(int i=0;i<projectionLists.size();i++)
        {
            ArrayList<String> p = projectionLists.get(i).getProjection();
            String lem = projectionLists.get(i).getLembda();
            if(p.equals(proj) && lem.equals(lembda))
            {
               pos = i;
                //break;
            }
        }
        return pos;
    }

    private static int findInAllProjection(ArrayList<ProjectionSupport>allProjection, ArrayList<String> proj)
    {
        int pos = -1;
        for(int i=0;i<allProjection.size();i++)
        {
            ArrayList<String> p = allProjection.get(i).getProjection();
            if(p.equals(proj))
            {
                pos = i;
                //break;
            }
        }
        return pos;
    }

}
