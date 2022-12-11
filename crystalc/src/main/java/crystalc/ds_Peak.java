package crystalc;

public class ds_Peak implements Comparable<ds_Peak>
{
    public double mz = 0;
    public double rt = 0;
    public double intensity = 0;
    public double noiselevel = 0;
    public int scanNum = -1;

    public int compareTo(ds_Peak p)
    {
        return Double.compare(mz, p.mz);
    }

}