package crystalc;

public class ds_XIC  implements Comparable<ds_XIC>
{
    public double Mz = 0;
    public double Rt = 0;
    public double Height = 0;
    public double Area = 0;
    public double StartRt = 999999;
    public double EndRt = 0;

    public int compareTo(ds_XIC xic)
    {
        return Double.compare(Mz, xic.Mz);
    }
}