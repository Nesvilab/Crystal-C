package crystalc;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class p_Operator
{
    public static double CalculateMz(double PepMass, double Z)
    {
        double CalMz = (PepMass + Z * 1.00794) / Z;
        return CalMz;
    }

    public static double CalculatePepMonoMass(String PepSeq, Map<String, Double> AAMonoMassMap)
    {
        double PepMonoMass = 0f;

        for(int i=0;i<PepSeq.length();i++)
        {
            PepMonoMass += AAMonoMassMap.get(String.valueOf(PepSeq.charAt(i)));
        }

        return PepMonoMass;
    }

    public static List<Double> CalculateTheoreticalIsotopeMz(int isoNum, double Mz, double Z)
    { //Calculate isotope peaks, including the monoisotope peak
        List<Double> TheoreticalIsotopeLi = new ArrayList<Double>();

        for(int i=0; i<isoNum; i++)
        {
            double isoMZ = Mz + i*(1/ Z);
            TheoreticalIsotopeLi .add(isoMZ);
        }

        return TheoreticalIsotopeLi;
    }
}

