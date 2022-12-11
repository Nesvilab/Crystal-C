package crystalc;

import java.io.IOException;
import java.util.*;
import umich.ms.datatypes.scan.IScan;
import umich.ms.datatypes.scancollection.IScanCollection;
import umich.ms.datatypes.spectrum.ISpectrum;
import umich.ms.fileio.exceptions.FileParsingException;
import umich.ms.ims.isodistr.impl.AveragineTemplateProvider;
import umich.ms.ims.isodistr.api.ISpectrumTemplate;

public class p_ProcessPSM
{
    public static ds_IsotopeCluster CalculateMonoisotopeMass(int PrecursorScanNum, double TargetMz, ds_Parameters inp, IScanCollection Scans, AveragineTemplateProvider templates, int Z, double MzRange, double RtRange,  Double[] RtAry) throws FileParsingException, IOException
    {
        p_PredictMonoisotope pm = new p_PredictMonoisotope(PrecursorScanNum, TargetMz, Z,  inp.MassTol, MzRange, RtRange, Scans, RtAry);
        pm.GetRepresentative();
        ds_IsotopeCluster IsotopeCluster = pm.IsotopeCluster;

        double expMass = TargetMz*Z-Z;
        if(IsotopeCluster.RepMz > 0)
        {
            double topXicMass = (IsotopeCluster.RepMz - inp.aa.ProtonMass) * Z;
            if (!templates.checkMass(topXicMass)) {
                System.out.println("Warmming: Unreasonable mass in scan number: " + PrecursorScanNum + " (precursor M/Z: "+IsotopeCluster.RepMz+" ; z: " +Z+ ") , possibly due to incorrect charge states or precursor m/z" );
                IsotopeCluster.PredictMonoMass = topXicMass;
                IsotopeCluster.PredictMz = (IsotopeCluster.PredictMonoMass + Z * inp.aa.ProtonMass) / (double)Z;
            }
            else
            {
                ISpectrumTemplate template = templates.forMass(topXicMass);
                IsotopeCluster.PredictMonoMass = topXicMass + template.dmTopToMono();
                IsotopeCluster.PredictMz = (IsotopeCluster.PredictMonoMass + Z * inp.aa.ProtonMass) / (double)Z;
            }
        }

        return IsotopeCluster;
    }

    public static String CheckMissedCleavagePeptides(ds_Parameters inp, ds_PSM psm, double TargetMass)
    {
        String result = "";
        //Check missed cleavage peptides using ProteinAceNo
        result = FindMissedCleavagePeptides(inp.fastaMap, psm.ProteinAceNo, psm.PepSeq, psm.TheoPepMass, TargetMass, inp.aa.AAMonoMassMap, psm.MassAbsThreshold, inp.CleavageSites, inp.CleavageInhibitors);

        //Check missed cleavage peptides using AlternativeProteins
        if((result.isEmpty()) && (!psm.AlternativeProteins.isEmpty()))
        {
            double MinMassDif = 10000000;
            String[] AlterProteinAry = psm.AlternativeProteins.split(",");
            for(String apStr : AlterProteinAry)
            {
                String outStr = FindMissedCleavagePeptides(inp.fastaMap, apStr.trim(), psm.PepSeq, psm.TheoPepMass, TargetMass, inp.aa.AAMonoMassMap, psm.MassAbsThreshold, inp.CleavageSites, inp.CleavageInhibitors);
                if(!outStr.isEmpty())
                {
                    String[] tStrAry = outStr.split("[$\\@]");
                    double MassDif = Double.parseDouble(tStrAry[1]);
                    if(MassDif < MinMassDif)
                    {
                        MinMassDif = MassDif;
                        result = outStr;
                    }
                }
            }
        }

        return result;
    }

    private static String FindMissedCleavagePeptides(Map<String, String> fastaMap, String ProteinAceNo, String PepSeq, double TheoPepMass, double ExpPepMass, Map<String, Double> AAMonoMassMap, double MassAbsThreshold, String CleavageSites, String CleavageInhibitors)
    {
        String result = "";

        String ProtSeq = "";
        if(fastaMap.get(ProteinAceNo)!=null){
            ProtSeq = fastaMap.get(ProteinAceNo);
        }
        else{
            System.out.println("Protein: "+ProteinAceNo+" not found in the protein sequence database. Please check that you are using the same database as the search engine.");
            System.exit(1);
        }
        if(ProtSeq.contains("*")){
            System.out.println("Warning: The protein sequences in your fasta file include '*'. Please remove it before running Crystal-C.");
            ProtSeq=ProtSeq.replace("*","");
        }

        int Index = ProtSeq.indexOf(PepSeq);
        if(Index < 0) //check I/L
        {	//Replace I/L with X
            String mProtSeq = ProtSeq.replaceAll("I", "X");
            mProtSeq = mProtSeq.replaceAll("L", "X");
            String mTargetPepSeq = PepSeq.replaceAll("I", "X");
            mTargetPepSeq =  mTargetPepSeq.replaceAll("L",  "X");
            Index = mProtSeq.indexOf(mTargetPepSeq);
        }

        if (Index < 0)
        {
            System.out.println("Error. Can't find target peptides: "+PepSeq+"@"+ProteinAceNo+"   in the protein sequence.");
            System.exit(1);
        }

        int lIndex = (Index - 1) >= 0 ? (Index - 1) : 0;
        int rIndex = (Index + PepSeq.length()) < ProtSeq.length() ? (Index + PepSeq.length()) : ProtSeq.length();

        String rPepSeq = "";
        for(int i = rIndex; i < ProtSeq.length(); i++) //Get right-side peptides
        {
            rPepSeq += ProtSeq.charAt(i);
            if((i+1)<ProtSeq.length())
            {
                if(CleavageSites.contains(String.valueOf(ProtSeq.charAt(i))) && !CleavageInhibitors.contains(String.valueOf(ProtSeq.charAt(i+1))))
                {
                    break;
                }
            }
        }

        String lPepSeq = "";
        if(lIndex == 0)
        {
            lPepSeq = String.valueOf(ProtSeq.charAt(lIndex));
        }
        else
        {
            for(int i = lIndex - 1; i >= 0;i--) //Get left-side peptides
            {
                lPepSeq = ProtSeq.charAt(i+1) + lPepSeq;
                if(CleavageSites.contains(String.valueOf(ProtSeq.charAt(i))) && (!CleavageInhibitors.contains(String.valueOf(ProtSeq.charAt(i+1))) || (i+1 == lIndex)))
                {
                    break;
                }
            }
        }

        //Calculate delta mass
        double rPepMonoMass = !rPepSeq.isEmpty()? p_Operator.CalculatePepMonoMass(rPepSeq, AAMonoMassMap) : 0f;
        double lPepMonoMass = !lPepSeq.isEmpty()? p_Operator.CalculatePepMonoMass(lPepSeq, AAMonoMassMap) : 0f;

        //Check if there is at least one missed cleavage peptides
        double lMassDif = Math.abs(lPepMonoMass + TheoPepMass - ExpPepMass);
        double rMassDif = Math.abs(rPepMonoMass + TheoPepMass - ExpPepMass);

        //Find the minimum mass difference; if rMassDif = lMassDif, use the lMassDif/lpepSeq first.
        if((lMassDif <= MassAbsThreshold) && (lMassDif <= rMassDif))
        {
            result = lPepSeq + PepSeq + "@" + lMassDif + "@" + (lPepMonoMass + TheoPepMass) + "$";
        }
        else if ((rMassDif <= MassAbsThreshold) && (rMassDif <= lMassDif))
        {
            result = PepSeq + rPepSeq + "@" + rMassDif + "@" + (rPepMonoMass + TheoPepMass) + "$";
        }

        //List all possible peptides
        if(lMassDif <= MassAbsThreshold)
        {
            //Find the previous and next amino acid
            int IndexOfPreviousAA = ((lIndex - lPepSeq.length()) >= 0) ? (lIndex - lPepSeq.length()) : -1;
            int IndexOfNextAA = (rIndex < ProtSeq.length()) ? rIndex : -1;
            String PreviousAA = (IndexOfPreviousAA > 0) ? String.valueOf(ProtSeq.charAt(IndexOfPreviousAA)) : "-";
            String NextAA = (IndexOfNextAA > 0) ? String.valueOf(ProtSeq.charAt(IndexOfNextAA)) : "-";

            result += lPepSeq + PepSeq + "~" + lMassDif + "~" + (lPepMonoMass + TheoPepMass) + "~" + PreviousAA + "~" + NextAA + "~" + ProteinAceNo + "~True" + "&";
        }
        if(rMassDif <= MassAbsThreshold)
        {
            //Find the previous and next amino acid
            int IndexOfPreviousAA = (lIndex >= 0) ? lIndex : -1;
            int IndexOfNextAA = ((rIndex + rPepSeq.length()) < ProtSeq.length()) ? (rIndex + rPepSeq.length()) : -1;
            String PreviousAA = (IndexOfPreviousAA > 0) ? String.valueOf(ProtSeq.charAt(IndexOfPreviousAA)) : "-";
            String NextAA = (IndexOfNextAA > 0) ? String.valueOf(ProtSeq.charAt(IndexOfNextAA)) : "-";

            result += PepSeq + rPepSeq + "~" + rMassDif + "~" + (rPepMonoMass + TheoPepMass) + "~" + PreviousAA + "~" + NextAA + "~" + ProteinAceNo + "~True" + "&";
        }

        return result;
    }

    public static String CheckSemiTrypticPeptides(String PepSeq, double TheoPepMass, double ExpPepMass, double MassAbsThreshold, Map<String, Double> AAMonoMassMap)
    {
        String bestResult = "";
        String result = "";

        double MinLMassDif = 1000000;
        String lPepSeq = "";
        for(int i=1 ; i<= PepSeq.length() ; i++) //check from left side of sequence
        {
            double lMassDif = Math.abs(TheoPepMass - p_Operator.CalculatePepMonoMass(PepSeq.substring(0, i), AAMonoMassMap) - ExpPepMass);
            if(lMassDif <= MassAbsThreshold)
            {
                lPepSeq = PepSeq.substring(i, PepSeq.length());
                MinLMassDif = lMassDif;
                break;
            }
        }

        double MinRMassDif = 1000000;
        String rPepSeq = "";
        for(int i = PepSeq.length() ; i>=0 ; i--)
        {
            double rMassDif = Math.abs(TheoPepMass - p_Operator.CalculatePepMonoMass(PepSeq.substring(i, PepSeq.length()), AAMonoMassMap) - ExpPepMass);
            if(rMassDif <= MassAbsThreshold)
            {
                rPepSeq = PepSeq.substring(0, i);
                MinRMassDif = rMassDif;
                break;
            }
        }

        //Find the peptide sequence with the minimum delta mass
        if ((Math.abs(MinLMassDif) <= MassAbsThreshold) && (Math.abs(MinLMassDif) <= MinRMassDif))
        {
            //bestResult = lPepSeq + "@" + MinLMassDif + "@"  + p_Operator.CalculatePepMonoMass(lPepSeq, AAMonoMassMap);
            bestResult = lPepSeq + "@" + MinLMassDif + "@"  + (MinLMassDif+ExpPepMass);
        }
        else if (Math.abs(MinRMassDif) <= MassAbsThreshold)
        {
            //bestResult = rPepSeq + "@" + MinRMassDif + "@"  + p_Operator.CalculatePepMonoMass(rPepSeq, AAMonoMassMap);
            bestResult = rPepSeq + "@" + MinRMassDif + "@"  + (MinRMassDif+ExpPepMass);
        }

        //List all possible peptides
        if(!lPepSeq.isEmpty())
        {
            //Find the previous and next amino acid
            String PreviousAA = (lPepSeq.length() < PepSeq.length())? String.valueOf(PepSeq.charAt(PepSeq.length() - lPepSeq.length() - 1))  :"-";
            String NextAA = "-";

            //result += lPepSeq + "~" + MinLMassDif + "~" + p_Operator.CalculatePepMonoMass(lPepSeq, AAMonoMassMap) + "~" + PreviousAA + "~" + NextAA + "~-~false" + "&";
            result += lPepSeq + "~" + MinLMassDif + "~" + (MinLMassDif+ExpPepMass) + "~" + PreviousAA + "~" + NextAA + "~-~false" + "&";
        }
        if(!rPepSeq.isEmpty())
        {
            String PreviousAA = "-";
            String NextAA = (rPepSeq.length() < PepSeq.length()) ? String.valueOf(PepSeq.charAt(rPepSeq.length())) : "-";

            //result += rPepSeq + "~" + MinRMassDif +"~" + p_Operator.CalculatePepMonoMass(rPepSeq, AAMonoMassMap) + "~" + PreviousAA + "~" + NextAA + "~-~false" +  "&";
            result += rPepSeq + "~" + MinRMassDif +"~" + (MinRMassDif+ExpPepMass) + "~" + PreviousAA + "~" + NextAA + "~-~false" +  "&";
        }

        if(!bestResult.isEmpty())
        {
            result = bestResult + "$" + result;
        }

        return result;
    }

    public static List<String> CalculateTheoreticalPrecursors(IScan PsmScan, int MinZ, int MaxZ, double TheoPepMass, double ProtonMass, double StartMz, double EndMz)
    {
        List<String> PossiblePrecursorLi = new ArrayList<String>();
        for(int i = MinZ; i <= MaxZ ; i++) //Set up charge states ranging from 1 to 6
        {
            double CalMZ = (TheoPepMass + i * ProtonMass)/i;
            if((CalMZ >= StartMz) && (CalMZ <= EndMz))
            {
                PossiblePrecursorLi.add(CalMZ +"_" + i);
            }
        }

        return PossiblePrecursorLi;
    }

    public static ds_Peak MatchIsotopes(ISpectrum PrecursorSpectrum, double DownPrecursorWindow, double UpPrecursorWindow, int IsoNum, double Mz, int Z, double MassTol) throws FileParsingException
    {
        List<ds_Peak> PreCanLi = GetPrecursorPeaks(PrecursorSpectrum, DownPrecursorWindow, UpPrecursorWindow);
        Collections.sort(PreCanLi);

        List<Double> TheoreticalIsotopeLi = p_Operator.CalculateTheoreticalIsotopeMz(IsoNum, Mz, Z);
        //List<ds_Peak> detectedIsoLi = new ArrayList<ds_Peak>(); //Search if these isotope peaks are in the precursor windows
        ds_Peak tp = null;
        double minMz = Mz * (1 - (MassTol));
        double maxMz = Mz * (1 + (MassTol));
        for (ds_Peak p : PreCanLi)
        {
            if((p.mz >= minMz) && (p.mz <= maxMz))
            {
                if(tp != null)
                {
                    if(tp.intensity < p.intensity)
                    {
                        tp = p;
                    }
                }
                else
                {
                    tp = p;
                }
            }
        }
        return tp;
    }

    public static List<ds_Peak> GetPrecursorPeaks(ISpectrum PrecursorSpectrum, double DownPrecursorWindow, double UpPrecursorWindow)
    {
        double[] mzAry = PrecursorSpectrum.getMZs();
        double[] intAry = PrecursorSpectrum.getIntensities();

        double noiselevel = 0;
        for(double intensity : intAry)
        {
            noiselevel+=intensity;
        }
        noiselevel = noiselevel/intAry.length;

        List<ds_Peak> pLi = new ArrayList<ds_Peak>();
        for(int i = 0 ; i < mzAry.length ; i++) //Transform data structure
        {
            ds_Peak p = new ds_Peak();
            p.mz = mzAry[i];
            p.intensity = intAry[i];
            p.noiselevel = noiselevel;
            pLi.add(p);
        }

        List<ds_Peak> PreCanLi = new ArrayList<ds_Peak>();
        for(ds_Peak p : pLi)
        {
            if((p.mz >= DownPrecursorWindow)&&(p.mz <= UpPrecursorWindow))
            {
                PreCanLi.add(p);
            }
        }

        return PreCanLi;
    }

    public static String CheckNonSpecificPeptides(ds_Parameters inp, ds_PSM psm)
    {
        String bestResult = "";
        String result = "";

        String ProtSeq = "";
        if(inp.fastaMap.get(psm.ProteinAceNo)!=null){
            ProtSeq = inp.fastaMap.get(psm.ProteinAceNo);
        }
        else{
            System.out.println("Protein: "+psm.ProteinAceNo+" not found in the protein sequence database. Please check that you are using the same database as the search engine.");
            System.exit(1);
        }
        if(ProtSeq.contains("*")){
            System.out.println("Warning: The protein sequences in your fasta file include '*'. Please remove it before running Crystal-C.");
            ProtSeq=ProtSeq.replace("*","");
            //System.exit(1);
        }

        int Index = ProtSeq.indexOf(psm.PepSeq);
        if(Index < 0) //check I/L
        {	//Replace I/L with X
            String mProtSeq = ProtSeq.replaceAll("I", "X");
            mProtSeq = mProtSeq.replaceAll("L", "X");
            String mTargetPepSeq = psm.PepSeq.replaceAll("I", "X");
            mTargetPepSeq =  mTargetPepSeq.replaceAll("L",  "X");
            Index = mProtSeq.indexOf(mTargetPepSeq);
        }

        if (Index < 0)
        {
            System.out.println("Error. Can't find target peptides: "+psm.PepSeq+"@"+psm.ProteinAceNo+"   in the protein sequence.");
            System.exit(1);
        }

        double MinLMassDif = 1000000;
        String lPepSeq = "";
        String lAddSeq = "";
        for(int i=Index ; i>=0 ; i--) //check from left side of sequence
        {
            double lMassDif = psm.TheoPepMass + p_Operator.CalculatePepMonoMass(ProtSeq.substring(i, Index), inp.aa.AAMonoMassMap) - psm.ObvPepMass;
            if(lMassDif>psm.MassAbsThreshold)
            {
                break;
            }
            else if(Math.abs(lMassDif) <= psm.MassAbsThreshold)
            {
                lPepSeq = ProtSeq.substring(i, Index)+psm.PepSeq;
                lAddSeq = ProtSeq.substring(i, Index);
                MinLMassDif = lMassDif;
                break;
            }
        }

        double MinRMassDif = 1000000;
        String rPepSeq = "";
        String rAddSeq = "";
        int startp = Index+psm.PepSeq.length();
        for(int i = startp ; i<ProtSeq.length() ; i++)
        {
            double rMassDif = psm.TheoPepMass + p_Operator.CalculatePepMonoMass(ProtSeq.substring(startp, i), inp.aa.AAMonoMassMap) - psm.ObvPepMass;
            if(rMassDif>psm.MassAbsThreshold)
            {
                break;
            }
            else if(Math.abs(rMassDif) <= psm.MassAbsThreshold)
            {
                rPepSeq = psm.PepSeq+ProtSeq.substring(startp, i);
                rAddSeq = ProtSeq.substring(startp, i);
                MinRMassDif = rMassDif;
                break;
            }
        }

        //Find the peptide sequence with the minimum delta mass
        if ((Math.abs(MinLMassDif) <= psm.MassAbsThreshold) && (Math.abs(MinLMassDif) <= MinRMassDif))
        {
            bestResult = lPepSeq + "@" + MinLMassDif + "@"  + (psm.TheoPepMass+p_Operator.CalculatePepMonoMass(lAddSeq, inp.aa.AAMonoMassMap));
        }
        else if (Math.abs(MinRMassDif) <= psm.MassAbsThreshold)
        {
            //bestResult = rPepSeq + "@" + MinRMassDif + "@"  + p_Operator.CalculatePepMonoMass(rPepSeq, inp.aa.AAMonoMassMap);
            bestResult = rPepSeq + "@" + MinRMassDif + "@"  + (psm.TheoPepMass+p_Operator.CalculatePepMonoMass(rAddSeq, inp.aa.AAMonoMassMap));
        }

        //List all possible peptides
        if(!lPepSeq.isEmpty())
        {
            //Find the previous and next amino acid
            String PreviousAA = ProtSeq.indexOf(lPepSeq)>0?String.valueOf(ProtSeq.charAt(ProtSeq.indexOf(lPepSeq)-1)):"-";
            String NextAA = "-";

            result += lPepSeq + "~" + MinLMassDif + "~" + (psm.TheoPepMass+p_Operator.CalculatePepMonoMass(lAddSeq, inp.aa.AAMonoMassMap))
                    + "~" + PreviousAA + "~" + NextAA + "~-~false" + "&";
        }
        if(!rPepSeq.isEmpty())
        {
            String PreviousAA = "-";
            String NextAA = (ProtSeq.indexOf(rPepSeq)+rPepSeq.length()<ProtSeq.length())?String.valueOf(ProtSeq.charAt(ProtSeq.indexOf(rPepSeq)+rPepSeq.length())):"-";

            result += rPepSeq + "~" + MinRMassDif +"~" + (psm.TheoPepMass+p_Operator.CalculatePepMonoMass(rAddSeq, inp.aa.AAMonoMassMap))
                    + "~" + PreviousAA + "~" + NextAA + "~-~false" +  "&";
        }

        if(!bestResult.isEmpty())
        {
            result = bestResult + "$" + result;
        }

        return result;
    }
}

