package crystalc;

import umich.ms.fileio.filetypes.pepxml.jaxb.standard.SpectrumQuery;

public class ds_PSM
{
    public String Str;
    public int PsmScanNum;
    public String PepSeq;
    public int ObvZ;
    public double ObvMz;
    public double DeltaMass;
    public double ObvPepMass;
    public double TheoPepMass;
    public String ProteinAceNo;
    public String AlternativeProteins;
    public double MassAbsThreshold;
    public SpectrumQuery spectrum;

    public String NewPepSeq;
    public double NewTheoPepMass;
    public int NewZ;
    public double NewObvMz;
    public double NewObvPepMass;
    public double NewObvAbund;
    public String Note;
    public String AllPossiblePepSeq;
    public boolean UsePredictedMz;

    public ds_IsotopeCluster IsotopeCluster = new ds_IsotopeCluster();
}
