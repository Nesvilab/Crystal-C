package crystalc;

import java.io.File;
import java.util.List;
import java.util.Map;

public class ds_Parameters
{
    //Input Information
    public int NoThread = (Runtime.getRuntime().availableProcessors() - 1);
    public File FastaF;
    public String RawDataDictionary;
    public String RawFileExtension;
    public String OutputFolder;

    //Process Parameters
    public int IsoNum;
    public int MaxZ;
    public int MinZ;
    public double PrecursorIsolationWindow;
    public double MassTol; //in ppm
    public boolean UsePredictedMz = false;
    public String CleavageSites = "";
    public String EnzymeName = "";
    public String CleavageInhibitors = "";
    public boolean CorrectIsotopeError = true;

    public Map<String, List<String>> FileMap;
    public Map<String, String> fastaMap;
    public Map<String, String> ProtIdMap;
    public ds_AminoAcid aa = new ds_AminoAcid();
}
