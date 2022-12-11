package crystalc;

import umich.ms.datatypes.scancollection.IScanCollection;
import umich.ms.datatypes.scan.IScan;
import umich.ms.datatypes.spectrum.ISpectrum;
import java.util.*;

public class p_PredictMonoisotope
{
    private int _PrecursorScanNum;
    private double _PrecursorMz;
    private int _Z;
    private double _MassTol;
    private double _MzTol;
    private double _RtTol;
    private IScanCollection _Scans;

    private List<Double> _TheoIsoLi = new ArrayList<Double>();
    private Double[] _RtAry;
    private double _MinMzRange;
    private double _MaxMzRange;
    private double _MinRtRange;
    private double _MaxRtRange;
    private TreeMap<Integer, ds_Peak[] > PeakMap = new TreeMap<Integer,ds_Peak[] >();

    public ds_IsotopeCluster IsotopeCluster = new ds_IsotopeCluster();

    public p_PredictMonoisotope(int PrecursorScanNum, double PrecursorMz, int Z, double MassTol, double MzTol, double RtTol, IScanCollection Scans, Double[] RtAry)
    {
        this._PrecursorScanNum = PrecursorScanNum;
        this._PrecursorMz = PrecursorMz;
        this._Z = Z;
        this._MassTol = MassTol;
        this._MzTol = MzTol;
        this._RtTol = RtTol;
        this._Scans = Scans;
        this._RtAry = RtAry;

        this._MinMzRange = this._PrecursorMz - this._MzTol;
        this._MaxMzRange = this._PrecursorMz + this._MzTol;
        this._MinRtRange = ((this._RtAry[this._PrecursorScanNum-1] - this._RtTol) > this._Scans.getRtMin()) ? (this._RtAry[this._PrecursorScanNum] - this._RtTol) : this._Scans.getRtMin();
        this._MaxRtRange = ((this._RtAry[this._PrecursorScanNum-1] + this._RtTol) <= this._Scans.getRtMax())? (this._RtAry[this._PrecursorScanNum] + this._RtTol) : this._Scans.getRtMax();
    }

    public void GetRepresentative()
    {
        //Generate theoretical isotopic peaks
        GenTheoIsoPeak();

        //Extract isotopic peaks from MS1 spectra
        ExtractPossibleIsotopes();

        //Project intensities and select the most intensive peak as the representative
        SelectRepresentative();
    }

    private void GenTheoIsoPeak()
    {
        int n = 4;
        for(int i=-n; i<=n; i++)
        {
            if(i!=0){
                double tMz = this._PrecursorMz + (1/(double)(this._Z))*i;
                this._TheoIsoLi.add(tMz);
            }
            else{
                this._TheoIsoLi.add(this._PrecursorMz);
            }
        }
        Collections.sort(this._TheoIsoLi);
    }

    private void ExtractPossibleIsotopes()
    {
        int lIndex = this._PrecursorScanNum - 1;
        int rIndex = this._PrecursorScanNum;

        int MaxCount = 3; //if the isotopic cluster is detected in far way from the max count.
        int lcount = 0;
        while((lIndex-1)>=0)//Search MS1 spectra before the precursor spectrum
        {
            if(this._RtAry[lIndex-1]!=null)
            {
                if(this._RtAry[lIndex-1] >= this._MinRtRange)
                {
                    if((this._Scans.getScanByNum(lIndex) != null) && (this._Scans.getScanByNum(lIndex).getMsLevel() == 1))
                    {
                        ds_Peak[] pAry = SearchIsotopes(lIndex);
                        pAry = CheckQuality(pAry);
                        if(pAry!=null){
                            PeakMap.put(lIndex, pAry);
                        }
                        else{
                            lcount +=1;
                        }
                        if(lcount>MaxCount){
                            break;
                        }
                    }
                }
                else
                {
                    break;
                }
            }

            lIndex -= 1;
        }

        int rcount = 0;
        while((rIndex-1) < this._RtAry.length) //Search MS1 spectra after the precursor spectrum
        {
            if(this._RtAry[rIndex-1]!=null){
                if(this._RtAry[rIndex-1] <= this._MaxRtRange)
                {
                    if((this._Scans.getScanByNum(rIndex) != null) && (this._Scans.getScanByNum(rIndex).getMsLevel() == 1))
                    {
                        ds_Peak[] pAry = SearchIsotopes(rIndex);
                        pAry = CheckQuality(pAry);
                        if(pAry!=null){
                            PeakMap.put(rIndex, pAry);
                        }
                        else{
                            rcount +=1;
                        }

                        if(rcount>MaxCount){
                            break;
                        }
                    }
                }
                else
                {
                    break;
                }
            }

            rIndex += 1;
        }
    }

    private ds_Peak[] CheckQuality(ds_Peak[] pAry)
    {
        if(pAry[4]!=null){ //check continuity
            int lindex = -1;
            for(int i=4; i>=0; i--){
                if(pAry[i]==null){
                    lindex=i+1;
                    break;
                }
            }
            int rindex = -1;
            for(int i=4; i<pAry.length; i++){
                if(pAry[i]==null){
                    rindex=i-1;
                    break;
                }
            }
            lindex = lindex<0? lindex=0: lindex;
            rindex = rindex<0? rindex=pAry.length: rindex;
            int count = rindex-lindex+1;
            if(count>1){
                for(int i=0;i<pAry.length;i++){
                    if((i<lindex)||(i>rindex)){
                        pAry[i]=null;
                    }
                }
            }
            else{
                pAry=null;
            }
        }
        else{
            pAry = null;
        }

        return pAry;
    }

    private ds_Peak[] SearchIsotopes(int index)
    {
        ds_Peak[] pAry = new ds_Peak[this._TheoIsoLi.size()];
        try
        {
            IScan scan = this._Scans.getScanByNum(index);
            ISpectrum spectrum = scan.fetchSpectrum();
            double[] MzAry = spectrum.getMZs();
            double[] IntAry = spectrum.getIntensities();
            for(int i=0; i<MzAry.length; i++)
            {
                if((MzAry[i] >= this._MinMzRange) && (MzAry[i] <= this._MaxMzRange))
                {
                    for(int j=0; j<this._TheoIsoLi.size(); j++)
                    {
                        double LowMzTol = this._TheoIsoLi.get(j) - this._TheoIsoLi.get(j) * this._MassTol;
                        double UpMzTol = this._TheoIsoLi.get(j) + this._TheoIsoLi.get(j) * this._MassTol;
                        if ((MzAry[i]>=LowMzTol) && (MzAry[i]<=UpMzTol))
                        {
                            if(pAry[j] != null){
                                if(pAry[j].intensity<IntAry[i]){
                                    ds_Peak p = new ds_Peak();
                                    p.mz = MzAry[i];
                                    p.intensity = IntAry[i];
                                    p.rt = scan.getRt() * 60;
                                    pAry[j] = p;
                                }
                            }
                            else{
                                ds_Peak p = new ds_Peak();
                                p.mz = MzAry[i];
                                p.intensity = IntAry[i];
                                p.rt = scan.getRt() * 60;
                                pAry[j] = p;
                            }
                            break;
                        }
                    }
                }
            }
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }

        return pAry;
    }

    private void SelectRepresentative()
    {
        ds_Peak[][] p2dAry = new ds_Peak[PeakMap.size()][9];
        int x=0;
        for(ds_Peak[] pAry : PeakMap.values()){
            for(int j=0; j<pAry.length;j++){
                if(pAry[j]!=null){
                    p2dAry[x][j]=pAry[j];
                }
            }
            x+=1;
        }

        //Select Representative
        TreeMap<Integer, Integer> tMap = new TreeMap<Integer, Integer>();
        for(int i=0;i<PeakMap.size();i++)
        {
            double MaxInt = 0;
            int index = -1;
            for(int j=0;j<9; j++)
            {
                if(p2dAry[i][j]!=null){
                    if (p2dAry[i][j].intensity>=MaxInt){
                        index = j;
                        MaxInt = p2dAry[i][j].intensity;
                    }
                }
            }
            if(tMap.containsKey(index))
            {
                int count = tMap.get(index);
                count += 1;
                tMap.remove(index);
                tMap.put(index, count);
            }
            else
            {
                tMap.put(index, 1);
            }
        }

        int MaxCount = 0;
        int rPindex = -1;
        for(int index : tMap.keySet()){
            int count = tMap.get(index);
            if(count>MaxCount){
                MaxCount  = count;
                rPindex = index;
            }
        }

        //Construct XIC
        double maxInt = 0;
        for(int j=0;j<9; j++)
        {
            ds_XIC xic = new ds_XIC();
            double sumInt = 0;
            for(int i=0;i<PeakMap.size();i++)
            {
                if(p2dAry[i][j]!=null){
                    if(xic.StartRt>p2dAry[i][j].rt)
                    {
                        xic.StartRt = p2dAry[i][j].rt;
                    }
                    if(xic.EndRt<p2dAry[i][j].rt)
                    {
                        xic.EndRt=p2dAry[i][j].rt;
                    }
                    if(p2dAry[i][j].intensity>xic.Height)
                    {
                        xic.Mz = p2dAry[i][j].mz;
                        xic.Rt = p2dAry[i][j].rt;
                        xic.Height = p2dAry[i][j].intensity;
                    }
                    sumInt += p2dAry[i][j].intensity;
                }
            }
            if(sumInt>0)
            {
                xic.Area = sumInt;
                IsotopeCluster.IsotopeXicLi.add(xic);
            }
            if(j==rPindex){
                IsotopeCluster.RepMz = xic.Mz;
            }
            if(j==4){
                IsotopeCluster.PrecursorXic = xic;
            }
        }
    }
}
