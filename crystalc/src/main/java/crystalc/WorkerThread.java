package crystalc;

import java.util.List;
import java.util.TreeMap;
import java.util.concurrent.Callable;

import umich.ms.datatypes.scan.IScan;
import umich.ms.datatypes.scancollection.IScanCollection;
import umich.ms.datatypes.spectrum.ISpectrum;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.AltProteinDataType;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.SearchHit;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.SpectrumQuery;
import umich.ms.ims.isodistr.impl.AveragineTemplateProvider;

public class WorkerThread implements Callable<ds_PSM>
{
    private ds_Parameters _inp;
    private ds_PSM _psm = new ds_PSM();
    private IScanCollection _scans;
    private AveragineTemplateProvider _templates;
    private Double[] _RtAry;

    public WorkerThread(ds_Parameters inp, SpectrumQuery Spectrum, IScanCollection scans, AveragineTemplateProvider templates, Double[] RtAry)
    {
        this._inp = inp;
        this._scans = scans;
        this._RtAry = RtAry;

        this._psm.spectrum = Spectrum;
        this._psm.PsmScanNum = Integer.parseInt(String.valueOf(Spectrum.getStartScan()));
        this._psm.ObvZ = Spectrum.getAssumedCharge();
        this._psm.ObvMz = p_Operator.CalculateMz(Spectrum.getPrecursorNeutralMass(), Spectrum.getAssumedCharge());
        //this._psm.ObvPepMass = (Spectrum.getUncalibratedPrecursorNeutralMass()!=null)?Spectrum.getUncalibratedPrecursorNeutralMass():Spectrum.getPrecursorNeutralMass();
        this._psm.ObvPepMass = Spectrum.getPrecursorNeutralMass();

        SearchHit searchHit = (!Spectrum.getSearchResult().isEmpty()) ? Spectrum.getSearchResult().get(0).getSearchHit().get(0) : null;
        this._psm.PepSeq = searchHit.getPeptide();
        this._psm.TheoPepMass = searchHit.getCalcNeutralPepMass();
        String[] strAry = searchHit.getProtein().split(" ");
        this._psm.ProteinAceNo = strAry[0];
        this._psm.DeltaMass = searchHit.getMassdiff();
        this._psm.UsePredictedMz = this._inp.UsePredictedMz;

        String AlproteinStr = "";
        List<AltProteinDataType> AlProteinLi = searchHit.getAlternativeProtein();
        for(AltProteinDataType AlProtein : AlProteinLi)
        {
            String[] sAry = AlProtein.getProtein().split(" ");
            AlproteinStr += sAry[0] + ",";
        }
        this._psm.AlternativeProteins = (AlproteinStr == "") ? "" : AlproteinStr.substring(0, AlproteinStr.lastIndexOf(','));
        this._psm.MassAbsThreshold = (this._psm.ObvPepMass * inp.MassTol);

        this._templates = templates;
    }

    @Override
    public ds_PSM call()
    {
        String resultStr = "";
        try
        {
            TreeMap<Integer, IScan> Num2ScanMap = this._scans.getMapNum2scan();
            double RtRange = 0.5;

            //Initialization
            this._psm.NewObvPepMass = this._psm.ObvPepMass;
            this._psm.NewObvMz = this._psm.ObvMz;
            this._psm.NewZ = this._psm.ObvZ;
            this._psm.NewPepSeq = this._psm.PepSeq;
            this._psm.NewTheoPepMass = this._psm.TheoPepMass;
            this._psm.NewObvAbund = this._psm.IsotopeCluster.PrecursorXic.Area;

            //Get precursor spectrum
            IScan PsmScan = this._scans.getScanByNum(this._psm.PsmScanNum);
            int PrecursorScanNum = PsmScan.getPrecursor().getParentScanNum()!=null?PsmScan.getPrecursor().getParentScanNum() : -1;

            //Get the precursor window in the parent spectrum
            double StartMz = (PsmScan.getPrecursor().getMzRangeStart() != null) ? PsmScan.getPrecursor().getMzRangeStart() : (this._psm.ObvMz - this._inp.PrecursorIsolationWindow);
            double EndMz =  (PsmScan.getPrecursor().getMzRangeEnd() != null) ? PsmScan.getPrecursor().getMzRangeEnd() : (this._psm.ObvMz + this._inp.PrecursorIsolationWindow);
            double TargetMz = (PsmScan.getPrecursor().getMzTargetMono()!=null) ? PsmScan.getPrecursor().getMzTargetMono(): PsmScan.getPrecursor().getMzTarget() ;

            double DeltaMass = -1;
            this._psm.IsotopeCluster = (PrecursorScanNum>=0)? p_ProcessPSM.CalculateMonoisotopeMass(PrecursorScanNum, TargetMz, this._inp, this._scans, this._templates, this._psm.ObvZ, 1.5, RtRange, this._RtAry) : null;   //Get the predicted monoisotope mass/Mz
            if(this._psm.UsePredictedMz && this._psm.IsotopeCluster!=null){
                DeltaMass = ((this._psm.ObvPepMass - this._psm.TheoPepMass) > (this._psm.IsotopeCluster.PredictMonoMass - this._psm.TheoPepMass)) ?
                        (this._psm.ObvPepMass - this._psm.TheoPepMass) :  (this._psm.IsotopeCluster.PredictMonoMass - this._psm.TheoPepMass) ;
            }
            else{
                DeltaMass = this._psm.ObvPepMass - this._psm.TheoPepMass;
            }

            if(Math.abs(DeltaMass) >= 57) //The minimum amino acid mass: 57.02146 (G)
            {
                if((DeltaMass > 0) && (!this._inp.EnzymeName.equalsIgnoreCase("nonspecific")))
                {
                    //region Check missed cleavage peptides
                    resultStr = p_ProcessPSM.CheckMissedCleavagePeptides(this._inp, this._psm, this._psm.ObvPepMass);

                    if(this._psm.UsePredictedMz) { //Using predicted Mz
                        resultStr = (resultStr != "") ?  (resultStr+"$(ExpMass)") : p_ProcessPSM.CheckMissedCleavagePeptides(this._inp, this._psm, this._psm.IsotopeCluster.PredictMonoMass);
                    }

                    if(resultStr != "")
                    {
                        String[] strAry = resultStr.split("[$\\@]");
                        this._psm.NewPepSeq = strAry[0];
                        this._psm.NewTheoPepMass = Double.parseDouble(strAry[2]);
                        this._psm.AllPossiblePepSeq = strAry[3];

                        if(this._psm.UsePredictedMz) { //Using predicted Mz
                            this._psm.NewObvMz = resultStr.contains("(ExpMass)") ? this._psm.NewObvMz : this._psm.IsotopeCluster.PredictMz;
                            resultStr = resultStr.contains("(ExpMass)") ? "Found Missed Cleavage Peptides (ExpMass)\t" : "Found Missed Cleavage Peptides (PredictMass)\t";
                        }
                        else{
                            resultStr = "Found Missed Cleavage Peptides (ExpMass)\t";
                        }

                        resultStr += this._psm.AllPossiblePepSeq;
                    }
                    //endregion
                }
                else if((DeltaMass > 0) && (this._inp.EnzymeName.equalsIgnoreCase("nonspecific")))
                {
                    //region Check nonspecific peptides
                    resultStr = p_ProcessPSM.CheckNonSpecificPeptides(this._inp, this._psm);
                    if(resultStr != "")
                    {
                        String[] strAry = resultStr.split("[$\\@]");
                        this._psm.NewPepSeq = strAry[0];
                        this._psm.NewTheoPepMass = Double.parseDouble(strAry[2]);
                        this._psm.AllPossiblePepSeq = strAry[3];

                        resultStr = "Found Nonspecific Peptides (ExpMass) \t";
                        resultStr += this._psm.AllPossiblePepSeq;
                    }
                    //endregion
                }
                else
                {
                    //region Check semi-enzymatic peptides
                    resultStr = p_ProcessPSM.CheckSemiTrypticPeptides(this._psm.PepSeq, this._psm.TheoPepMass, this._psm.ObvPepMass, this._psm.MassAbsThreshold, this._inp.aa.AAMonoMassMap);

                    if(this._psm.UsePredictedMz) { //Using predicted Mz
                        resultStr =  (resultStr != "") ?  (resultStr+"$(ExpMass)") : p_ProcessPSM.CheckSemiTrypticPeptides(this._psm.PepSeq, this._psm.TheoPepMass,
                                this._psm.IsotopeCluster.PredictMonoMass, this._psm.MassAbsThreshold, this._inp.aa.AAMonoMassMap);
                    }

                    if(resultStr != "")
                    {
                        String[] strAry = resultStr.split("[$\\@]");
                        this._psm.NewPepSeq = strAry[0];
                        this._psm.NewTheoPepMass = Double.parseDouble(strAry[2]);
                        this._psm.AllPossiblePepSeq = strAry[3];

                        if(this._psm.UsePredictedMz) { //Using predicted Mz
                            this._psm.NewObvMz = resultStr.contains("(ExpMass)") ? this._psm.NewObvMz : this._psm.IsotopeCluster.PredictMz;
                            resultStr = resultStr.contains("(ExpMass)") ? "Found Semi-enzymatic Peptides (ExpMass) \t" : "Found Semi-enzymatic Peptides (PredictMass) \t";
                        }
                        else{
                            resultStr = "Found Semi-enzymatic Peptides (ExpMass) \t";
                        }

                        resultStr += this._psm.AllPossiblePepSeq;
                    }
                    //endregion
                }
            }

            if ((resultStr == "") && (PrecursorScanNum>=0))
            {
                IScan PrecursorScan = this._scans.getScanByNum(PrecursorScanNum);

                //region Check co-fragmentation
                List<String> PossiblePrecursorLi = p_ProcessPSM. CalculateTheoreticalPrecursors(PsmScan, this._inp.MinZ, this._inp.MaxZ, this._psm.TheoPepMass, this._inp.aa.ProtonMass, StartMz, EndMz);
                if(PossiblePrecursorLi.size() == 1)
                {
                    String[] strAry = PossiblePrecursorLi.get(0).split("_"); //get m/z; z.
                    double TheoreticalMonoMz = Double.parseDouble(strAry[0]);
                    int TheoreticalZ = Integer.parseInt(strAry[1]);
                    List<Double> TheoreticalIsotopeLi = p_Operator.CalculateTheoreticalIsotopeMz(this._inp.IsoNum, TheoreticalMonoMz, TheoreticalZ);

                    //region 1. Check if PsmScan.getPrecursor().getMzTarget belong to the isotopic pattern of the identified peptide
                    boolean isFit = false;
                    for(double isoMz : TheoreticalIsotopeLi)
                    {
                        double minMz = isoMz * (1 - this._inp.MassTol);
                        double maxMz = isoMz * (1 + this._inp.MassTol);
                        if((TargetMz >= minMz) && (TargetMz <= maxMz))
                        {
                            isFit = true;
                            break;
                        }
                    }
                    //endregion

                    //region 2. Search target peaks
                    double downRT = (PrecursorScan.getRt() - RtRange)>0 ? (PrecursorScan.getRt() - RtRange) : 0;
                    ds_Peak tp = null;
                    int foundMS1ScanNum = -1;
                    for(int i = this._psm.PsmScanNum; i >= 1; i--)
                    {
                        IScan scan = Num2ScanMap.get(i);
                        if(scan!=null){
                            if(scan.getRt() >= downRT) //Search previous MS1 spectrum within the rtTol
                            {
                                if (scan.getMsLevel() == 1)
                                {
                                    ISpectrum Ms1Spectrum = scan.fetchSpectrum();
                                    tp = p_ProcessPSM.MatchIsotopes(Ms1Spectrum, StartMz, EndMz, this._inp.IsoNum, TheoreticalMonoMz, TheoreticalZ, this._inp.MassTol);
                                    if(tp!=null){
                                        foundMS1ScanNum = i;
                                        break;
                                    }
                                }
                            }
                            else
                            {
                                break;
                            }
                        }

                    }
                    //endregion

                    //region 3. Annotate
                    if(tp!=null)
                    {
                        this._psm.NewObvMz =  (isFit && this._psm.ObvZ==TheoreticalZ) ? this._psm.ObvMz : tp.mz;
                        this._psm.NewObvPepMass =  (isFit && this._psm.ObvZ==TheoreticalZ)? this._psm.ObvPepMass:
                                tp.mz*TheoreticalZ-TheoreticalZ*this._inp.aa.ProtonMass;
                        this._psm.NewZ = TheoreticalZ;

                        this._psm.NewObvAbund = this._psm.IsotopeCluster.PrecursorXic.Area;
                        this._psm.UsePredictedMz = false;

                        resultStr = isFit ? "Found isotopic peak - " : "Found  another isotopic peak - ";
                        resultStr +=  ((foundMS1ScanNum == PrecursorScanNum) ? "Precursor Spectrum" : "Previous MS1 Spectrum") +" - "
                                + ((TheoreticalZ == this._psm.ObvZ) ? "Same Z": "Different Z");
                    }
                    else
                    {
                        resultStr = "Can't find any possible precursors.";
                    }
                    //endregion
                }
                else
                {
                    resultStr = "No possible precursors can be found in the precursor window";
                }
                //endregion
            }

            this._psm.Note = PrecursorScanNum + "\t" + resultStr;
        }
        catch(Exception e)
        {
            System.out.println("Error at: " + this._psm.spectrum.getSpectrum());
            e.printStackTrace();
            System.exit(1);
        }
        return this._psm;
    }
}
