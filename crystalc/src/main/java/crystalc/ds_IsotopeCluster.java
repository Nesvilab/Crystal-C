package crystalc;

import java.util.ArrayList;
import java.util.List;

public class ds_IsotopeCluster {
    public double PredictMonoMass;
    public double PredictMz;
    public double RepMz = 0;
    public ds_XIC PrecursorXic = new ds_XIC();
    public List<ds_XIC> IsotopeXicLi = new ArrayList<ds_XIC>();
}
