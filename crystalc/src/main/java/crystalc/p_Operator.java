/*
 * Licensed to the Apache Software Foundation (ASF) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The ASF licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing,
 * software distributed under the License is distributed on an
 * "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
 * KIND, either express or implied.  See the License for the
 * specific language governing permissions and limitations
 * under the License.
 */

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

