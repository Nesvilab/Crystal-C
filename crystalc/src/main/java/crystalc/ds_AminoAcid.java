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

import java.util.HashMap;
import java.util.Map;

public class ds_AminoAcid
{
    public Map<String, Double> AAMonoMassMap = new HashMap<String, Double>(); //Key: Amino Acid; Value: Monoisotopic Mass
    public Map<String, Double> AAAvgMassMap = new HashMap<String, Double>(); //Key: Amino Acid; Value: Average Mass
    public double ProtonMass = 1.00794;

    public void init()
    {
        //Monoisotopic Mass
        AAMonoMassMap.put("A", 71.03711d);
        AAMonoMassMap.put("B", 0d);
        AAMonoMassMap.put("C", 103.00919d);
        AAMonoMassMap.put("D", 115.02694d);
        AAMonoMassMap.put("E", 129.04259d);
        AAMonoMassMap.put("F", 147.06841d);
        AAMonoMassMap.put("G", 57.02146d);
        AAMonoMassMap.put("H", 137.05891d);
        AAMonoMassMap.put("I", 113.08406d);
        AAMonoMassMap.put("J", 0d);
        AAMonoMassMap.put("K", 128.09496d);
        AAMonoMassMap.put("L", 113.08406d);
        AAMonoMassMap.put("M", 131.04049d);
        AAMonoMassMap.put("N", 114.04293d);
        AAMonoMassMap.put("O", 114.07931d);
        AAMonoMassMap.put("P", 97.05276d);
        AAMonoMassMap.put("Q", 128.05858d);
        AAMonoMassMap.put("R", 156.10111d);
        AAMonoMassMap.put("S", 87.03203d);
        AAMonoMassMap.put("T", 101.04768d);
        AAMonoMassMap.put("U", 0d);
        AAMonoMassMap.put("V", 99.06841d);
        AAMonoMassMap.put("W", 186.07931d);
        AAMonoMassMap.put("X", 0d);
        AAMonoMassMap.put("Y", 163.06333d);
        AAMonoMassMap.put("Z", 0d);
        AAMonoMassMap.put("[", 10000.00000d);

        //Average Mass
        AAAvgMassMap.put("A", 71.0788d);
        AAAvgMassMap.put("B", 114.5962d);
        AAAvgMassMap.put("C", 103.1388d);
        AAAvgMassMap.put("D", 115.0886d);
        AAAvgMassMap.put("E", 129.1155d);
        AAAvgMassMap.put("F", 147.1766d);
        AAAvgMassMap.put("G", 57.0519d);
        AAAvgMassMap.put("H", 137.1411d);
        AAAvgMassMap.put("I", 113.1594d);
        AAAvgMassMap.put("J", 0.0000d);
        AAAvgMassMap.put("K", 128.1741d);
        AAAvgMassMap.put("L", 113.1594d);
        AAAvgMassMap.put("M", 131.1926d);
        AAAvgMassMap.put("N", 114.1038d);
        AAAvgMassMap.put("O", 114.1472d);
        AAAvgMassMap.put("P", 97.1167d);
        AAAvgMassMap.put("Q", 128.1307d);
        AAAvgMassMap.put("R", 156.1875d);
        AAAvgMassMap.put("S", 87.0782d);
        AAAvgMassMap.put("T", 101.1051d);
        AAAvgMassMap.put("U", 0.0000d);
        AAAvgMassMap.put("V", 99.1326d);
        AAAvgMassMap.put("W", 186.2132d);
        AAAvgMassMap.put("X", 113.1594d);
        AAAvgMassMap.put("Y", 163.1760d);
        AAAvgMassMap.put("Z", 128.6231d);
        AAAvgMassMap.put("[", 10000.00000d);
    }
}