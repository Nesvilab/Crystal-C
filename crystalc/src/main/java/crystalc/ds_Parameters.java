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
