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

import java.io.*;
import java.util.*;

import umich.ms.datatypes.LCMSDataSubset;
import umich.ms.datatypes.scancollection.IScanCollection;
import umich.ms.fileio.exceptions.FileParsingException;
import umich.ms.fileio.filetypes.mzml.MZMLFile;
import umich.ms.fileio.filetypes.mzxml.MZXMLFile;
import umich.ms.datatypes.scancollection.impl.ScanCollectionDefault;
import umich.ms.fileio.filetypes.AbstractLCMSDataSource;
import umich.ms.fileio.filetypes.thermo.ThermoRawFile;

public class p_ReadData
{
    public static ds_Parameters ParseParameters(File ParamFile)
    {
        ds_Parameters inp = new ds_Parameters();
        try
        {
            BufferedReader br = new BufferedReader(new FileReader(ParamFile.getAbsolutePath()));
            String line = "";
            while ((line = br.readLine()) != null)
            {
                if (!line.equalsIgnoreCase("") && line.length()>0 && !line.startsWith("#"))
                {
                    String type = line.split("=")[0].trim();
                    String value = line.split("=")[1].trim();
                    value = value.contains("#")? value.substring(0, value.indexOf("#")).trim() : value.trim();
                    switch (type)
                    {
                        case "thread":
                        {
                            inp.NoThread =  (Integer.parseInt(value)<=0) ? (Runtime.getRuntime().availableProcessors() - 1):Integer.parseInt(value);
                            break;
                        }
                        case "fasta":
                        {
                            inp.FastaF = new File(value);
                            break;
                        }
                        case "raw_file_location":
                        {
                            inp.RawDataDictionary = value;
                            break;
                        }
                        case "raw_file_extension":
                        {
                            inp.RawFileExtension = value;
                            break;
                        }
                        case "output_location":
                        {
                            inp.OutputFolder = value;
                            break;
                        }
                        case "precursor_charge":
                        {
                            String[] strAry = value.split(" ");
                            inp.MinZ = Integer.parseInt(strAry[0]);
                            inp.MaxZ = Integer.parseInt(strAry[1]);
                            break;
                        }
                        case "isotope_number":
                        {
                            inp.IsoNum = Integer.parseInt(value);
                            break;
                        }
                        case "precursor_mass":
                        {
                            inp.MassTol = Double.parseDouble(value)/1000000;
                            break;
                        }
                        case "precursor_isolation_window":
                        {
                            inp.PrecursorIsolationWindow = Double.parseDouble(value);
                            break;
                        }
                        case "correct_isotope_error":
                        {
                            inp.CorrectIsotopeError = Boolean.parseBoolean(value);
                            break;
                        }
                    }
                }
            }
            br.close();
        }
        catch (Exception e)
        {
            e.printStackTrace();
            System.exit(1);
        }

        return inp;
    }

    public static List<Map> LoadFasta(File fastaF)
    {
        List<Map> MapLi = new ArrayList<Map>();
        Map<String, String> fastaMap = new HashMap<String, String>(); //Key: z; Value: ds_Peak
        Map<String, String> ProtIdMap = new HashMap<String, String>();
        try
        {
            BufferedReader br = new BufferedReader(new FileReader(fastaF));
            String line = br.readLine();
            String[] KeyAry = line.split(" |\t");
            String KeyStr = KeyAry[0].replace(">", "");
            line = line.replace(">", "");
            ProtIdMap.put(KeyStr, line);

            String ValueStr = "";
            while ((line = br.readLine()) != null)
            {
                if(line.contains(">"))
                {
                    KeyStr = KeyStr.replace(">", "");
                    fastaMap.put(KeyStr, ValueStr);

                    line = line.replace(">", "");
                    String tmp = line.split(" |\t")[0];
                    ProtIdMap.put(tmp, line);

                    KeyAry = line.split(" |\t");
                    KeyStr = KeyAry[0];
                    ValueStr = "";
                }
                else
                {
                    ValueStr += line;
                }
            }
            br.close();
            KeyStr = KeyStr.replace(">", "");
            fastaMap.put(KeyStr, ValueStr);
        }
        catch (Exception e)
        {
            e.printStackTrace();
            System.exit(1);
        }

        MapLi.add(fastaMap);
        MapLi.add(ProtIdMap);

        return MapLi;
    }

    public static IScanCollection LoadRawFile(String RawDataDictionary, String fName, String FileExtension) throws FileParsingException
    {
        AbstractLCMSDataSource<?> source = null;
        if(fName.contains("."))
        {
            fName=fName.substring(0,fName.indexOf("."));
        }
        if("mzXML".equalsIgnoreCase(FileExtension))
        {
            source = new MZXMLFile(RawDataDictionary + File.separator  + fName + "." + FileExtension);
        }
        else if("mzML".equalsIgnoreCase(FileExtension))
        {
            source = new MZMLFile(RawDataDictionary + File.separator + fName + "." + FileExtension);
        }
        else if("raw".equalsIgnoreCase(FileExtension))
        {
            source = new ThermoRawFile(RawDataDictionary + File.separator + fName + "." + FileExtension);
        }
        else
        {
            return null;
        }
        source.setExcludeEmptyScans(false);
        source.setNumThreadsForParsing(null);
        ScanCollectionDefault scans = new ScanCollectionDefault();
        scans.setDataSource(source);
        scans.loadData(LCMSDataSubset.WHOLE_RUN);

        return scans;
    }
}

