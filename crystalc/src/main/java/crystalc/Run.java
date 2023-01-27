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
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import javax.xml.stream.XMLStreamException;
import umich.ms.datatypes.scan.IScan;
import umich.ms.datatypes.scancollection.IScanCollection;
import umich.ms.datatypes.spectrum.ISpectrum;
import umich.ms.fileio.exceptions.FileParsingException;
import umich.ms.fileio.filetypes.pepxml.PepXmlParser;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.*;
import umich.ms.ims.isodistr.impl.AveragineTemplateProvider;

public class Run
{
    public static void main(String[] args) throws FileParsingException, XMLStreamException, IOException, InterruptedException, ExecutionException
    {
        Locale.setDefault(Locale.ROOT);

        File ParamFile = new File(args[0]);

        List<File> FileLi = new ArrayList<File>();
        for (int i = 1 ; i < args.length ; i++){
            FileLi.add(new File(args[i]));
        }

        ds_Parameters inp = p_ReadData.ParseParameters(ParamFile); //Parse parameter file
        inp.aa.init();
        List<Map> MapLi  = p_ReadData.LoadFasta(inp.FastaF); //Load Fasta file
        inp.fastaMap = MapLi.get(0);
        inp.ProtIdMap = MapLi.get(1);
        AveragineTemplateProvider templates = new AveragineTemplateProvider(200, 10000);

        for(File InputFile : FileLi){
            String FileExtension = InputFile.getName().substring(InputFile.getName().lastIndexOf('.'));
            if(FileExtension.equalsIgnoreCase(".pepXML") && !InputFile.getName().endsWith("_c.pepXML")){
                ProcessPepXML(inp, InputFile, templates);
            }
        }
    }

    private static void ProcessPepXML(ds_Parameters inp, File InputFile, AveragineTemplateProvider templates) throws FileParsingException, InterruptedException, ExecutionException, IOException, XMLStreamException
    {
        long start = System.currentTimeMillis();

        ExecutorService executor = Executors.newFixedThreadPool(inp.NoThread);
        List<Future<ds_PSM>> rStrLi = new ArrayList<Future<ds_PSM>>();
        List<ds_PSM> outputLi = new ArrayList<ds_PSM>();

        Path path = Paths.get(InputFile.getAbsolutePath());
        MsmsPipelineAnalysis analysis = PepXmlParser.parse(path);
        List<MsmsRunSummary> runSummaries = analysis.getMsmsRunSummary();
        for (MsmsRunSummary runSummary : runSummaries)
        {
            IScanCollection scans = p_ReadData.LoadRawFile(inp.RawDataDictionary, path.getFileName().toString(), inp.RawFileExtension);
            Double[] RtAry = GetRT(scans);
            SampleEnzyme se = runSummary.getSampleEnzyme();
            List<Specificity> SpecificityLi = se.getSpecificity();
            inp.CleavageSites = SpecificityLi.get(0).getCut();
            inp.EnzymeName = se.getName();
            inp.CleavageInhibitors = SpecificityLi.get(0).getNoCut();

            List<SpectrumQuery> spectrumQueries = runSummary.getSpectrumQuery();
            for (SpectrumQuery sq : spectrumQueries){
                WorkerThread worker = new WorkerThread(inp, sq, scans, templates, RtAry);
                Future<ds_PSM> result = (Future<ds_PSM>) executor.submit(worker);
                rStrLi.add(result);
            }
        }

        for(Future<ds_PSM> result : rStrLi){
            outputLi.add(result.get());
        }
        executor.shutdown();

       p_WriteResults.WritePepXML(InputFile, outputLi, analysis, inp);
       //p_WriteResults.WriteTxt(InputFile, rStrLi);
       //InputFile.delete(); //Delete original pepXML

        long end = System.currentTimeMillis();
        NumberFormat formatter = new DecimalFormat("#0.00000");
        System.out.println(InputFile.getName() +" Done!   The execution time is " + formatter.format((end - start) / (1000d * 60)) + " min.");
    }

    private static Double[] GetRT(IScanCollection scans)
    {
        TreeMap<Integer, IScan> Num2ScanMap = scans.getMapNum2scan();
        Double[] RtAry = new Double[Num2ScanMap.size()+2];
        try
        {
            RtAry[0] = 0.0;
            for(int num : Num2ScanMap.keySet()){
                IScan scan = Num2ScanMap.get(num);
                RtAry[num] = scan.getRt();
            }
        }
        catch(Exception e)
        {
            System.out.println("Please check if there are MS1 spectra in raw or mzML files.");
            e.printStackTrace();
            System.exit(1);
        }

        return RtAry;
    }
}