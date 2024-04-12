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
import java.nio.channels.FileChannel;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

import javax.xml.XMLConstants;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathFactory;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.*;

public class p_WriteResults
{
    static DecimalFormat df = new DecimalFormat("#.#####");
    static double ProtonMass = 1.00794; //need to think how to design the program

    public static void WriteTxt(File InputFile, List<Future<ds_PSM>> rStrLi) throws IOException, InterruptedException, ExecutionException
    {
        BufferedWriter wr = new BufferedWriter(new FileWriter(InputFile.getAbsolutePath().replace(".pepXML", "_c.txt")));
        wr.write("Spectrum \t Peptide \t Precursor Neutral Mass \t Theoretical Mass \t Delta Mass \t Charge \t Retention \tPredict MonoMass \t Predict Mz " +
                "\t NewObvMz \t NewObvPepMass \t NewTheoPepMass \t New Delta Mass \t New Z \t Obv Mz Abundance \t Num of XIC \t Xics \t Precursor ScanNum \t Note");
        wr.newLine();
        for(Future<ds_PSM> result : rStrLi)
        {
            ds_PSM psm = result.get();

            String XicStr = "";
            if(psm.IsotopeCluster!=null)
            {
                if(psm.IsotopeCluster.IsotopeXicLi.size()>0)
                {
                    Collections.sort(psm.IsotopeCluster.IsotopeXicLi);
                    for(ds_XIC xic : psm.IsotopeCluster.IsotopeXicLi)
                    {
                        XicStr += xic.Mz + "_" + xic.Rt + "_" + xic.Height +"_"+ xic.StartRt + "_" + xic.EndRt + ";";
                    }
                }
            }
            SearchHit searchHit = (!psm.spectrum.getSearchResult().isEmpty()) ? psm.spectrum.getSearchResult().get(0).getSearchHit().get(0) : null;

            double NewTheoPepM = psm.NewTheoPepMass;
            if (psm.AllPossiblePepSeq != null) {
                String[] PepStrAry = psm.AllPossiblePepSeq.split("&");
                String[] SubStrAry = PepStrAry[0].split("~");
                NewTheoPepM = psm.NewObvPepMass -Double.parseDouble(SubStrAry[1]);
            }

            wr.write(psm.spectrum.getSpectrum() + "\t" + searchHit.getPeptide() + "\t" + psm.spectrum.getPrecursorNeutralMass() + "\t"
                    + psm.TheoPepMass  + "\t" + (psm.spectrum.getPrecursorNeutralMass() - psm.TheoPepMass) + "\t"
                    + psm.spectrum.getAssumedCharge() + "\t" + psm.spectrum.getRetentionTimeSec() + "\t"
                    + (psm.IsotopeCluster!=null?psm.IsotopeCluster.PredictMonoMass:" ") + "\t" + (psm.IsotopeCluster!=null?psm.IsotopeCluster.PredictMz:" ") + "\t" + psm.NewObvMz + "\t"
                    + df.format(psm.NewObvPepMass) +"\t" + df.format(NewTheoPepM) +"\t"+ df.format(psm.NewObvPepMass - NewTheoPepM) + "\t" + psm.NewZ
                    + "\t" + psm.NewObvAbund + "\t" + (psm.IsotopeCluster!=null?psm.IsotopeCluster.IsotopeXicLi.size():" ") + "\t" + XicStr + "\t" + psm.Note);
            wr.newLine();
        }
        wr.close();
    }

    public static void WritePepXML(File InputFile, List<ds_PSM> outputLi, MsmsPipelineAnalysis analysis,  ds_Parameters inp) throws IOException, XMLStreamException, InterruptedException, ExecutionException
    {
        try{
            Path path = Paths.get(inp.OutputFolder +File.separator+InputFile.getName().replaceAll(".pepXML", "_c.pepXML"));
            OutputStream os = new BufferedOutputStream(Files.newOutputStream(path));

            XMLOutputFactory outputFactory = XMLOutputFactory.newDefaultFactory();
            XMLStreamWriter writer = null;
            writer = outputFactory.createXMLStreamWriter(os, "utf-8");
            writer.writeStartDocument("utf-8", "1.0");
            writer.writeCharacters(System.getProperty("line.separator"));

            writer.writeProcessingInstruction("xml-stylesheet", "type=\"text/xsl\" href=\"myXsl.xsl\"");
            writer.writeCharacters(System.getProperty("line.separator"));

            writer.writeStartElement("msms_pipeline_analysis");
            writer.writeAttribute("date", analysis.getDate().toString());
            writer.writeAttribute("xmlns", "http://regis-web.systemsbiology.net/pepXML");
            writer.writeAttribute("summary_xml", analysis.getSummaryXml());
            writer.writeAttribute("xsi:schemaLocation", "http://sashimi.sourceforge.net/schema_revision/pepXML/pepXML_v118.xsd");
            writer.writeAttribute("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");
            writer.writeCharacters(System.getProperty("line.separator"));

            List<MsmsRunSummary> MsmsRunSummaryLi = analysis.getMsmsRunSummary();
            for(MsmsRunSummary mrs : MsmsRunSummaryLi)
            {
                writer.writeStartElement("msms_run_summary");
                writer.writeAttribute("base_name", mrs.getBaseName());
                writer.writeAttribute("raw_data_type", mrs.getRawDataType());
                if(mrs.getComment()!=null){
                    writer.writeAttribute("comment", mrs.getComment());
                }
                writer.writeAttribute("raw_data", mrs.getRawData());
                writer.writeCharacters(System.getProperty("line.separator"));

                //Sample Enzyme
                SampleEnzyme se = mrs.getSampleEnzyme();
                writer.writeStartElement("sample_enzyme");
                writer.writeAttribute("name", se.getName());
                writer.writeCharacters(System.getProperty("line.separator"));
                List<Specificity> SpecificityLi = se.getSpecificity();
                for(Specificity sf : SpecificityLi)
                {
                    writer.writeEmptyElement("specificity");
                    writer.writeAttribute("cut", sf.getCut());
                    writer.writeAttribute("no_cut", sf.getNoCut());
                    writer.writeAttribute("sense", sf.getSense());
                    writer.writeCharacters(System.getProperty("line.separator"));
                }
                writer.writeEndElement(); //corresponding to sample_enzyme
                writer.writeCharacters(System.getProperty("line.separator"));

                //Search Summary
                List<SearchSummary> SearchSummaryLi = mrs.getSearchSummary();
                for(SearchSummary ss : SearchSummaryLi)
                {
                    writer.writeStartElement("search_summary");
                    writer.writeAttribute("base_name", ss.getBaseName());
                    writer.writeAttribute("precursor_mass_type", ss.getPrecursorMassType().value().toString());
                    writer.writeAttribute("search_engine", ss.getSearchEngine().value().toString());
                    writer.writeAttribute("search_engine_version", ss.getSearchEngineVersion());
                    writer.writeAttribute("fragment_mass_type", ss.getFragmentMassType().value().toString());
                    writer.writeAttribute("search_id", String.valueOf(ss.getSearchId()));
                    writer.writeCharacters(System.getProperty("line.separator"));

                    //Search Database
                    SearchDatabase sd = ss.getSearchDatabase();
                    writer.writeEmptyElement("search_database");
                    writer.writeAttribute("local_path", sd.getLocalPath());
                    writer.writeAttribute("type", sd.getType());
                    writer.writeCharacters(System.getProperty("line.separator"));

                    //Enzymatic Search Constraint
                    EnzymaticSearchConstraint esc = ss.getEnzymaticSearchConstraint();
                    writer.writeEmptyElement("enzymatic_search_constraint");
                    writer.writeAttribute("enzyme", esc.getEnzyme());
                    writer.writeAttribute("min_number_termini", esc.getMinNumberTermini().toString());
                    writer.writeAttribute("max_num_internal_cleavages", esc.getMaxNumInternalCleavages().toString());
                    writer.writeCharacters(System.getProperty("line.separator"));

                    //Amino Acid Modification
                    List<AminoacidModification> AminoacidModificationLi = ss.getAminoacidModification();
                    for(AminoacidModification aam : AminoacidModificationLi)
                    {
                        writer.writeEmptyElement("aminoacid_modification");
                        writer.writeAttribute("aminoacid", aam.getAminoacid());
                        writer.writeAttribute("massdiff", String.valueOf(aam.getMassdiff()));
                        if(aam.getProteinTerminus() != null) {
                            writer.writeAttribute("protein_terminus", aam.getProteinTerminus());
                        }
                        writer.writeAttribute("mass", String.valueOf(aam.getMass()));
                        if(aam.getVariable() != null){
                            writer.writeAttribute("variable", aam.getVariable());
                        }
                        writer.writeCharacters(System.getProperty("line.separator"));
                    }

                    //Terminal modification
                    List<TerminalModification> TerminalModificationLi = ss.getTerminalModification();
                    for(TerminalModification tm : TerminalModificationLi)
                    {
                        writer.writeEmptyElement("terminal_modification");
                        writer.writeAttribute("massdiff", String.valueOf(tm.getMassdiff()));
                        writer.writeAttribute("protein_terminus", String.valueOf(tm.getProteinTerminus()));
                        writer.writeAttribute("mass", String.valueOf(tm.getMass()));
                        writer.writeAttribute("terminus", String.valueOf(tm.getTerminus()));
                        writer.writeAttribute("variable", String.valueOf(tm.getVariable()));

                        writer.writeCharacters(System.getProperty("line.separator"));
                    }

                    //Parameter
                    List<NameValueType> ParamNameLi = ss.getParameter();
                    for(NameValueType nt : ParamNameLi)
                    {
                        writer.writeEmptyElement("parameter");
                        writer.writeAttribute("name", nt.getName());
                        writer.writeAttribute("value", nt.getValueStr());
                        writer.writeCharacters(System.getProperty("line.separator"));
                    }

                    writer.writeEndElement();
                    writer.writeCharacters(System.getProperty("line.separator"));
                }

                //Spectrum Query
                for(ds_PSM psm : outputLi)
                {
                    //double ObvMass = ((psm.spectrum.getUncalibratedPrecursorNeutralMass()!=null) && (psm.NewObvPepMass!=psm.spectrum.getUncalibratedPrecursorNeutralMass()))?psm.NewObvPepMass: psm.spectrum.getPrecursorNeutralMass();
                    double ObvMass = psm.NewObvPepMass;

                    writer.writeStartElement("spectrum_query");
                    writer.writeAttribute("start_scan", String.valueOf(psm.spectrum.getStartScan()));
                    if(psm.spectrum.getUncalibratedPrecursorNeutralMass()!=null){
                        writer.writeAttribute("uncalibrated_precursor_neutral_mass", df.format(psm.spectrum.getUncalibratedPrecursorNeutralMass()));
                    }
                    writer.writeAttribute("assumed_charge", String.valueOf(psm.NewZ));
                    if(psm.spectrum.getSpectrumNativeID()!=null){
                        writer.writeAttribute("spectrumNativeID", psm.spectrum.getSpectrumNativeID());
                    }
                    writer.writeAttribute("spectrum", psm.spectrum.getSpectrum());
                    writer.writeAttribute("end_scan", String.valueOf(psm.spectrum.getEndScan()));
                    writer.writeAttribute("index", String.valueOf(psm.spectrum.getIndex()));
                    writer.writeAttribute("precursor_neutral_mass",df.format(ObvMass));
                    writer.writeAttribute("retention_time_sec", String.valueOf(psm.spectrum.getRetentionTimeSec()));
                    writer.writeCharacters(System.getProperty("line.separator"));

                    //Search Result
                    writer.writeStartElement("search_result");
                    writer.writeCharacters(System.getProperty("line.separator"));

                    List<SearchResult> SearchResultLi = psm.spectrum.getSearchResult();
                    for(SearchResult sr : SearchResultLi)
                    {
                        List<SearchHit> SearchHitLi = sr.getSearchHit();
                        for(int i = 0;  i < SearchHitLi.size(); i++)
                        {
                            SearchHit sh = SearchHitLi.get(i);
                            if(i==0) //only modify top hit
                            {
                                if (psm.AllPossiblePepSeq != null)
                                {
                                    String[] PepStrAry = psm.AllPossiblePepSeq.split("&");
                                    String[] SubStrAry = PepStrAry[0].split("~");
                                    try
                                    {
                                        double deltaMass = ObvMass-Double.parseDouble(SubStrAry[2]);
                                        WriteSearchHit(writer, sh, sh.getPeptide(), SubStrAry[0], df.format(deltaMass), df.format(Double.valueOf(SubStrAry[2])),
                                                SubStrAry[3], SubStrAry[4], SubStrAry[5], SubStrAry[6]);
                                    }
                                    catch(Exception e)
                                    {
                                        e.printStackTrace();
                                        System.exit(1);
                                    }
                                }
                                else
                                {
                                    double deltaMass = ObvMass-psm.NewTheoPepMass;
                                    WriteSearchHit(writer, sh, sh.getPeptide(), psm.NewPepSeq, df.format(deltaMass), df.format(psm.NewTheoPepMass), "-", "-", "-", "False");
                                }
                            }
                            else
                            {
                                WriteSearchHit(writer, sh, sh.getPeptide(), sh.getPeptide(), df.format(sh.getMassdiff()), df.format(sh.getCalcNeutralPepMass()),  "-", "-", "-",  "False");
                            }
                        }
                    }

                    writer.writeEndElement(); //corresponding to search_result
                    writer.writeCharacters(System.getProperty("line.separator"));

                    writer.writeEndElement(); //corresponding to spectrum_query
                    writer.writeCharacters(System.getProperty("line.separator"));
                }

                writer.writeEndElement();  //corresponding to msms_run_summary
                writer.writeCharacters(System.getProperty("line.separator"));
            }

            writer.writeEndElement();
            writer.writeEndDocument();
            writer.close();
        }
        catch (Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }

    private static void WriteSearchHit(XMLStreamWriter writer, SearchHit sh, String OldPepSeq, String NewPepSeq, String NewDeltaMass, String NewPepMass, String peptide_prev_aa, String peptide_next_aa, String ProtId, String AddNumMissedCleavage) throws XMLStreamException
    {
        try{
            writer.writeStartElement("search_hit");
            writer.writeAttribute("peptide", NewPepSeq);
            writer.writeAttribute("massdiff", NewDeltaMass);
            writer.writeAttribute("calc_neutral_pep_mass", NewPepMass);
            writer.writeAttribute("peptide_next_aa", (peptide_next_aa.equals("-")) ? sh.getPeptideNextAa() : peptide_next_aa);
            writer.writeAttribute("num_missed_cleavages", (!AddNumMissedCleavage.equals("True")) ? String.valueOf(sh.getNumMissedCleavages()) : String.valueOf(sh.getNumMissedCleavages() + 1));
            writer.writeAttribute("num_tol_term", (NewPepSeq.length()<sh.getPeptide().length())? "1" : String.valueOf(sh.getNumTolTerm()));
            if(sh.getProteinDescr()!=null){
                writer.writeAttribute("protein_descr", sh.getProteinDescr());
            }
            writer.writeAttribute("num_tot_proteins", String.valueOf(sh.getNumTotProteins()));
            writer.writeAttribute("tot_num_ions", String.valueOf(sh.getTotNumIons()));
            writer.writeAttribute("hit_rank", String.valueOf(sh.getHitRank()));
            writer.writeAttribute("num_matched_ions", String.valueOf(sh.getNumMatchedIons()));
            writer.writeAttribute("protein", (ProtId.equals("-")) ? sh.getProtein() : ProtId);
            writer.writeAttribute("peptide_prev_aa", (peptide_prev_aa.equals("-")) ? sh.getPeptidePrevAa() : peptide_prev_aa);
            writer.writeAttribute("is_rejected", String.valueOf(sh.getIsRejected()));
            writer.writeCharacters(System.getProperty("line.separator"));

            //alternative protein
            if(sh.getAlternativeProtein() != null)
            {
                List<AltProteinDataType> altpLi = sh.getAlternativeProtein();
                for(AltProteinDataType altp : altpLi)
                {
                    writer.writeEmptyElement("alternative_protein");
                    if(altp.getProteinDescr()!=null){
                        writer.writeAttribute("protein_descr", altp.getProteinDescr());
                    }
                    writer.writeAttribute("protein", altp.getProtein());
                    writer.writeAttribute("peptide_prev_aa", altp.getPeptidePrevAa());
                    writer.writeAttribute("peptide_next_aa", altp.getPeptideNextAa());
                    writer.writeAttribute("num_tol_term", altp.getNumTolTerm().toString());
                    writer.writeCharacters(System.getProperty("line.separator"));
                }
            }

            //modification_info
            if(sh.getModificationInfo() != null)
            {
                ModificationInfo mi = sh.getModificationInfo();
                writer.writeStartElement("modification_info");

                if(mi.getModNtermMass() != null)
                    writer.writeAttribute("mod_nterm_mass", String.valueOf(mi.getModNtermMass()));
                if(mi.getModCtermMass() != null)
                    writer.writeAttribute("mod_cterm_mass", String.valueOf(mi.getModCtermMass()));
                if(mi.getModifiedPeptide() != null)
                    writer.writeAttribute("modified_peptide", String.valueOf(mi.getModifiedPeptide().replace(OldPepSeq,NewPepSeq)));

                writer.writeCharacters(System.getProperty("line.separator"));

                List<ModAminoacidMass> ModAminoacidMassLi = mi.getModAminoacidMass();
                for(ModAminoacidMass maam : ModAminoacidMassLi)
                {
                    int index = CheckModPos(String.valueOf(sh.getPeptide()), NewPepSeq, maam.getPosition());

                    if  (index!=99999999) {
                        writer.writeEmptyElement("mod_aminoacid_mass");
                        writer.writeAttribute("mass", String.valueOf(maam.getMass()));
                        writer.writeAttribute("position", String.valueOf(index));
                        writer.writeCharacters(System.getProperty("line.separator"));
                    }
                }
                writer.writeEndElement();
                writer.writeCharacters(System.getProperty("line.separator"));
            }

            //search_score
            List<NameValueType> NameValueTypeLi = sh.getSearchScore();
            for(NameValueType nvt : NameValueTypeLi){
                writer.writeEmptyElement("search_score");
                writer.writeAttribute("name", nvt.getName());
                writer.writeAttribute("value", nvt.getValueStr());
                writer.writeCharacters(System.getProperty("line.separator"));
            }

            //PTM result
            if(sh.getPtmResult() != null)
            {
                List<PtmResult> ptmLi = sh.getPtmResult();
                for(PtmResult ptm : ptmLi)
                {
                    writer.writeEmptyElement("ptm_result");
                    if (ptm.getLocalization() != null) {
                        writer.writeAttribute("localization", ptm.getLocalization());
                    }
                    if (ptm.getBestPositions() != null) {
                        writer.writeAttribute("best_positions", ptm.getBestPositions());
                    }
                    if (ptm.getNumberBestPositions() != null) {
                        writer.writeAttribute("number_best_positions", ptm.getNumberBestPositions().toString());
                    }
                    if (ptm.getBestScoreWithPtm() != null) {
                        writer.writeAttribute("best_score_with_ptm", ptm.getBestScoreWithPtm().toString());
                    }
                    if (ptm.getScoreWithoutPtm() != null) {
                        writer.writeAttribute("score_without_ptm", ptm.getScoreWithoutPtm().toString());
                    }
                    if (ptm.getLocalizationPeptide() != null) {
                        writer.writeAttribute("localization_peptide", ptm.getLocalizationPeptide());
                    }
                    if (ptm.getSecondBestScoreWithPtm() != null) {
                        writer.writeAttribute("second_best_score_with_ptm", ptm.getSecondBestScoreWithPtm().toString());
                    }
                    if (ptm.getPtmMass() != null) {
                        writer.writeAttribute("ptm_mass", ptm.getPtmMass().toString());
                    }
                    if (ptm.getBestPosition() != null) {
                        writer.writeAttribute("best_position", ptm.getBestPosition().toString());
                    }
                    if (ptm.getScoreAllUnshifted() != null) {
                        writer.writeAttribute("score_all_unshifted", ptm.getScoreAllUnshifted().toString());
                    }
                    if (ptm.getScoreBestPosition() != null) {
                        writer.writeAttribute("score_best_position", ptm.getScoreBestPosition().toString());
                    }
                    if (ptm.getScoreSecondBestPosition() != null) {
                        writer.writeAttribute("score_second_best_position", ptm.getScoreSecondBestPosition().toString());
                    }
                    if (ptm.getScoreShiftedAllPositions() != null) {
                        writer.writeAttribute("score_shifted_all_positions", ptm.getScoreShiftedAllPositions().toString());
                    }
                    if (ptm.getPositionScores() != null) {
                        writer.writeAttribute("position_scores", ptm.getPositionScores());
                    }
                    if (ptm.getShiftedOnlyPositionScores() != null) {
                        writer.writeAttribute("shifted_only_position_scores", ptm.getShiftedOnlyPositionScores());
                    }
                    if (ptm.getShiftedOnlyPositionIons() != null) {
                        writer.writeAttribute("shifted_only_position_ions", ptm.getShiftedOnlyPositionIons());
                    }
                    if (ptm.getIonsBestPosition() != null) {
                        writer.writeAttribute("ions_best_position", ptm.getIonsBestPosition().toString());
                    }
                    if (ptm.getIonsSecondBestPosition() != null) {
                        writer.writeAttribute("ions_second_best_position", ptm.getIonsSecondBestPosition().toString());
                    }
                    if (ptm.getIonsShiftedAllPositions() != null) {
                        writer.writeAttribute("ions_shifted_all_positions", ptm.getIonsShiftedAllPositions().toString());
                    }
                    if (ptm.getScoreShiftedBestPosition() != null) {
                        writer.writeAttribute("score_shifted_best_position", ptm.getScoreShiftedBestPosition().toString());
                    }
                    if (ptm.getIonsShiftedBestPosition() != null) {
                        writer.writeAttribute("ions_shifted_best_position", ptm.getIonsShiftedBestPosition().toString());
                    }
                    if (ptm.getIonsAllUnshifted() != null) {
                        writer.writeAttribute("ions_all_unshifted", ptm.getIonsAllUnshifted().toString());
                    }
                    writer.writeCharacters(System.getProperty("line.separator"));
                }
            }

            writer.writeEndElement(); //corresponding to search_hit
            writer.writeCharacters(System.getProperty("line.separator"));
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }

    private static int CheckModPos(String OrgPepSeq, String NewPepSeq, int OrgPos)
    {
        int index = 99999999;

        if(OrgPepSeq.length()<=NewPepSeq.length()) //missed cleavage peptides
        {
            String ModPep = OrgPepSeq.substring(0, OrgPos);
            if(NewPepSeq.contains(ModPep))
            {
                index = NewPepSeq.indexOf(ModPep) + ModPep.length();
            }
        }
        else  if(OrgPepSeq.length()>NewPepSeq.length()) //semi-tryptic peptides
        {
            int nsIndex = OrgPepSeq.indexOf(NewPepSeq); // the start index of new pep seq in the org pep seq
            int neIndex = nsIndex+NewPepSeq.length();
            int gap = OrgPos-nsIndex;
            index = ((gap<=neIndex) && (gap>0))? gap: 99999999;
        }

        return index;
    }
}
