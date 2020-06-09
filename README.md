
# Introduction
  Shotgun proteomics using liquid chromatography coupled to mass spectrometry (LC-MS) is commonly used to identify peptides containing post-translational modifications. With the emergence of fast database search tools such as MSFragger, the approach of enlarging precursor mass tolerances during the search (termed “open search”) has been increasingly used for comprehensive characterization of post-translational and chemical modifications of protein samples. However, not all mass shifts detected using the open search strategy represent true modifications, as artifacts exist from sources such as unaccounted missed cleavages or peptide co-fragmentation (chimeric MS/MS spectra). Here, we present Crystal-C, a computational tool that detects and removes such artifacts from open search results. Our analysis using Crystal-C shows that, in a typical shotgun proteomics data set, the number of such observations is relatively small. Nevertheless, removing these artifacts helps to simplify the interpretation of the mass shift histograms, which in turn should improve the ability of open search-based tools to detect potentially interesting mass shifts for follow-up investigation.


# Parameters
* thread = -1                             # Number of threads. "-1" means that Crystal-C automatically uses (total number of threads - 1) in your computer for processing. <br />
* fasta = D:\test.fasta                   # Protein Fasta File <br />
* raw_file_location = D:\test             # The dictionary where the raw data locates
* raw_file_extension = mzML               # The file extension of raw data
* output_location = D:\Test               # The folder for the newly generated pepXML files

* precursor_charge = 1 6                  # The precursor charge state range
* isotope_number = 3                      # Number of theoretical isotope peaks need to be generated
* precursor_mass = 20                     # Precursor mass tolerance
* precursor_isolation_window = 0.7        # Precursor Isolation Window 



# How to Download
Download the latest version [here](https://github.com/Nesvilab/Crystal-C/releases/latest)


# How to Cite
•	Chang HY, Kong AT, da Veiga Leprevost F, Avtonomov DM, Haynes SE, Nesvizhskii AI. Crystal-C: A Computational Tool for Refinement of Open Search Results. J Proteome Res. 2020 ([Manuscript](https://pubs.acs.org/doi/abs/10.1021/acs.jproteome.0c00119)).

For other tools developed by the Nesvizhskii lab, see our website www.nesvilab.org.


# Commands
`java -jar Crystal-C.jar Crystal-C.params PepXMLs`

(For single pepXML)

`java -Xmx8g -jar Crystal-C.jar Crystal-C.params test.pepXML`

(For multiple pepXMLs)

`java -Xmx8g -jar Crystal-C.jar Crystal-C.params *.pepXML`
