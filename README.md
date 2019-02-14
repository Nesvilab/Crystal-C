
# Introduction
Open database search, which enlarges precursor tolerances (e.g. to 500 Da) to include more possible modified peptide candidates, is a practical and efficient approach for elucidating the “dark matter” of proteomics, including chemical and biological peptide modifications. Nevertheless, this strategy considers only unmodified peptide sequences during protein database search, resulting in large mass shifts between precursor neutral masses and identified peptide masses. These mass shifts may correspond to true post-translational modifications (PTMs) or artifacts from missed cleavages, in-source decay, and co-fragmentation (i.e., chimeric spectra). We hereby implemented a bioinformatics tool that detects and corrects possible artifacts, improving the interpretability of the open search results.


# Parameters
* Thread = -1   # Number of threads. "-1" means that Crystal-C automatically uses (thread number - 1) in your computer for processing.
* Fasta = D:\test.fasta   # Protein Fasta File
* RawDataDictionary = D:\test   # The dictionary where the raw data locates
* RawFileExtension = mzML   # The file extension of raw data
* MaxZ = 6   # Maximum precursor charge state
* IsoNum = 3   # Number of theoretical isotope peaks need to be generated
* MassTol = 20   # Precursor mass tolerance
* PrecursorIsolationWindow = 0.7   # Precursor Isolation Window 
* OutputFolder = D:\Test   # The folder for the newly generated pepXML files


# Commands
`java -jar Crystal-C.jar Crystal-C.params PepXMLs`

(For single pepXML)

`java -Xmx8g -jar Crystal-C.jar Crystal-C.params test.pepXML`

(For multiple pepXMLs)

`java -Xmx8g -jar Crystal-C.jar Crystal-C.params *.pepXML`
