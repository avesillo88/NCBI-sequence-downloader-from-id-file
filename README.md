# NCBI sequence downloader from ID file

## Description

The aim of this programme is to download amino acid or nucleotide sequences from a list of identifiers stored in a file. For this purpose, the NCBI E-utilities Efecth and Epost are used.  

It is coded in an attempt to extend the work of [NCBI Mass Sequence Downloader](https://github.com/StuntsPT/NCBI_Mass_Downloader), which does not support the use of files with identifiers to be downloaded.

To carry out this programme, the Epost E-utility is first used in a general way until all sequences are downloaded or an error arises. If an error is received, a second download attempt is started with the E-utility EFetch.

This is because during the creation of the code, I encountered different problems:

1. **Epost does not always return all requested sequences**. For example, if you ask for 200 streams, depending on the time of day, the server may return a lower number of sequences. Therefore, at the end of the download, there are IDs whose sequences have not been loaded into the output file. As long as all the sequences have not been downloaded, the programme continues to update the IDs that remain to be downloaded.

2. **Some sequences cannot be downloaded by Epost but can be downloaded by Efetch**. EPost works with a smaller number of databases than EFetch. For example, ['EAI4316178.1'](https://www.ncbi.nlm.nih.gov/search/all/?term=EAI4316178.1) protein is found in the 'Identifical Protein Group' database. If you try to download this sequence with EPost, it will return an empty response, whereas EFetch will respond with the correct sequence. Thus, if an empty response is perceived from Epost, the download shall be attempted with EFetch. This shall only happen when only the proteins that cannot be downloaded with EPost are left.

## Dependencies
- Python 3
- [Requests](https://github.com/psf/requests)

## References

Sayers E. The E-utilities In-Depth: Parameters, Syntax and More. 2009 May 29 [Updated 2022 Nov 30]. In: Entrez Programming Utilities Help [Internet]. Bethesda (MD): National Center for Biotechnology Information (US); 2010-. Available from: https://www.ncbi.nlm.nih.gov/books/NBK25499/
