"""
Created on Wed Jan 17 2022
This script will generate a csv reporting the number of PTMs and Genes
@author: Taojunfeng Su
"""
import PTM_predictor_module as PTMp
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--query_xml', type=str, default=r'~\..\uniprot-proteome_UP000000625.xml',
                        help='Where is your query xml file?')
    parser.add_argument('--basicInfo_csv', type=str, default=r'~\..\Basic_infor.csv',
                        help='Where do you want to store the csv report')
    args = parser.parse_args()

    def get_PTMs_from_UniprotXML(PathUniprotXML:str, PathBasicInfor:str, Encoding="utf_8", FeatureType=['modified residue', 'lipidation']):

        """
        Parameters
        ----------
        PathUniprotPTM : str, directory of a csv file containing Accession numbers and PTMs info
            DESCRIPTION. The default is PathUniprotPTM.
        Encoding : TYPE, optional
            DESCRIPTION. The default is "utf_8".
        PathBasicInfor : str, directory of the XML file containing number of genes and PTM
            DESCRIPTION. The default is PathUniprotXML.

        Returns
        -------
        None.

        """
        Gene_count = 0
        PTM_count = 0
        for record in PTMp.BioUniprotIO_iterparser(PathUniprotXML):
            if record.id:
                Gene_count += 1
            for ft in record.features:
                if ft.type in FeatureType:
                    PTM_count += 1

        import codecs
        import csv
        with codecs.open(PathBasicInfor, 'w', Encoding) as Basic_Infor:
            BC_writer = csv.writer(Basic_Infor, quoting=csv.QUOTE_MINIMAL)
            BC_writer.writerow(["No. of PTM", PTM_count])
            BC_writer.writerow(["No. of Gene", Gene_count])

    PathUniprotXML = args.query_xml
    PathBasicInfor = args.basicInfo_csv

    get_PTMs_from_UniprotXML(PathUniprotXML, PathBasicInfor, Encoding = "utf_8", FeatureType = ['modified residue', 'lipidation'])


if __name__ == '__main__':
    main()


