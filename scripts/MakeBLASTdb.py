"""
Created on Wed Aug 04 2021

@author: Taojunfeng Su
"""

#!/usr/bin/python
import argparse
import PTM_predictor_module as PTMp


# Example usage
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', type=str, help='Where is your fasta file?')
    parser.add_argument('--result_name', type=str, help='What is the name of the result file? Put the species name is recommended')

    args = parser.parse_args()
    PTMp.create_blast_database(input_fasta = args.fasta, db_type="prot", db_name=args.result_name)
if __name__ == '__main__':
    main()