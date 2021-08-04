
import argparse
from xml.etree import ElementTree as ET
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--base_xml', type=str, help='Where is your base xml file?')
    parser.add_argument('--new_xml', type=str, help='Where is your name xml file?')
    parser.add_argument('--result_name', type=str, help='What is the name of the result file?')

    args = parser.parse_args()
    # disable uniprot xml prefix
    ET.register_namespace('', "http://uniprot.org/uniprot")
    base_tree = ET.parse(args.base_xml)
    new_tree = ET.parse(args.new_xml)

    base_entry = base_tree.getroot()
    new_entry = new_tree.getroot()
    base_entry.append(new_entry)
    base_tree.write(args.result_name)
if __name__ == '__main__':
    main()