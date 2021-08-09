
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

    #remove copyright tag for the base xml
    for entry in base_entry.findall('{http://uniprot.org/uniprot}copyright'):
        base_entry.remove(entry)

    # to remove the outer tag (uniport) for new xml
    # take out the first entry which contains all protein related information and the second entry which is copyright node
    # and append sequentially.
    tobe_added_entry = new_entry[0]
    base_entry.append(tobe_added_entry)
    base_entry.append(new_entry[1])

    base_tree.write(args.result_name)
if __name__ == '__main__':
    main()
