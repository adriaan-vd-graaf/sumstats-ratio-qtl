import pandas as pd
import xml.etree.ElementTree as ET

def parse_hmdb_xml(xml_file, user_ids):
    """
    Parses the HMDB XML file iteratively to find the classification for the user's IDs.

    Args:
        xml_file (str): The path to the local HMDB XML file.
        user_ids (set): A set of HMDB IDs to check against the database.

    Returns:
        dict: A dictionary where keys are HMDB IDs and values are their classification info.
    """
    print("\nParsing the HMDB XML database. This may take a few minutes...")
    classified_metabolites = {}
    # The namespace is required to find tags correctly in the HMDB XML
    namespace = '{http://www.hmdb.ca}'

    # Use iterparse for memory-efficient parsing of the large XML file
    # We listen for the 'end' event on the 'metabolite' tag
    context = ET.iterparse(xml_file, events=('end',))

    processed_count = 0
    found_user_ids = set()

    for event, elem in context:
        # When a </metabolite> tag is reached, process the element
        if elem.tag == f'{namespace}metabolite':
            processed_count += 1
            if processed_count % 20000 == 0:
                print(f"  ...processed {processed_count} metabolites.")

            accession_elem = elem.find(f'{namespace}accession')
            if accession_elem is not None and accession_elem.text in user_ids:
                hmdb_id = accession_elem.text
                found_user_ids.add(hmdb_id)

                # Find the taxonomy elements
                taxonomy_path = f'{namespace}taxonomy'
                super_class_elem = elem.find(f'{taxonomy_path}/{namespace}super_class')
                class_elem = elem.find(f'{taxonomy_path}/{namespace}class')
                sub_class_elem = elem.find(f'{taxonomy_path}/{namespace}sub_class')

                # Store the classification, using 'N/A' if a tag is missing
                classified_metabolites[hmdb_id] = {
                    'super_class': super_class_elem.text if super_class_elem is not None else 'N/A',
                    'class': class_elem.text if class_elem is not None else 'N/A',
                    'sub_class': sub_class_elem.text if sub_class_elem is not None else 'N/A',
                }

                print(f"  -> Classified: {hmdb_id}")

            # Clear the element from memory to keep memory usage low
            elem.clear()

            # Optimization: stop parsing if all user IDs have been found
            if len(found_user_ids) == len(user_ids):
                print("\nAll user-provided IDs have been found. Stopping parse early.")
                break

    print(f"Finished parsing. Total metabolites processed: {processed_count}")
    return classified_metabolites
