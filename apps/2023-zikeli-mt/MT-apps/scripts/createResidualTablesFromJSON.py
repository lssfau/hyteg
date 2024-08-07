import json
import os
import sys
from typing import Union, List, Optional, Any
from pathlib import Path


def generate_tables(json_data: Union[dict, List[dict]], parent_label: Optional[str] = None) -> dict[str:List[List[str]]]:
    """
    Recursively generates tables for each node in the JSON tree.

    Args:
        json_data (dict or list): JSON data representing a node or list of nodes.
        parent_label (str): Label of the parent node. Default is None.

    Returns:
        list: List of tables, where each table is represented as a list of rows.
    """
    tables: dict[str:List[List[str]]] = {}

    # Generate table for current level
    if parent_label is not None:
        for node in json_data:
            entry_name = parent_label + '_' + node['groupLabel']
            if entry_name not in tables:
                tables[entry_name] = []
            x_value = node['label'].split("_")[-1]
            if node['data'] is None:
                node['data'] = 'nan'
            tables[entry_name].append([parent_label, node['label'], node['groupLabel'], str(x_value), str(node['data'])])

    # Generate tables for children
    for child in json_data:
        if 'children' in child:
            tables.update(generate_tables(child['children'], child['label']))

    return tables


def write_tables(output_folder: Path, tables: dict[str:List[List[str]]]) -> None:
    """
    Writes the tables to tab-separated files in the specified output folder.
    Also, creates a new folder with the name of the root node since usually more than one run is done.

    Args:
        output_folder (str): Path to the output folder.
        tables (list): List of tables, where each table is represented as a list of rows.
    """
    for key, table in tables.items():
        table_path: str = os.path.join(output_folder, f"{key}.dat")
        with open(table_path, "w") as f:
            f.write('\t'.join(["parent_label", "label", "group_label", "x-value", "y-value"]) + '\n')
            for row in table:
                f.write('\t'.join(row) + '\n')


def main(json_file: Path, output_folder: Path) -> None:
    """
    Main function to generate tables from a JSON file.

    Args:
        json_file (str): Path to the JSON file containing tree structure.
        output_folder (str): Path to the output folder. Default is "output_pgfplot".
    """

    output_folder = Path.joinpath(output_folder, json_file.stem)

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    with open(json_file, "r") as file:
        json_data: Union[List[dict], dict] = json.load(file)

    tables: List[List[List[str]]] = generate_tables([json_data])
    write_tables(output_folder, tables)


if __name__ == "__main__":
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print("Usage: python script.py <json_file> [output_folder]")
        sys.exit(1)

    json_file: Path = Path(sys.argv[1])
    output_folder: Path = Path("output_pgfplot")
    if len(sys.argv) == 3:
        output_folder = Path(sys.argv[2])

    main(json_file, output_folder)
