import os
import sys
import numpy as np

def read_pgfplot_file(file_path):
    """
    Reads a pgfplot printable file and returns its contents as a list of lists.
    """
    data = []
    with open(file_path, 'r') as file:
        for i, line in enumerate(file):
            if line.strip() and not line.startswith('#') and i != 0:
                row = line.strip().split()
                # Drop columns 1 to 4 (indexing starts from 0)
                row = [float(val) for val in row[5:-1]]
                data.append(row)
    return data

def read_pgfplot_files(folder_path):
    """
    Reads pgfplot printable files from a folder and stores the contents
    in a dictionary. Keys are filenames and values are lists of rows.
    """
    data_dict = {}
    for filename in os.listdir(folder_path):
        if filename.endswith('.dat'):
            file_path = os.path.join(folder_path, filename)
            data_dict[filename] = read_pgfplot_file(file_path)
    return data_dict

if __name__ == "__main__":
    if 2 <= len(sys.argv) <= 3:
        folder_path = sys.argv[1]
        reference_file = ""
        if len(sys.argv) == 3:
            reference_file = sys.argv[2]
        else:
            reference_file = os.path.join(folder_path, "3D_R-fp64_S-fp64_SolvePoisson_GMG.dat")
        print(sys.argv[0])
        print(folder_path)
        print(reference_file)
        print()
        pgfplot_data = read_pgfplot_files(folder_path)
        pgfplot_data["reference"] = read_pgfplot_file(reference_file)
        new_data = pgfplot_data.copy()
        reference_data = np.array(pgfplot_data["reference"])
        reference_length = len(reference_data)
        for filename, values in pgfplot_data.items():
            if filename != "reference":
                # Pad the array with NaN values to match the dimension of reference_data
                if len(values) > reference_length:
                    values = values[:reference_length]
                padded_values = np.pad(np.array(values), ((0, reference_length - len(values)), (0, 0)), mode='constant', constant_values=np.nan)
                new_data[filename] = np.divide(padded_values, reference_data)
                with open(os.path.join("./", filename), 'w') as file:
                    file.write("Level Max-Error Point-Wise-Error Relative-Error Point-Wise-Residual Relative-Residual Res-Conv-Rate Runtime Solver-Iterations\n")
                    for lvl, row in enumerate(new_data[filename]):
                        file.write(str(lvl+1) + ' ' + ' '.join(map(str, row)) + '\n')
    else:
        print("Usage: python script_name.py folder_path [relativeFile]")
