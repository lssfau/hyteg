import numpy as np
import sys
import os
from pathlib import Path
from vtk import vtkXMLUnstructuredGridReader


def read_vtu(filename):
    reader = vtkXMLUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.Update()
    data = reader.GetOutput()
    scalar_data = data.GetPointData().GetArray("res")
    values = []
    for i in range(scalar_data.GetNumberOfTuples()):
        values.append(scalar_data.GetTuple(i)[0])
    return np.array(values)


def count_underflow(data):
    # underflow_map = {
    #     "normal positive":    0,
    #     "subnormal positive": 0,
    #     "flush to zero":      0,
    #     "subnormal negative": 0,
    #     "normal negative":    0,
    # }
    underflow_array = [0., 0., 0., 0., 0., 0., 0., 0]
    # for value in data:
    #     if value >= 6.10e-5:
    #         underflow_map["normal positive"] += 1
    #     elif value >= 5.96e-8:
    #         underflow_map["subnormal positive"] += 1
    #     elif value <= - 6.10e-5:
    #         underflow_map["subnormal negative"] += 1
    #     elif value <= - 5.96e-8:
    #         underflow_map["normal negative"] += 1
    #     else:
    #         underflow_map["flush to zero"] += 1
    #     count += 1
    for value in data:
        if abs(value) >= 6.10e-5:
            underflow_array[0] += 1.
        elif 5.96e-8 <= abs(value) < 6.10e-5:
            underflow_array[1] += 1.
        elif 0. < abs(value) < 5.96e-8:
            underflow_array[2] += 1.
        elif value < - 0.:
            underflow_array[4] += 1.
        elif value <= - 6.10e-5:
            underflow_array[5] += 1.
        elif value <= - 5.96e-8:
            underflow_array[6] += 1.
        else:
            underflow_array[3] += 1.
        underflow_array[-1] += 1

    for i in range(7):
        underflow_array[i] /= underflow_array[-1]

#     print(f'''
#         "normal positive":    {(underflow_array["normal positive"]   /count):.4f}
#         "subnormal positive": {(underflow_array["subnormal positive"]/count):.4f}
#         "flush to zero":      {(underflow_array["flush to zero"]     /count):.4f}
#         "subnormal negative": {(underflow_array["subnormal negative"]/count):.4f}
#         "normal negative":    {(underflow_array["normal negative"]   /count):.4f}
# ''')
#     print(f"{(underflow_array[0]/count):.4f}\t"
#           f"{(underflow_array[1]/count):.4f}\t"
#           f"{(underflow_array[2]/count):.4f}\t"
#           f"{(underflow_array[3]/count):.4f}\t"
#           f"{(underflow_array[4]/count):.4f}\t")
    return underflow_array


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <directory>", file=sys.stderr)
        sys.exit(1)
    directory: Path = Path(sys.argv[1])

    underflow_dictionary = {}
    try:
        for filename in os.listdir(directory):
            if filename.endswith(".vtu"):
                filepath = Path(os.path.join(directory, filename))
                file_information = filepath.stem.split("_")
                level = file_information[0]
                residual_precision = file_information[2][-2:]
                inner_solve_precision = file_information[3][-2:]
                iteration = file_information[-1][2:]
                run_configuration: str = f"R-fp{residual_precision}_I-fp{inner_solve_precision}_level-{level}"
                # print(f"\n{run_configuration}: {iteration}")
                try:
                    scalar_values = read_vtu(filepath)
                    if run_configuration not in underflow_dictionary:
                        underflow_dictionary[run_configuration] = []
                    underflow_dictionary[run_configuration].append(count_underflow(scalar_values))
                except Exception as e:
                    print("An error occurred:", e, file=sys.stderr)
    except Exception as e:
        print("An error occurred:", e, file=sys.stderr)

    for filename, file_content in underflow_dictionary.items():
        Path(os.path.join(directory, "residualDistribution")).mkdir(parents=True, exist_ok=True)
        with open(os.path.join(directory, "residualDistribution", filename + ".dat"), mode='w+') as f:
            f.write("normal_positive\tsubnormal_positive\tflush_to_zero_positive\t"
                    "zero\t"
                    "flush_to_zero_negative\tsubnormal_negative\tnormal_negative\t"
                    "number_of_DoFs")
            for lines in file_content:
                f.write("\n")
                for elements in lines:
                    f.write(f"{elements:.4f}\t")

            f.write("\n")
