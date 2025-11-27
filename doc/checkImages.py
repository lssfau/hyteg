import os
import re
import sys


def check():
    in_folder = [file for file in os.listdir('images')]

    # use this to print the files and copy the output to doxygen.config EXTRA_HTML_FILES
    # for f in sorted(in_folder):
    #     print("@hyteg_SOURCE_DIR@/doc/images/" + f + " \\")

    with open("doxygen.config") as f:
        s = f.read()
    in_docu = re.findall(r"images/(\S*\.png)", s)
    in_docu += re.findall(r"images/(\S*\.jpg)", s)

    for f in in_docu:
        if f not in in_folder:
            print(f + " does not exist in images directory")
            sys.exit(1)

    for f in in_folder:
        if f not in in_docu:
            print(f + " is not in doxygen config. Add it to HTML_EXTRA_FILES. See this script for help")
            sys.exit(1)


if __name__ == "__main__":
    check()
