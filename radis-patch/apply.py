# Radis patch script
# Radis version: 0.16.2
import os
import shutil
import sys
from pathlib import Path

# Check if running on Windows or Linux/MacOS
WINDOWS = "nt" # Windows
POSIX = "posix" # Linux/MacOS
OS = os.name

# Get python packages path depending on OS
path = ""
if OS == WINDOWS:
    path = "..\.venv\Lib\site-packages\\radis\\"
elif OS == POSIX:
    path = "../.venv/lib/python3.10/site-packages/radis/"
else:
    print(f"Unknown OS: {OS}.")
    sys.exit(1)

# The source directory is the directory containing this script
SOURCE_DIR = Path(__file__).resolve().parent
DEST_DIR_LBL = SOURCE_DIR / f"{path}lbl"
DEST_DIR_API = SOURCE_DIR / f"{path}api"

LBL_FILES = ["base.py", "broadening.py", "loader.py"]

API_FILES = ["dbmanager.py", "tools.py", "hdf5.py"]

# Resolve destination directories
DEST_DIR_LBL = DEST_DIR_LBL.resolve()
DEST_DIR_API = DEST_DIR_API.resolve()

# Exit if the destination directories do not exist
if not DEST_DIR_LBL.is_dir():
    print(f"Destination directory for radis/lbl does not exist: {DEST_DIR_LBL}")
    print("Have you created the virtual environment and installed radis?")
    sys.exit(1)

if not DEST_DIR_API.is_dir():
    print(f"Destination directory for radis/api does not exist: {DEST_DIR_API}")
    print("Have you created the virtual environment and installed radis?")
    sys.exit(1)

# Move the files from the source to the destination
for file in LBL_FILES:
    src_file = SOURCE_DIR / file
    if src_file.is_file():
        shutil.copy(src_file, DEST_DIR_LBL)
        print(f"Moved {file} to {DEST_DIR_LBL}")
    else:
        print(f"File {file} does not exist in the source directory.")

for file in API_FILES:
    src_file = SOURCE_DIR / file
    if src_file.is_file():
        shutil.copy(src_file, DEST_DIR_API)
        print(f"Moved {file} to {DEST_DIR_API}")
    else:
        print(f"File {file} does not exist in the source directory.")

# Print a success message
print("Radis patch applied successfully.")
