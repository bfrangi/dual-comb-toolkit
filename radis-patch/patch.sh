#!/bin/bash

# Radis patch script
# Radis version: 0.16.2

# The source directory is the directory containing this script
SOURCE_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
DEST_DIR_LBL="$SOURCE_DIR/../.venv/lib/python3.10/site-packages/radis/lbl"
DEST_DIR_API="$SOURCE_DIR/../.venv/lib/python3.10/site-packages/radis/api"

LBL_FILES=(
    "base.py"
    "broadening.py"
    "loader.py"
)
API_FILES=(
    "dbmanager.py"
    "tools.py"
    "hdf5.py"
)

# Exit if the destination directories do not exist
if [ ! -d "$DEST_DIR_LBL" ]; then
    echo "Destination directory for radis/lbl does not exist: $DEST_DIR_LBL"
    echo "Have you created the virtual environment and installed radis?"
    exit 1
fi
if [ ! -d "$DEST_DIR_API" ]; then
    echo "Destination directory for radis/api does not exist: $DEST_DIR_API"
    echo "Have you created the virtual environment and installed radis?"
    exit 1
fi
    
# Move the files from the source to the destination
for file in "${LBL_FILES[@]}"; do
    if [ -f "$SOURCE_DIR/$file" ]; then
        cp "$SOURCE_DIR/$file" "$DEST_DIR_LBL"
        echo "Moved $file to $DEST_DIR_LBL"
    else
        echo "File $file does not exist in the source directory."
    fi
done

for file in "${API_FILES[@]}"; do
    if [ -f "$SOURCE_DIR/$file" ]; then
        cp "$SOURCE_DIR/$file" "$DEST_DIR_API"
        echo "Moved $file to $DEST_DIR_API"
    else
        echo "File $file does not exist in the source directory."
    fi
done


# Print a success message
echo "Radis patch applied successfully."
