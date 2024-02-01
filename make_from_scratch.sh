#!/bin/bash

# Define the list of files
FILES_TO_REMOVE="CFPAssoc.m
CFPs.m
CFPTable.m
fnTerms.m
chenDeltas.m
ThreeBodyTables.m
ThreeBodyTable.m
SpinSpinTable.m
SOOandECSOTable.m
ElectrostaticTable.m
ReducedSOOandECSOLSTable.m
ReducedT22Table.m
ReducedV1kTable.m
SpinOrbitTable.m
ReducedUkTable.m
CrystalFieldTable_f1.zip
CrystalFieldTable_f2.zip
CrystalFieldTable_f3.zip
CrystalFieldTable_f4.zip
CrystalFieldTable_f4.zip"

# Directory from which files will be removed
DIRECTORY="./data"

# Function to confirm deletion
confirm_deletion() {
    echo "Are you sure you want to delete the following files from $DIRECTORY?"
    for FILE in $FILES_TO_REMOVE; do
        echo "$FILE"
    done

    read -p "Type 'yes' to proceed: " CONFIRM
    if [ "$CONFIRM" = "yes" ]; then
        return 0
    else
        echo "Deletion cancelled."
        return 1
    fi
}

# Call the confirmation function
if confirm_deletion; then
    # Iterate through the file list and remove each one
    cd "$DIRECTORY"
    for FILE in $FILES_TO_REMOVE; do
        # Check if the file exists
        if [ -f "$FILE" ]; then
            echo "Removing $FILE"
            rm "$FILE"
        else
            echo "File not found: $FILE"
        fi
    done

    echo "File removal completed."
else
    echo "Operation aborted."
fi
