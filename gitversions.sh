#!/bin/bash

# Check if a filename was provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <filename>"
    exit 1
fi

filename=$1
dirname="${filename}_versions"
commit_messages_file="$dirname/commit_messages.txt"

# Create the directory if it doesn't exist
mkdir -p "$dirname"

# Initialize the commit messages file
echo -n "" > "$commit_messages_file"

# Retrieve the commit hashes, dates, and messages for the given file
# Use 'END_OF_COMMIT' as a delimiter to identify the end of a commit message
git log --pretty=format:'%H%n%ad%n%B%nEND_OF_COMMIT' --date=iso -- "$filename" | {
    commit_hash=""
    commit_date=""
    commit_message=""
    reading_message=false

    while IFS= read -r line || [[ -n "$line" ]]; do
        if [[ "$line" == "END_OF_COMMIT" ]]; then
            # Process the commit
            # Format the date and time for the filename
            formatted_date=$(echo "$commit_date" | awk '{print $1}')
            version_filename="${formatted_date}_${filename}"

            # Save the file version
            git show "$commit_hash:$filename" > "$dirname/$version_filename"

            # Save the commit message
            echo "[$version_filename]" >> "$commit_messages_file"
            echo "$commit_message" >> "$commit_messages_file"
            echo "" >> "$commit_messages_file"

            # Reset for the next commit
            commit_hash=""
            commit_date=""
            commit_message=""
            reading_message=false
        elif [[ -z "$commit_hash" ]]; then
            commit_hash="$line"
        elif [[ -z "$commit_date" ]]; then
            commit_date="$line"
        else
            if [[ -n "$line" ]]; then
                commit_message+="$line"$'\n'
            fi
        fi
    done
}

echo "All versions of '$filename' have been saved in '$dirname'."
echo "Commit messages have been saved in '$commit_messages_file'."

