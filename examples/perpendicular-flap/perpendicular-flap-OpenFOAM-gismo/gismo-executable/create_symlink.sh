#!/bin/bash

# Path to the source file (relative to the script's directory)
SOURCE="../../../../../../build/bin/perpendicular-flap-vertex-gismo"

# Path to the symbolic link
LINK="./gismo-executable"

# Check if the source file exists
if [ -e "$SOURCE" ]; then
  # Check if the symbolic link already exists
  if [ -L "$LINK" ]; then
    echo "Symbolic link already exists. Updating the link."
    rm "$LINK" # Remove the existing symbolic link
  elif [ -e "$LINK" ]; then
    echo "A file or directory named '$LINK' already exists. Please remove it manually."
    exit 1
  fi

  # Create the symbolic link
  ln -s "$SOURCE" "$LINK"
  echo "Symbolic link created: $LINK -> $SOURCE"
else
  echo "Source file '$SOURCE' does not exist. Please check the path."
  exit 1
fi
