#!/bin/bash

# Path to the source file (relative to the script's directory)
SOURCE_FLUID="../../../../../../build/bin/vertical-beam-vertex-vertex-fluid"
SOURCE_SOLID="../../../../../../build/bin/vertical-beam-vertex-vertex-solid"

# Path to the symbolic link
LINK_FLUID="./gismo-executable-dirichlet"
LINK_SOLID="./gismo-executable-neumann"

# Check if the source file exists
if [ -e "$SOURCE_FLUID" ]; then
  # Check if the symbolic link already exists
  if [ -L "$LINK_FLUID" ]; then
    echo "Symbolic link already exists. Updating the link."
    rm "$LINK_FLUID" # Remove the existing symbolic link
  elif [ -e "$LINK_FLUID" ]; then
    echo "A file or directory named '$LINK_FLUID' already exists. Please remove it manually."
    exit 1
  fi

  # Create the symbolic link
  ln -s "$SOURCE_FLUID" "$LINK_FLUID"
  echo "Symbolic link created: $LINK_FLUID -> $SOURCE_FLUID"
else
  echo "Source file '$SOURCE_FLUID' does not exist. Please check the path."
  exit 1
fi


# Check if the source file exists
if [ -e "$SOURCE_SOLID" ]; then
  # Check if the symbolic link already exists
  if [ -L "$LINK_SOLID" ]; then
    echo "Symbolic link already exists. Updating the link."
    rm "$LINK_SOLID" # Remove the existing symbolic link
  elif [ -e "$LINK_SOLID" ]; then
    echo "A file or directory named '$LINK_SOLID' already exists. Please remove it manually."
    exit 1
  fi

  # Create the symbolic link
  ln -s "$SOURCE_SOLID" "$LINK_SOLID"
  echo "Symbolic link created: $LINK_SOLID -> $SOURCE_SOLID"
else
  echo "Source file '$SOURCE_SOLID' does not exist. Please check the path."
  exit 1
fi
