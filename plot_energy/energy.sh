#!/bin/bash

# Input/output files
LOGFILE="log.lammps"
OUTFILE="energy.dat"

# Find line numbers
start_line=$(grep -n "^ *Step" "$LOGFILE" | tail -1 | cut -d: -f1)
end_line=$(grep -n "^Loop time" "$LOGFILE" | tail -1 | cut -d: -f1)

# Check if both were found
if [[ -z "$start_line" || -z "$end_line" ]]; then
    echo "ERROR: Could not find both 'Step' and 'Loop time' in $LOGFILE"
    exit 1
fi

# Make sure there’s data between them
if (( start_line + 1 >= end_line )); then
    echo "ERROR: No data found between 'Step' and 'Loop time'"
    exit 1
fi

# Extract data (excluding 'Step' and 'Loop time' lines)
sed -n "$((start_line + 1)),$((end_line - 1))p" "$LOGFILE" > "$OUTFILE"

echo "Extracted data from lines $((start_line + 1)) to $((end_line - 1)) into $OUTFILE"

num_atoms=$(grep "^Loop time" log.lammps | tail -1 | awk '{for(i=1;i<=NF;i++) if($i=="atoms") print $(i-1)}')

python plot.py $num_atoms
