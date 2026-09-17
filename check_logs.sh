#!/bin/bash
# Summarise recent pipeline job logs.
# Usage: bash check_logs.sh [N_hours]   (default: last 24 hours)
# Prints a PASS/FAIL line per log file plus the first error found in failing jobs.

HOURS=${1:-24}
MINUTES=$(( HOURS * 60 ))
LOG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/logs"

# Patterns that indicate a fatal problem
ERROR_PAT='EXITING because of FATAL|FATAL ERROR|RuntimeError:|error: \*\*\* JOB|CANCELLED AT|Traceback \(most recent|^Error in |slurmstepd.*error'

# Pattern that the job scripts print on clean completion
OK_PAT='\] .* done$'

echo "Logs modified in the last ${HOURS}h"
echo "------------------------------------------------------------"

pass=0; fail=0; running=0
while IFS= read -r f; do
    name=$(basename "$f")
    first_error=$(grep -m1 -E "$ERROR_PAT" "$f" 2>/dev/null || true)
    has_ok=$(grep -cE "$OK_PAT" "$f" 2>/dev/null || echo 0)

    if [[ -n "$first_error" ]]; then
        echo "FAIL     $name"
        echo "         $(echo "$first_error" | sed 's/^[[:space:]]*//')"
        (( fail++ ))
    elif [[ "$has_ok" -gt 0 ]]; then
        echo "PASS     $name"
        (( pass++ ))
    else
        echo "RUNNING? $name"
        (( running++ ))
    fi
done < <(find "$LOG_DIR" -name "*.out" -mmin "-${MINUTES}" | sort)

echo "------------------------------------------------------------"
echo "PASS: $pass   FAIL: $fail   RUNNING/UNKNOWN: $running"
