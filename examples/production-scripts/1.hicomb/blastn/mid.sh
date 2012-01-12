#!/bin/sh
if test $# -lt 2; then
    echo "$0: insufficient arguments on the command line." >&2
    echo "usage: $0 startlinenum numlines [filename]" >&2
    exit 1
fi
tail -n +$1 $3 | head -n $2
exit $?
