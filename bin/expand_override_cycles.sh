#!/usr/bin/env bash

set -euo pipefail

raw="$1"; shift
read_lens=("$@")

declare -a oc_parts expanded
IFS=';' read -ra oc_parts <<< "$raw"
expanded=()
for i in "${!oc_parts[@]}"; do
    part="${oc_parts[$i]}"
    if [ "$i" -ge "${#read_lens[@]}" ]; then
        echo "ERROR: mask '${raw}' has more segments than the run has reads (${#read_lens[@]}). Check RunInfo.xml and the OverrideCycles template." >&2
        exit 1
    fi
    len="${read_lens[$i]}"
    if [[ "$part" == *'*' ]]; then
        base="${part%\*}"
        used=$(echo "$base" | { grep -oE '[0-9]+' || true; } | awk '{s+=$1}END{print s+0}')
        rest=$((len - used))
        if [ "$rest" -gt 0 ]; then
            expanded+=("${base}${rest}")
        elif [ "$rest" -eq 0 ]; then
            expanded+=("${base%[A-Z]}")
        else
            echo "ERROR: read $((i+1)) has only ${len} cycles but mask '${part}' requires at least ${used} cycles. Check RunInfo.xml and OverrideCycles mask." >&2
            exit 1
        fi
    else
        expanded+=("$part")
    fi
done
(IFS=';'; echo "${expanded[*]}")
