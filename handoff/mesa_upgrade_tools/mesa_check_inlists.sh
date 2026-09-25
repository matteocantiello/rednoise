#!/bin/bash
# mesa_check_inlists.sh NEW_MESA_DIR inlist [inlist ...]
# Controls set in the given inlists that do not exist in the new version's defaults files.
new=$1; shift
known=$(cat $new/star/defaults/{controls,star_job,pgstar}.defaults $new/eos/defaults/eos.defaults \
            $new/kap/defaults/kap.defaults $new/star/defaults/*_dev.defaults 2>/dev/null \
        | sed 's/!.*//' | grep -oE '^ *[A-Za-z_][A-Za-z0-9_]*' | tr -d ' ' | sort -u)
for f in "$@"; do
  sed 's/!.*//' "$f" | grep -oE '^ *[A-Za-z_][A-Za-z0-9_]*(\([^)]*\))? *=' | sed -E 's/\(.*//; s/[ =]//g' | sort -u \
  | comm -23 - <(echo "$known") | sed "s#^#$f: unknown in new MESA: #"
done
