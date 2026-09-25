#!/bin/bash
# mesa_defaults_diff.sh OLD_MESA_DIR NEW_MESA_DIR
# Lists controls removed / added / with changed default values between two MESA installs.
old=$1; new=$2
files="star/defaults/controls.defaults star/defaults/star_job.defaults star/defaults/pgstar.defaults eos/defaults/eos.defaults kap/defaults/kap.defaults"
extract() { # name|value, comments stripped, one per control (array index kept)
  [ -f "$1" ] && sed 's/!.*//' "$1" | grep -E '^ *[A-Za-z_][A-Za-z0-9_]*(\([^)]*\))? *=' | sed -E 's/^ *//; s/ *= */|/; s/ *$//' | sort -u -t'|' -k1,1
}
for f in $files; do
  o=$(mktemp); n=$(mktemp)
  extract "$old/$f" > $o; extract "$new/$f" > $n
  # namelists moved between files (e.g. kap/eos split out of controls) show up as removed here and added there
  echo "=== $f"
  echo "--- removed (in old, not in new):"; join -t'|' -v1 $o $n | cut -d'|' -f1 | tr '\n' ' '; echo
  echo "--- added:";   join -t'|' -v2 $o $n | cut -d'|' -f1 | wc -l | xargs echo "  count:"
  echo "--- default changed (name | old | new):"; join -t'|' $o $n | awk -F'|' '$2!=$3'
  rm -f $o $n
done
for f in star/defaults/history_columns.list star/defaults/profile_columns.list; do
  echo "=== $f: column names removed"
  diff <(sed 's/!.*//' $old/$f | awk '{print $1}' | grep -v '^$' | sort -u) <(sed 's/!.*//' $new/$f | awk '{print $1}' | grep -v '^$' | sort -u) | grep '^<' | tr '\n' ' '; echo
done
