#!/usr/bin/env bash
# convert_raws_from_s3.sh
# Stream each .raw from S3 -> convert -> QC (MS2 present?) -> upload mzML -> delete both.
# Never keeps more than one raw + one mzML on local disk.
set -uo pipefail   # not -e: one bad file shouldn't abort the batch

# ---- .NET / CoreCLR memory caps (fixes 0x8007000E on small boxes; harmless on big ones) ----
export DOTNET_gcServer=0                    # workstation GC: uses a fraction of server-GC memory
export DOTNET_GCHeapHardLimit=0xC0000000    # ~3 GB hard cap on managed heap (hex bytes)
export DOTNET_gcConcurrent=0                 # no background GC heap -> lower startup reservation

# ---- config ----
# Accept the S3 raw input prefix and mzML output prefix as two arguments.
#   $1 = where the .raw files live      e.g. s3://bucket/PXD013649
#   $2 = where the .mzML files go       e.g. s3://bucket/PXD013649/mzml_v2
if [ "$#" -lt 2 ]; then
  echo "Usage: $0 <S3_RAW_PREFIX> <S3_MZML_PREFIX>" >&2
  echo "  e.g. $0 s3://bucket/PXD013649 s3://bucket/PXD013649/mzml_v2" >&2
  exit 1
fi
S3_RAW_PREFIX="${1%/}"    # strip any trailing slash
S3_MZML_PREFIX="${2%/}"   # strip any trailing slash
WORKDIR="${WORKDIR:-/tmp/raw_convert}"
LOG="${WORKDIR}/convert.log"
MIN_MS2="${MIN_MS2:-1}"
RECONVERT="${RECONVERT:-0}"   # set to 1 to ignore existing mzMLs and overwrite

mkdir -p "$WORKDIR"/{raw,mzml}
: > "$LOG"
say(){ echo "[$(date +%H:%M:%S)] $*" | tee -a "$LOG"; }

say "GC config: gcServer=$DOTNET_gcServer heapLimit=$DOTNET_GCHeapHardLimit"
say "Free memory at start:"; free -h | tee -a "$LOG"

VENDOR="${VENDOR:-thermo}"
EXT="${EXT:-}"
if [ -z "$EXT" ]; then case "$VENDOR" in bruker) EXT=".d";; sciex) EXT=".wiff";; *) EXT=".raw";; esac; fi
say "vendor = ${VENDOR}  ext = ${EXT}"
if [ "$VENDOR" = "bruker" ]; then LISTPAT='\.d\.zip$'; else LISTPAT="${EXT//./\\.}$"; fi
say "Listing ${VENDOR} (${EXT}) under ${S3_RAW_PREFIX}/"
mapfile -t KEYS < <(aws s3 ls "${S3_RAW_PREFIX}/" | awk '{print $NF}' | grep -iE "$LISTPAT")
N=${#KEYS[@]}
say "Found ${N} input file(s)."
[ "$N" -eq 0 ] && { say "No input files found — check the prefix / --vendor (${VENDOR})."; exit 1; }

# MS2-present QC, tolerant of cvParam attribute order (Thermo vs tdf2mzml)
count_ms2(){ grep -cE '(value="2"[^/>]*name="ms level")|(name="ms level"[^/>]*value="2")' "$1" 2>/dev/null || echo 0; }

ok=0; bad=0; skipped=0
for i in "${!KEYS[@]}"; do
  key="${KEYS[$i]}"
  if [ "$VENDOR" = "bruker" ]; then base="${key%.d.zip}"; else base="${key%$EXT}"; fi
  in_local="${WORKDIR}/raw/${key}"
  mzml_local="${WORKDIR}/mzml/${base}.mzML"
  mzml_s3="${S3_MZML_PREFIX}/${base}.mzML"
  say "[$((i+1))/${N}] ${key}"

  if [ "$RECONVERT" -ne 1 ] && aws s3 ls "$mzml_s3" >/dev/null 2>&1; then
    say "   already converted in S3 — skipping"; skipped=$((skipped+1)); continue
  fi

  if ! aws s3 cp "${S3_RAW_PREFIX}/${key}" "$in_local" 2>>"$LOG"; then
    say "   ERROR: download failed; df:"; df -h "$WORKDIR" | tee -a "$LOG"
    bad=$((bad+1)); rm -rf "$in_local"; continue
  fi

  conv_ok=1
  if [ "$VENDOR" = "bruker" ]; then
    # unzip sample.d.zip -> sample.d, then convert with tdf2mzml (Bruker timsTOF)
    if ! unzip -q -o "$in_local" -d "${WORKDIR}/raw/" 2>>"$LOG"; then
      say "   ERROR: unzip failed"; conv_ok=0
    fi
    ddir="${WORKDIR}/raw/${base}.d"
    [ -d "$ddir" ] || ddir="$(find "${WORKDIR}/raw" -maxdepth 2 -type d -name '*.d' | head -1)"
    if [ "$conv_ok" -eq 1 ] && { [ -z "$ddir" ] || [ ! -d "$ddir" ]; }; then
      say "   ERROR: no .d directory after unzip"; conv_ok=0
    fi
    if [ "$conv_ok" -eq 1 ]; then
      # tdf2mzml (https://github.com/mafreitas/tdf2mzml); adjust flags to your build if needed
      if ! tdf2mzml.py -i "$ddir" -o "$mzml_local" 2>>"$LOG"; then
        say "   ERROR: tdf2mzml failed — last log:"; tail -4 "$LOG" | sed 's/^/       /'; conv_ok=0
      fi
    fi
    rm -rf "$ddir"
  elif [ "$VENDOR" = "sciex" ]; then
    # Sciex .wiff -> msconvert (ProteoWizard + vendor libs). Needs the .wiff.scan companion too.
    if ! msconvert "$in_local" --mzML --outdir "${WORKDIR}/mzml/" 2>>"$LOG"; then
      say "   ERROR: msconvert failed — last log:"; tail -4 "$LOG" | sed 's/^/       /'; conv_ok=0
    fi
  else
    if ! ThermoRawFileParser -i "$in_local" -o "${WORKDIR}/mzml/" -f 2 2>>"$LOG"; then
      say "   ERROR: conversion failed — last log lines:"; tail -4 "$LOG" | sed 's/^/       /'; conv_ok=0
    fi
  fi
  if [ "$conv_ok" -ne 1 ]; then bad=$((bad+1)); rm -rf "$in_local" "$mzml_local"; continue; fi

  ms2=$(count_ms2 "$mzml_local")
  if [ "${ms2:-0}" -lt "$MIN_MS2" ]; then
    say "   ERROR: MS2=${ms2} (<${MIN_MS2}) — NOT uploading, flagged bad"
    bad=$((bad+1)); rm -rf "$in_local" "$mzml_local"; continue
  fi
  say "   OK: MS2=${ms2}"

  if aws s3 cp "$mzml_local" "$mzml_s3" 2>>"$LOG"; then
    say "   uploaded -> ${mzml_s3}"; ok=$((ok+1))
  else
    say "   ERROR: upload failed (keeping local mzML)"; bad=$((bad+1)); rm -rf "$in_local"; continue
  fi
  rm -rf "$in_local" "$mzml_local"
done

say "DONE. converted=${ok}  skipped=${skipped}  bad=${bad}  of ${N}"
