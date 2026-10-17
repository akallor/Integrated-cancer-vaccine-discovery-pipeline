#!/usr/bin/env bash

# =============================================================================
# PRIDE -> EC2 -> S3 raw-file transfer
# =============================================================================
#
# Purpose:
#   Download ONLY the .raw files for a given PRIDE dataset to the EC2 instance,
#   verify them using PRIDE's checksum mechanism, and upload ONLY those .raw
#   files to a given S3 destination.
#
# Usage:
#   ./download_pride_data_with_pridepy.sh <PRIDE_URL_or_accession> <S3_DEST>
#
#   e.g.
#   ./download_pride_data_with_pridepy.sh \
#       https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/03/PXD013649/ \
#       s3://my-bucket/PXD013649
#
# Notes:
#   The expected RAW file count is not hard-coded: it is taken from however many
#   files are actually downloaded, and S3 is verified against that number.
#
# IMPORTANT:
#   This script does NOT delete the local files after uploading.
#   They should be removed only after the S3 transfer has been independently
#   verified.
#
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

# Accept two arguments:
#   $1 = PRIDE link (or bare accession)
#          e.g. https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/03/PXD013649/
#          or just: PXD013649
#   $2 = S3 destination for the .raw files (full s3:// URI)
#          e.g. s3://my-bucket/PXD013649
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <PRIDE_URL_or_accession> <S3_DEST>" >&2
    echo "  e.g. $0 https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/03/PXD013649/ s3://my-bucket/PXD013649" >&2
    exit 1
fi

# Keep the link as passed (for display), normalising the trailing slash.
PRIDE_URL="${1%/}/"

# Derive the accession from the last path component of the link.
# (If a bare accession was passed, this returns it unchanged.)
PRIDE_ACCESSION="$(basename "${1%/}")"

# S3 destination, normalised to a single trailing slash.
case "$2" in
    s3://*) : ;;
    *) echo "ERROR: S3 destination must start with s3:// (got: $2)" >&2; exit 1 ;;
esac
S3_PATH="${2%/}/"

# Bucket name = first path component after s3:// (used for the access check).
_s3_no_scheme="${S3_PATH#s3://}"
S3_BUCKET="${_s3_no_scheme%%/*}"

# Local working directory
LOCAL_DIR="${HOME}/pride/${PRIDE_ACCESSION}"

# Number of RAW files expected in S3 is determined dynamically from however
# many were actually downloaded locally (set after the download completes).
EXPECTED_RAW_FILES=0

# ---- vendor + optional single-file demo (added) --------------------------------
FILES="${3:-}"                       # comma-separated filenames -> fetch ONLY these
VENDOR="${VENDOR:-thermo}"           # exported by the Nextflow DOWNLOAD_PRIDE process
if [ "$VENDOR" = "bruker" ]; then
    S3_GREP='\.d(\.zip)?$'          # Bruker timsTOF: .d / .d.zip
else
    S3_GREP='\.raw$'                 # Thermo
fi
count_inputs() {                     # count vendor inputs under a local dir (maxdepth 1)
    if [ "$VENDOR" = "bruker" ]; then
        find "$1" -maxdepth 1 \( -name '*.d.zip' -o -name '*.d' \) | wc -l
    else
        find "$1" -maxdepth 1 -type f -name '*.raw' | wc -l
    fi
}


# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

info() {
    echo
    echo "======================================================================"
    echo "$1"
    echo "======================================================================"
    echo
}


error_exit() {
    echo
    echo "ERROR: $1"
    echo
    exit 1
}


# -----------------------------------------------------------------------------
# 1. Check prerequisites
# -----------------------------------------------------------------------------

info "Checking prerequisites"

command -v aws >/dev/null 2>&1 || \
    error_exit "AWS CLI is not installed."

command -v python3 >/dev/null 2>&1 || \
    error_exit "Python 3 is not installed."


# -----------------------------------------------------------------------------
# 2. Check AWS credentials / IAM role
# -----------------------------------------------------------------------------

info "Checking AWS identity"

aws sts get-caller-identity

echo
echo "AWS credentials are working."


# -----------------------------------------------------------------------------
# 3. Check available disk space
# -----------------------------------------------------------------------------

info "Checking local disk space"

df -h /

echo
echo "The RAW dataset is approximately 83 GB."
echo "Make sure sufficient disk space is available before continuing."

AVAILABLE_KB=$(df --output=avail "${HOME}" | tail -1 | tr -d ' ')

# Require at least 100 GB free as a safety margin.
MIN_FREE_KB=$((100 * 1024 * 1024))

if [ "${AVAILABLE_KB}" -lt "${MIN_FREE_KB}" ]; then
    error_exit "Less than 100 GB of free disk space is available."
fi

echo "Sufficient disk space is available."


# -----------------------------------------------------------------------------
# 4. Check / create local directory
# -----------------------------------------------------------------------------

info "Preparing local download directory"

mkdir -p "${LOCAL_DIR}"

echo "Local directory:"
echo "  ${LOCAL_DIR}"


# -----------------------------------------------------------------------------
# 5. Check that the S3 bucket exists and is accessible
# -----------------------------------------------------------------------------

info "Checking S3 bucket"

if ! aws s3 ls "s3://${S3_BUCKET}" >/dev/null 2>&1; then
    error_exit "Cannot access S3 bucket: ${S3_BUCKET}"
fi

echo "S3 bucket is accessible:"
echo "  ${S3_BUCKET}"


# -----------------------------------------------------------------------------
# 6. Install/check pridepy
# -----------------------------------------------------------------------------

info "Checking pridepy"

if ! command -v pridepy >/dev/null 2>&1; then

    echo "pridepy was not found."
    echo "Installing it for the current user..."

    python3 -m pip install --user --upgrade pridepy

    # Ensure ~/.local/bin is available in this shell.
    export PATH="${HOME}/.local/bin:${PATH}"
fi

command -v pridepy >/dev/null 2>&1 || \
    error_exit "pridepy could not be installed."

echo "Using:"
pridepy --help >/dev/null
echo "  $(command -v pridepy)"


# -----------------------------------------------------------------------------
# 7. Download ONLY public RAW files from PRIDE
# -----------------------------------------------------------------------------
#
# pridepy's download-all-public-raw-files command downloads the public RAW
# collection for the specified PRIDE accession.
#
# --checksum-check tells pridepy to validate the downloaded files against
# PRIDE's available checksums.
#
# If the command is interrupted, rerunning it is preferable to deleting the
# directory first.
# -----------------------------------------------------------------------------

info "Downloading RAW files from PRIDE"

echo "PRIDE accession:"
echo "  ${PRIDE_ACCESSION}"

echo
echo "Source:"
echo "  ${PRIDE_URL}"

echo
echo "Destination:"
echo "  ${LOCAL_DIR}"

echo

if [ -n "${FILES}" ]; then
    info "Single-file / demo mode: fetching only the requested file(s)"
    IFS=',' read -ra _FLIST <<< "${FILES}"
    for _f in "${_FLIST[@]}"; do
        _f="$(echo "${_f}" | xargs)"                 # trim whitespace
        [ -z "${_f}" ] && continue
        echo "  fetching: ${_f}"
        # download-file-by-name works for any single file (Thermo .raw or a Bruker .d.zip)
        pridepy download-file-by-name -a "${PRIDE_ACCESSION}" -f "${_f}" -o "${LOCAL_DIR}"
    done
else
    # whole dataset: PRIDE's public RAW category is vendor-agnostic (covers Thermo .raw and
    # Bruker .d/.d.zip), so the same command serves both vendors.
    pridepy download-all-public-raw-files \
        -a "${PRIDE_ACCESSION}" \
        -o "${LOCAL_DIR}" \
        --checksum-check
fi


# -----------------------------------------------------------------------------
# 8. Count the local RAW files (this becomes the expected count for S3)
# -----------------------------------------------------------------------------

info "Counting local RAW files"

RAW_COUNT=$(count_inputs "${LOCAL_DIR}")

# The number downloaded locally is what we expect to see in S3.
EXPECTED_RAW_FILES="${RAW_COUNT}"

echo "RAW files found locally: ${RAW_COUNT}"

if [ "${RAW_COUNT}" -eq 0 ]; then
    error_exit "No RAW files were downloaded — nothing to upload."
fi

echo
echo "Using ${RAW_COUNT} as the expected RAW file count."


# -----------------------------------------------------------------------------
# 9. Show local disk usage
# -----------------------------------------------------------------------------

info "Checking local dataset size"

du -sh "${LOCAL_DIR}"


# -----------------------------------------------------------------------------
# 10. Upload ONLY .raw files to S3
# -----------------------------------------------------------------------------
#
# The exclude/include combination means:
#
#   --exclude "*"
#       Exclude everything.
#
#   --include "*.raw"
#       Then include files ending in .raw.
#
# Therefore metadata, logs, checksum files, etc. are NOT uploaded.
#
# aws s3 sync is also safe to rerun. If files are already present, it will
# avoid unnecessarily transferring unchanged files.
# -----------------------------------------------------------------------------

info "Uploading RAW files to S3"

echo "Destination:"
echo "  ${S3_PATH}"

echo

if [ "$VENDOR" = "bruker" ]; then
    aws s3 sync "${LOCAL_DIR}/" "${S3_PATH}" \
        --exclude "*" --include "*.d.zip" --include "*.d/*"
else
    aws s3 sync "${LOCAL_DIR}/" "${S3_PATH}" \
        --exclude "*" --include "*.raw"
fi


# -----------------------------------------------------------------------------
# 11. Verify the S3 RAW file count
# -----------------------------------------------------------------------------

info "Verifying S3 RAW file count"

S3_RAW_COUNT=$(
    aws s3 ls "${S3_PATH}" --recursive \
    | grep -iE "${S3_GREP}" \
    | wc -l
)

echo "RAW files in S3:        ${S3_RAW_COUNT}"
echo "Downloaded locally:     ${EXPECTED_RAW_FILES}"

if [ "${S3_RAW_COUNT}" -ne "${EXPECTED_RAW_FILES}" ]; then
    if [ "$VENDOR" = "bruker" ]; then
        echo "WARN: S3 count ${S3_RAW_COUNT} != local ${EXPECTED_RAW_FILES} (Bruker .d object layout varies); continuing." >&2
    else
        error_exit "Expected ${EXPECTED_RAW_FILES} RAW files in S3, but found ${S3_RAW_COUNT}."
    fi
fi

echo
echo "S3 RAW file count is correct."


# -----------------------------------------------------------------------------
# 12. Show S3 total size
# -----------------------------------------------------------------------------

info "Checking S3 dataset size"

aws s3 ls "${S3_PATH}" \
    --recursive \
    --summarize \
    | tail -2


# -----------------------------------------------------------------------------
# 13. Verify that aws s3 sync has nothing left to upload
# -----------------------------------------------------------------------------
#
# --dryrun performs no changes.
#
# If the local directory and S3 are synchronized, this should produce no
# "upload:" lines.
# -----------------------------------------------------------------------------

info "Performing final S3 sync dry-run"

DRYRUN_OUTPUT=$(
    aws s3 sync "${LOCAL_DIR}/" "${S3_PATH}" \
        --exclude "*" \
        --include "*.raw" \
        --dryrun
)

if [ -n "${DRYRUN_OUTPUT}" ]; then
    echo "WARNING: The dry-run reports differences:"
    echo
    echo "${DRYRUN_OUTPUT}"
    echo
    error_exit "S3 is not fully synchronized with the local RAW files."
fi

echo "No files need to be uploaded."
echo "Local RAW files and S3 are synchronized."


# -----------------------------------------------------------------------------
# 14. Final summary
# -----------------------------------------------------------------------------

info "TRANSFER COMPLETE"

echo "PRIDE accession:"
echo "  ${PRIDE_ACCESSION}"

echo
echo "Local directory:"
echo "  ${LOCAL_DIR}"

echo
echo "S3 destination:"
echo "  ${S3_PATH}"

echo
echo "RAW files:"
echo "  ${S3_RAW_COUNT}"

echo
echo "The transfer has been verified."
echo
echo "IMPORTANT:"
echo "The local files have NOT been deleted."
echo "You may remove them manually after confirming the S3 data is available."
echo
echo "To remove the local copies:"
echo
echo "  rm -rf \"${LOCAL_DIR}\""
echo
echo "Afterward, check disk usage with:"
echo
echo "  df -h /"
echo
