#!/usr/bin/env bash
# download_repo.sh -- generic repository downloader for Module 2.
#
# Resolves ANY ProteomeXchange-family dataset -- PRIDE (PXD/RPXD), jPOST (JPST),
# MassIVE (MSV), iProX (IPX), PDC/CPTAC (PDC) -- OR a direct file-listing URL to its
# raw files, keeps only those matching the vendor extension, and uploads them to S3.
#
# It uses pridepy (>= v0.0.16), which natively downloads PRIDE / MassIVE / jPOST / iProX
# / ProteomeXchange datasets: for jPOST it goes through the PROXI listing endpoint at
# repository.jpostdb.org and transfers from ftp.jpostdb.org. This is why a jPOST entry
# whose *PRIDE mirror* shows nothing (e.g. JPST003393) still resolves here -- we hit the
# jPOST-native path, not the empty PRIDE record.
#
#   VENDOR / EXT env pick the extension (.raw / .d / .wiff); EXT overrides VENDOR.
#   usage: download_repo.sh <accession|URL> <s3://dest/raw> [comma,sep,filenames]
set -euo pipefail

REPO="${1:?repository accession or URL}"
S3_DEST="${2:?s3:// destination prefix}"
FILES="${3:-}"
VENDOR="${VENDOR:-thermo}"
EXT="${EXT:-}"
if [ -z "$EXT" ]; then
    case "$VENDOR" in bruker) EXT=".d";; sciex) EXT=".wiff";; *) EXT=".raw";; esac
fi

# Emit a message to BOTH stderr and stdout: the calling process tees stdout into
# download.summary, so the failure reason reaches the run log and the summary file.
fail() { echo "$1" >&2; echo "$1"; exit 1; }

# ---- ensure pridepy is available (same bootstrap as the PRIDE arm) ----
if ! command -v pridepy >/dev/null 2>&1; then
    echo "pridepy not found; installing for the current user ..." >&2
    python3 -m pip install --user --upgrade pridepy >/dev/null 2>&1 || true
    export PATH="${HOME}/.local/bin:${PATH}"
fi
command -v pridepy >/dev/null 2>&1 \
    || fail "pridepy is required to resolve repositories but could not be installed"

LOCAL="./repo_dl"; mkdir -p "$LOCAL"

count_inputs() {   # vendor-aware count of downloaded inputs
    case "$VENDOR" in
        bruker) find "$1" -maxdepth 4 \( -name '*.d' -o -name '*.d.zip' \) 2>/dev/null | wc -l ;;
        *)      find "$1" -maxdepth 4 -type f -iname "*${EXT}" 2>/dev/null | wc -l ;;
    esac
}

# ---- ProteomeXchange-family accession from a URL or a bare token ----
ACC=$(echo "$REPO" | grep -oiE '(RPXD|PXD|JPST|MSV|IPX|PDC)[0-9]+' | head -1 || true)

if [ -n "$ACC" ]; then
    echo "Resolving ${ACC} via pridepy (PRIDE / jPOST / MassIVE / iProX / PDC) ..." >&2
    case "$ACC" in
        PDC*)
            pridepy download-pdc-files -a "$ACC" -o "$LOCAL" || true
            ;;
        *)
            if [ -n "$FILES" ]; then
                # single-file / demo subset (mirrors --pride_files)
                IFS=',' read -ra WANT <<< "$FILES"
                for f in "${WANT[@]}"; do
                    f=$(echo "$f" | xargs); [ -z "$f" ] && continue
                    echo "  fetching: $f" >&2
                    pridepy download-file-by-name -a "$ACC" -f "$f" -o "$LOCAL" || true
                done
            else
                # whole public RAW collection; retry without --checksum-check for repos
                # (jPOST/MassIVE) that don't always expose checksums
                pridepy download-all-public-raw-files -a "$ACC" -o "$LOCAL" --checksum-check \
                  || pridepy download-all-public-raw-files -a "$ACC" -o "$LOCAL" \
                  || true
            fi
            ;;
    esac
else
    # ---- no PX-family accession: treat REPO as a direct file / listing URL ----
    echo "$REPO" | grep -qiE '^(ftp|https?)://' \
        || fail "No files exist in the provided repository link"
    EXT_RE="${EXT//./\\.}"
    if echo "$REPO" | grep -qiE "${EXT_RE}(\.zip)?([\"?& ]|$)"; then
        URLS=("$REPO")                                   # REPO is itself a file URL
    else
        mapfile -t URLS < <(curl -fsSL "${REPO%/}/" 2>/dev/null \
            | grep -oiE 'href="[^"]+"' | sed 's/href="//; s/"//' \
            | grep -iE "${EXT_RE}(\.zip)?$" | sed "s|^|${REPO%/}/|" | sort -u || true)
    fi
    { [ "${#URLS[@]}" -eq 0 ] || [ -z "${URLS[0]:-}" ]; } \
        && fail "No files exist in the provided repository link"
    for u in "${URLS[@]}"; do
        [ -z "$u" ] && continue
        echo "  $u" >&2
        ( cd "$LOCAL" && { curl -fsSLO "$u" || wget -q "$u"; } ) || echo "  [warn] failed: $u" >&2
    done
fi

N=$(count_inputs "$LOCAL")
[ "$N" -eq 0 ] && fail "No files exist in the provided repository link"

# ---- upload only the vendor files to S3 ----
if [ "$VENDOR" = "bruker" ]; then
    aws s3 sync "$LOCAL/" "${S3_DEST%/}/" --exclude "*" --include "*.d.zip" --include "*.d/*"
else
    aws s3 sync "$LOCAL/" "${S3_DEST%/}/" --exclude "*" --include "*${EXT}"
fi
SIZE=$(du -sh "$LOCAL" 2>/dev/null | cut -f1 || echo "?")
echo "Files downloaded from repository (${ACC:-$REPO}): Total number ${N} and total size ${SIZE}"
