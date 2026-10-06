#!/usr/bin/env bash
# Download the Hawkes-Process-Toolkit (Xu, Farajtabar & Zha, ICML 2016) at a
# pinned commit into external/ and apply the small patch used in this project.
#
# The toolkit is distributed under its own licence (GPL-3.0) and is NOT part of
# this repository; it is fetched from the original GitHub repository.
#
# Usage (from the repository root):  bash setup/install_toolkit.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
COMMIT="2548a41c7418b8edef3261ab4479cee4e8eaf071"   # upstream master, 14 Feb 2018
URL="https://github.com/HongtengXu/Hawkes-Process-Toolkit/archive/${COMMIT}.tar.gz"
DEST="${REPO_ROOT}/external/Hawkes-Process-Toolkit"
PATCH="${REPO_ROOT}/setup/patches/thap_empty_sequence_guard.patch"

if [ -d "${DEST}" ]; then
  echo "Toolkit already present in ${DEST} (delete it to reinstall)."
  exit 0
fi

mkdir -p "${REPO_ROOT}/external"
TMP="$(mktemp -d)"
trap 'rm -rf "${TMP}"' EXIT

echo "Downloading Hawkes-Process-Toolkit @ ${COMMIT:0:7} ..."
curl -sSL "${URL}" -o "${TMP}/thap.tar.gz"
tar -xzf "${TMP}/thap.tar.gz" -C "${TMP}"
mv "${TMP}/Hawkes-Process-Toolkit-${COMMIT}" "${DEST}"

echo "Applying patch: $(basename "${PATCH}")"
patch -d "${DEST}" -p1 < "${PATCH}"

echo "Done. MATLAB code adds this folder via matlab/setup_paths.m"
