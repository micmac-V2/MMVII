# syntax=docker/dockerfile:1
FROM ubuntu:24.04

LABEL org.opencontainers.image.title="MMVII build environment"
LABEL org.opencontainers.image.description="Toolchain to build MMVII, its Python API and its PDF documentation"
LABEL org.opencontainers.image.source="https://github.com/micmac-v2/MMVII"
LABEL org.opencontainers.image.licenses="CECILL-B"

# ARG and not ENV: apt must not ask anything while the image is built, but the
# variable has no reason to survive in the containers started from it.
ARG DEBIAN_FRONTEND=noninteractive

# ---- Toutes les dependances necessaires pour compiler MMVII ----
#
# Two compilers are meant to work. clang-18 with libomp-18-dev is the pair the
# GitHub workflow builds MMVII with, so what this image produces matches the CI;
# g++ comes from build-essential and works as an alternative, its OpenMP being
# libgomp. The version of libomp is not a detail: it must match the clang it serves,
# clang -fopenmp looking for its omp.h under its own llvm-<version> directory. With
# a mismatched one, find_package(OpenMP) fails and MMVII is silently built without
# OpenMP -- measured with clang-20 against libomp-dev, which is the 18 on noble.
#
# lld and ccache are not conveniences: CMakeLists.txt uses the first through
# -fuse-ld=lld when the linker accepts it, and installs the second by itself as
# RULE_LAUNCH_COMPILE when it finds it.
#
# The python3/pybind11/doxygen group is what apib11 needs: pybind11-dev for
# find_package(pybind11), python3-dev for Python.h, doxygen for the docstrings
# extracted by makedoc.py, and numpy/setuptools/wheel so that "pip wheel
# --no-build-isolation" works without fetching anything. python3-venv is there
# because Ubuntu 24.04 refuses "pip install" in the system environment (PEP 668):
# install the wheel in a venv, or pass --break-system-packages knowingly.
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    git \
    ccache \
    lld \
    gdb \
    pkg-config \
    libproj-dev \
    proj-bin \
    libgdal-dev \
    gdal-bin \
    libxerces-c-dev \
    libomp-18-dev \
    ca-certificates \
    vim \
    ninja-build \
    clang-18 \
    texlive-latex-base \
    texlive-latex-recommended \
    texlive-latex-extra \
    texlive-fonts-recommended \
    texlive-pictures \
    texlive-science \
    pandoc \
    python3 \
    python3-dev \
    python3-pip \
    python3-venv \
    python3-numpy \
    python3-setuptools \
    python3-wheel \
    pybind11-dev \
    doxygen \
    locales \
    && rm -rf /var/lib/apt/lists/*

# Locales. Without them, every shell of the container greets the user with
# "setlocale: LC_ALL: cannot change locale", because both Docker and
# Singularity/Apptainer pass the locale of the host into a container that has none.
# C.UTF-8 needs no generation, glibc builds it in.
#
# LANG is set to C.UTF-8 rather than to a national locale on purpose: a build must
# not depend on the language of whoever runs it -- the decimal separator, the
# collation order of a sort and the format of the compiler messages all follow it.
# The two national locales are there for whoever asks for them explicitly.
RUN printf 'en_US.UTF-8 UTF-8\nfr_FR.UTF-8 UTF-8\n' > /etc/locale.gen \
 && locale-gen
ENV LANG=C.UTF-8

# Fallback identity, used only when the working directory is root-owned -- that is,
# when nothing is mounted over /src. Ubuntu 24.04 already ships a user at uid 1000,
# so it is replaced rather than added. Singularity/Apptainer ignores all of this:
# it runs as the invoking user, whoever that is.
RUN userdel --remove ubuntu \
 && useradd --create-home --shell /bin/bash --uid 1000 mmvii

# Repertoire de travail : votre depot sera monte ici au lancement du conteneur
WORKDIR /src
RUN chown mmvii:mmvii /src

# Nothing in this image runs as root, and the container still writes in the mounted
# directory with the right owner -- without the caller having to pass --user.
#
# The two requirements pull in opposite directions: writing in a bind mount needs the
# uid of the host owner, which is unknown until the container starts, and a container
# with no --user starts as root. So the entry point starts as root, reads the owner of
# the working directory, and immediately drops to it. Nothing else runs privileged.
#
# This cannot be made absolute: whoever runs the container can always bypass it with
# --entrypoint, or come back as root with "docker exec -u 0". What the image controls
# is its own default, and that default is never root.
#
# Singularity/Apptainer needs none of it -- it never runs as root and binds the home
# and current directories itself -- and the script detects it as the ordinary case of
# already being a non-root user.
COPY <<'EOF' /usr/local/bin/mmvii-entrypoint
#!/bin/bash
# Entry point: drop root, then run the requested command.
set -eu

if [ "$(id -u)" = 0 ] && [ "${MMVII_ALLOW_ROOT:-0}" != 1 ] ; then
    # Owner of the directory the container works in: the bind mount of the user in
    # the normal case, the mmvii user when nothing is mounted over /src.
    uid=$(stat -c %u . 2>/dev/null || echo 0)
    gid=$(stat -c %g . 2>/dev/null || echo 0)
    if [ "$uid" = 0 ] ; then uid=1000 ; gid=1000 ; fi
    # --clear-groups and not --init-groups: the uid of the host has no reason to
    # exist in /etc/passwd here. Re-executing this very script means the setup below
    # is done once, with the final identity.
    exec setpriv --reuid="$uid" --regid="$gid" --clear-groups "$0" "$@"
fi

# HOME is inherited from the identity we no longer have (/root), or from a Docker
# invocation that named none. Under Singularity/Apptainer it is the home of the user,
# bind-mounted and writable, and must be left alone -- hence the test on writability
# rather than an unconditional assignment.
if [ ! -w "${HOME:-}" ] ; then
    export HOME=/tmp
fi

# ccache in the build directory, so that the cache follows the tree it belongs to and
# survives the container without needing a volume of its own. A "rm -rf build" throws
# it away with the rest, which is the intent.
export CCACHE_DIR="${CCACHE_DIR:-$PWD/${MMVII_BUILD_DIR:-build}/ccache}"

exec "$@"
EOF
RUN chmod +x /usr/local/bin/mmvii-entrypoint

# Optional convenience, not the default command: a full build of the MMVII source
# tree found in the current directory. "full" is the complete cycle, not a plain
# build -- it erases the generated symbolic-derivative code, builds a bootstrap MMVII,
# runs GenCodeSymDer and rebuilds. It therefore writes in the *source* tree
# (src/GeneratedCodes and bin/), so the mounted directory must be writable, which is
# what the entry point above is for.
#
# "full" and not "full_P2007": the latter adds libP2007.a, which MMVII does not need
# -- it links the object files directly -- and whose ar takes a long time. Ask for it
# through MMVII_TARGET only to export the library, which is what apib11 links against.
#
#     docker run --rm -v "$PWD:$PWD" -w "$PWD" <image> mmvii-build
#     singularity run <image>.sif mmvii-build
#
# MMVII_BUILD_DIR, MMVII_TARGET and MMVII_JOBS change the build directory, the target
# and the parallelism; the arguments of the script go to the cmake configuration step.
COPY <<'EOF' /usr/local/bin/mmvii-build
#!/bin/bash
# Full build of the MMVII source tree found in the current directory.
set -euo pipefail

if [ ! -f CMakeLists.txt ] || ! grep -q '^project(MMVII' CMakeLists.txt ; then
    echo "mmvii-build: $PWD does not look like an MMVII source tree." >&2
    echo 'Mount your checkout and run from it: -v "$PWD:$PWD" -w "$PWD"' >&2
    exit 2
fi

BUILD_DIR=${MMVII_BUILD_DIR:-build}
TARGET=${MMVII_TARGET:-full}
JOBS=${MMVII_JOBS:-$(nproc)}
export CCACHE_DIR="${CCACHE_DIR:-$PWD/$BUILD_DIR/ccache}"

# Impose the generator on a fresh build directory only: forcing -G on an existing
# CMakeCache.txt configured with another generator is a hard CMake error, so an
# already configured tree (or one configured by hand) keeps working.
GENERATOR=()
if [ ! -f "$BUILD_DIR/CMakeCache.txt" ] ; then
    GENERATOR=(-G Ninja)
fi

set -x
cmake -S . -B "$BUILD_DIR" "${GENERATOR[@]}" "$@"
cmake --build "$BUILD_DIR" -j "$JOBS" --target "$TARGET"
EOF
RUN chmod +x /usr/local/bin/mmvii-build

# The entry point only drops privileges and hands over, so a command line given to
# "docker run" is executed as is and the default below is a plain interactive shell.
# cmake, ninja, git and the two scripts above are run by hand from it.
ENTRYPOINT ["/usr/local/bin/mmvii-entrypoint"]
CMD ["/bin/bash"]
