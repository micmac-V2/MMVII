# MMVII runtime image, based on Ubuntu.
#
# This Dockerfile does *not* compile MMVII: it packages the binary and the data
# directories already produced by the Linux build job. Its build context is the
# staging directory prepared by .github/workflows/build_mmvii.yml:
#
#     <context>/Dockerfile
#     <context>/MMVII/bin/MMVII
#     <context>/MMVII/MMVII-RessourceDir/
#     <context>/MMVII/MMVII-LocalParameters/
#     <context>/MMVII/MMVII-TestDir/
#     <context>/MMVII/bash-completion/
#     <context>/runtime-packages.txt
#
# The binary is dynamically linked, so it must have been built on the same Ubuntu
# release as the base image. The workflow therefore passes the release of its own
# runner as UBUNTU_VERSION, and the image follows it instead of being frozen here.
# The default below only serves a manual build, and the smoke test at the end of
# this file fails the image build rather than producing a container dying at run
# time. The runtime packages are derived from the binary by that same job, so no
# release-specific package name is written here either.
#
# Usage: the container runs a shell, so that a whole session of MMVII commands
# costs one container instead of one per command, and so that the completion is
# available. Giving a command line runs it directly instead.
#     docker run --rm <image> MMVII                        # list the MMVII commands
#     docker run --rm <image> MMVII Bench 1                # run the benches
#
# Working on the files of the user needs three things, and --user is not a detail:
# without it the container runs as the uid below, which cannot write in a bind mount
# unless the host account happens to have that same uid, and MMVII dies on its log
# file. Mounting under the *host* paths then makes an absolute path designate the same
# file on both sides. See README.md, which documents this for the user.
#     docker run --rm -it --user "$(id -u):$(id -g)" -e HOME="$HOME" \
#         -v "$HOME:$HOME" -v "$PWD:$PWD" -w "$PWD" <image> [MMVII <cmd> <args>]
#
# Singularity/Apptainer needs none of that -- it runs as the invoking user and binds
# the current and home directories itself:
#     singularity run mmvii.sif                            # a shell, with completion
#     singularity run mmvii.sif MMVII <cmd> <args>         # a single command

ARG UBUNTU_VERSION=24.04
FROM ubuntu:${UBUNTU_VERSION}

LABEL org.opencontainers.image.title="MMVII"
LABEL org.opencontainers.image.description="MMVII photogrammetry command line tool"
LABEL org.opencontainers.image.source="https://github.com/micmac-v2/MMVII"
LABEL org.opencontainers.image.licenses="CECILL-B"

# Runtime packages: the list is not written here because it is release-specific
# (libgdal34t64, libproj25, libxerces-c3.2t64... on 24.04). The build job derives
# it from the very binary being packaged, with ldd and dpkg -S, so it names the
# libraries MMVII is actually linked against -- the OpenMP runtime of the compiler
# that built it included, without having to install both. gdal-data and proj-data
# arrive as dependencies of the GDAL and PROJ packages, so their resource files
# are present.
#
# The few packages named here are not needed by MMVII itself, they serve the
# interactive shell the container runs: bash-completion and python3 for the command
# completion (see below), less, vim and nano to read and edit the files of a project,
# file, tree, unzip and xz-utils to inspect and unpack its data, curl, wget and the
# ssh client to bring that data in. ca-certificates, already required, is what makes
# their TLS work.
COPY runtime-packages.txt /tmp/runtime-packages.txt
RUN apt-get update \
 && DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends --yes \
        bash-completion \
        ca-certificates \
        curl \
        file \
        less \
        nano \
        openssh-client \
        python3 \
        tree \
        unzip \
        vim \
        wget \
        xz-utils \
        $(cat /tmp/runtime-packages.txt) \
 && rm -rf /var/lib/apt/lists/* /tmp/runtime-packages.txt

# Ubuntu 24.04 already ships a user with uid 1000: replace it by an MMVII one,
# so that a bind-mounted directory of a standard Linux account stays writable.
RUN userdel --remove ubuntu \
 && useradd --create-home --shell /bin/bash --uid 1000 mmvii \
 && mkdir /data \
 && chown mmvii:mmvii /data

# MMVII deduces the root of its installation from the canonical path of its own
# executable (<directory of the executable>/..), so this layout is a contract:
# bin/ and the data directories must remain siblings.
#
# Left owned by root, hence read only for the user running the container. Nothing
# writes in the installation: the benches work in the system temporary directory and
# the profile under $HOME. The only exception, EditCamDataBase asked for its LocalMVVI
# or GlobalMMVII database, cannot work here and could not anyway: this image is also
# meant to be converted to a Singularity/Apptainer SIF, whose filesystem is read only
# and which runs as the invoking user, so no ownership inside the image would help.
# Being root owned makes that constraint explicit under Docker too, and turns a
# regression writing in the installation into a visible failure -- the image build
# runs MMVII after USER mmvii below, and the CI runs the benches in the container.
COPY MMVII /opt/MMVII

# Resolving the symlink gives back /opt/MMVII/bin/MMVII, so the installation
# root is still found when MMVII is invoked by its short name.
RUN ln -s /opt/MMVII/bin/MMVII /usr/local/bin/MMVII

# Gives the *host* shell an MMVII command and its completion, without installing
# anything there. Run without argument, it prints the command the user has to type;
# it serves Docker and Singularity/Apptainer, detecting which one runs it.
RUN ln -s /opt/MMVII/bash-completion/mmvii-bash-init /usr/local/bin/mmvii-bash-init

# Name under which this image is published, the only way mmvii-bash-init can name it
# in the MMVII function it generates: a Docker container cannot know the image it started
# from, its hostname being the container id. Singularity/Apptainer needs nothing here,
# the path of the SIF being exposed in the container. Empty for a manual build, which
# then falls back on an argument or on a placeholder.
ARG MMVII_IMAGE=
ENV MMVII_IMAGE=${MMVII_IMAGE}

# Command completion for an interactive shell. The generated mmvii-completion is
# self-contained -- make_completion.sh inlines the python script in it as a here
# document -- so it can be linked from anywhere. /etc/bash_completion.d/ is read
# by the bash-completion package for every user, root included, which is why the
# completion is not installed in the .bashrc of the mmvii user only. That directory
# is read by bash-completion but created by no package, hence the mkdir.
RUN mkdir -p /etc/bash_completion.d \
 && ln -s /opt/MMVII/bash-completion/mmvii-completion /etc/bash_completion.d/mmvii

USER mmvii
WORKDIR /data
ENV HOME=/home/mmvii

# Fail the image build on a missing or incompatible shared library, rather than
# at the first run of the container.
RUN ! ldd /opt/MMVII/bin/MMVII | grep 'not found'

# Smoke test: lists the available commands. It also creates the user profile
# (~/.config/MMVII) now, so that no first-run write is needed later.
RUN MMVII > /dev/null

# A command line given to docker run replaces this CMD and is executed as is, so
# both a shell session and a one-shot "docker run <image> MMVII <cmd>" work. An
# ENTRYPOINT would forbid the latter: bash would take the command as a script to
# read. The shell is interactive only under -it, which is what makes it read the
# .bashrc bringing in the completion.
CMD ["/bin/bash"]
