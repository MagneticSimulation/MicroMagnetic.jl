#!/usr/bin/env bash
# Build the MicroMagnetic.jl CUDA image with a baked sysimage.
#
# The sysimage must be built on a machine with an NVIDIA GPU (kernels are
# baked for that GPU's compute capability). `docker build` cannot access
# GPUs, so this script builds in three steps:
#
#   1. docker build  - dependencies, packages, workload (no GPU needed)
#   2. docker run --gpus all  - build /opt/micromagnetic.so in a container
#   3. docker commit - turn the container into the final image with the
#      sysimage wired up as the default entrypoint
#
# Usage:
#   ./docker/build_cuda.sh [TAG]
#
# Example with registry mirrors:
#   docker build --build-arg CUDA_BASE=docker.m.daocloud.io/nvidia/cuda:12.6.3-devel-ubuntu22.04 ...

set -euo pipefail

TAG=${1:-ghcr.io/magneticsimulation/micromagnetic.jl:latest}
STAGE_TAG=micromagnetic:cuda-sysimage-stage
CONTAINER=micromagnetic-sysimage-build

cd "$(dirname "$0")/.."

echo "==> [1/4] building dependencies image (no GPU needed)"
JULIA_BIN=$(readlink -f "$(command -v julia)")
JULIA_DIR=$(cd "$(dirname "$JULIA_BIN")/.." && pwd)
echo "    using host Julia install: $JULIA_DIR"

# Seed the depot artifacts from the last fully-built final image when
# available (artifacts are content-addressed; identical versions are
# reused), keeping the build immune to mirror hiccups for the multi-GB CUDA
# artifacts. NB: seed from $TAG, not the current stage image — the stage is
# rebuilt by recent (possibly incomplete) attempts.
SEED_DIR=$(mktemp -d /tmp/mm_seed_artifacts.XXXXXX)
if docker image inspect "$TAG" > /dev/null 2>&1; then
    echo "    seeding depot artifacts from $TAG"
    docker create --name mm-seed "$TAG" > /dev/null
    docker cp mm-seed:/usr/local/share/julia/artifacts/. "$SEED_DIR/"
    docker rm mm-seed > /dev/null
fi

# Package/artifact downloads occasionally stall (slow mirror, flaky
# upstream); failed layers are not cached, so retry the whole build.
build_ok=0
for attempt in 1 2 3; do
    # 45 min hard cap per attempt: a hung build (stalled download inside a
    # RUN) never exits on its own, which would bypass the retry logic.
    if timeout 2700 docker build -f docker/Dockerfile.cuda \
        --build-context "julia=$JULIA_DIR" --build-context "seed=$SEED_DIR" \
        -t "$STAGE_TAG" .; then
        build_ok=1
        break
    fi
    rc=$?
    [ "$rc" = 124 ] && echo "    attempt $attempt timed out"
    echo "    attempt $attempt failed, retrying in 30 s"
    sleep 30
done
[ "$build_ok" = 1 ] || { echo "ERROR: docker build failed after 3 attempts"; exit 1; }
rm -rf "$SEED_DIR"

echo "==> [2/4] building sysimage inside a GPU container"
# The probe after the sysimage build exercises the full user path (GPU sim,
# plotting, movie export); runtime artifact downloads it triggers (MKL,
# FFmpeg, ...) are baked into the image by the commit step.
sysimg_ok=0
for attempt in 1 2; do
    docker rm -f "$CONTAINER" 2>/dev/null || true
    if timeout 3600 docker run --name "$CONTAINER" --gpus all "$STAGE_TAG" \
        bash -c "julia /opt/build_sysimage.jl && julia -J /opt/micromagnetic.so /opt/smoke.jl"; then
        sysimg_ok=1
        break
    fi
    echo "    sysimage attempt $attempt failed, retrying in 30 s"
    sleep 30
done
[ "$sysimg_ok" = 1 ] || { echo "ERROR: sysimage build failed after 2 attempts"; exit 1; }

echo "==> [3/4] committing final image"
docker commit \
  --change 'ENTRYPOINT ["julia", "-J", "/opt/micromagnetic.so"]' \
  --change 'CMD []' \
  "$CONTAINER" "$TAG"
# :latest doubles as the CUDA variant; keep the explicit :cuda alias in sync
docker tag "$TAG" "${TAG%:*}:cuda"
docker rm "$CONTAINER"

echo "==> [4/4] smoke test"
docker run --rm --gpus all "$TAG" -e '
    using MicroMagnetic, CUDA
    @assert CUDA.functional()
    println("SMOKE_OK ", CUDA.name(CUDA.device()))
'

echo "==> done: $TAG"
docker images "$TAG"
