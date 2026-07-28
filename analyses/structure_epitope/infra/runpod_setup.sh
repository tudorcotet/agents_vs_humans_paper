#!/usr/bin/env bash
# Working pod environment for Boltz-2 v2.2.1 WITH cuEquivariance kernels on RunPod A40 (driver CUDA 12.8).
# Base image: runpod/pytorch:2.4.0-py3.11-cuda12.4.1-devel-ubuntu22.04
# Verified 2026-07-25: torch 2.8.0+cu128 cuda=True, cuEquivariance kernel active, boltz folds OK.
set -euo pipefail
pip install --quiet boltz cuequivariance-torch cuequivariance-ops-torch-cu12 pyyaml gemmi
# pin torch/torchvision to a cu128 build the driver (12.8) supports (pip otherwise pulls torch cu130 -> cuda False)
pip install --quiet --force-reinstall torch==2.8.0 torchvision==0.23.0 --index-url https://download.pytorch.org/whl/cu128
pip install --quiet "numpy<2.2"   # numba needs numpy<=2.1
python3 -c "import torch;assert torch.cuda.is_available();from cuequivariance_torch.primitives.triangle import triangle_multiplicative_update;print('boltz+kernels env OK', torch.__version__)"
