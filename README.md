# Fast Partial-Sum Update Method for Polar Codes and PAC Codes

This repository implements the proposed **Fast Partial-Sum Update (FPSU)** algorithm for Polar codes and PAC codes.

The project mainly includes:

- FPSU-based Polar Successive Cancellation Decoding, FPSU Polar SCD
- FPSU-based PAC Successive Cancellation Decoding, FPSU PAC SCD
- Conventional Polar and PAC decoding baselines
- Decoding performance comparison
- Total number of partial-sum updates comparison

## Overview

The proposed FPSU algorithm reduces redundant partial-sum updates during decoding by exploiting the structural properties of Polar and PAC codes. This project evaluates both decoding performance and partial-sum update complexity, and compares the proposed FPSU-based decoders with conventional decoding methods.

## Citation

If you find this algorithm or implementation useful for your research, please cite our paper:

```bibtex
@article{li2026fpsu,
  title   = {Fast Partial-Sum Update Method for Polar Codes and PAC Codes},
  author  = {Li, Xun},
  journal = {IEEE Communications Letters},
  year    = {2026},
  doi     = {10.1109/LCOMM.2026.3683509}
}