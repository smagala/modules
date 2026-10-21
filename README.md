```md
Org: NCEZID
Contact Email: ncezid_shareit@cdc.gov
Status: Active
Keywords: bioinformatics
Version: N/A
Contract#: 47QFCA23F0058
```

# ![amd-ph-core/modules](docs/images/amd_logo.png)

[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

A repository for hosting [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) module files containing tool-specific process definitions and their associated documentation for the AMD Platform for Genomic Surveillance.

## Table of contents

- [Using existing modules](#using-existing-modules)
- [Available modules](#available-modules)
- [Help](#help)
- [Citation](#citation)

## Using existing modules

The module files hosted in this repository define a set of processes for software tools such as `varpipe/coverageanalysis`, `picard/buildbamindex`, `gatk4/mutect2`, etc. This allows you to share and add common functionality across multiple pipelines in a modular fashion.

We use the `nf-core/tools` package with a custom `--git-remote` to manage modules from this repository. This includes using `git` commit hashes to track changes for reproducibility purposes, and to download and install all of the relevant module files.

1. Install the latest version of [`nf-core/tools`](https://github.com/nf-core/tools#installation) (`>=2.0`)
2. List the available modules:

   ```console
   nf-core modules --git-remote https://github.com/amd-ph-core/modules.git list remote
   ```

3. Install a module in your pipeline directory:

   ```console
   nf-core modules --git-remote https://github.com/amd-ph-core/modules.git install varpipe/coverageanalysis
   ```

4. Import the module in your Nextflow script:

   ```nextflow
   #!/usr/bin/env nextflow

   nextflow.enable.dsl = 2

   include { VARPIPE_COVERAGEANALYSIS } from './modules/amd-ph-core/varpipe/coverageanalysis/main'
   ```

5. Remove a module from the pipeline repository:

   ```console
   nf-core modules --git-remote https://github.com/amd-ph-core/modules.git remove varpipe/coverageanalysis
   ```

6. Check that a locally installed module is up-to-date:

   ```console
   nf-core modules --git-remote https://github.com/amd-ph-core/modules.git lint varpipe/coverageanalysis
   ```

## Available modules

### Contamination control

- `clockwork/minimap2` - Minimap2 alignment for contamination detection
- `clockwork/removecontam` - Remove contaminating reads

### BAM processing

- `picard/buildbamindex` - Build BAM index
- `picard/samformatconverter` - Convert between SAM/BAM formats
- `picard/sortsam` - Sort SAM/BAM files

### TB varpipe analysis

- `varpipe/bwamem` - BWA-MEM alignment for TB variant pipeline
- `varpipe/coverageanalysis` - Coverage statistics
- `varpipe/createannotations` - Create variant annotations
- `varpipe/interpretation` - Variant interpretation
- `varpipe/lineage` - Lineage classification
- `varpipe/parseannotations` - Parse variant annotations
- `varpipe/pdf` - Generate PDF report
- `varpipe/report` - Generate analysis report
- `varpipe/structuralvariants` - Structural variant analysis
- `varpipe/tar` - Archive pipeline outputs

### Virus annotation

- `daisribosome` - Annotate supported virus genomes and protein products with DAIS-ribosome

## Available subworkflows

- `fastq_clockwork_decontaminate` - Read decontamination using Clockwork
- `varpipe_core` - Top-level TB varpipe analysis workflow
- `varpipe_processbam` - BAM processing for varpipe (sort, index, mark duplicates)
- `varpipe_variantanalysis` - Variant calling, filtering, and annotation

> **Note:** The tbvarpipe pipeline also uses standard modules from [nf-core/modules](https://github.com/nf-core/modules) (e.g., `bwa/index`, `samtools/*`, `gatk4/*`, `snpeff/snpeff`, `trimmomatic`, `cat/fastq`, `multiqc`). Those are installed separately via `nf-core modules install`.

## Help

If you have any questions or issues please open an issue on this repository.

## Citation

If you use modules from this repository, please cite:

> The nf-core framework for community-curated bioinformatics pipelines.
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> Nat Biotechnol. 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

## Public Domain Standard Notice

This repository constitutes a work of the United States Government and is not
subject to domestic copyright protection under 17 USC § 105. This repository is in
the public domain within the United States, and copyright and related rights in
the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).
All contributions to this repository will be released under the
[MIT No Attribution (MIT-0)](LICENSE) license. By submitting a pull request you
are agreeing to comply with this waiver of copyright interest.

## License Standard Notice

This repository is licensed under the [MIT-0](LICENSE) license.

This source code is free: you can redistribute it and/or modify it under the terms
of the MIT-0 license.

This source code is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the [LICENSE](LICENSE) for more details.

The source code forked from other open source projects will inherit its license.

## Privacy Standard Notice

This repository contains only non-sensitive, publicly available data and
information. All material and community participation is covered by the
[Disclaimer](DISCLAIMER.md).
For more information about CDC's privacy policy, please visit [http://www.cdc.gov/other/privacy.html](https://www.cdc.gov/other/privacy.html).

## Contributing Standard Notice

Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo)
and submitting a pull request. (If you are new to GitHub, you might start with a
[basic tutorial](https://help.github.com/articles/set-up-git).) By contributing
to this project, you grant a world-wide, royalty-free, perpetual, irrevocable,
non-exclusive, transferable license to all users under the terms of the
[MIT-0](LICENSE) license.

Please see the [contributing guidelines](.github/CONTRIBUTING.md) for more information.

All comments, messages, pull requests, and other submissions received through
CDC including this GitHub page may be subject to applicable federal law, including
but not limited to the Federal Records Act, and may be archived. Learn more at
[http://www.cdc.gov/other/privacy.html](https://www.cdc.gov/other/privacy.html).

## Records Management Standard Notice

This repository is not a source of government records, but is a copy to increase
collaboration and collaborative potential. All government records will be
published through the [CDC web site](http://www.cdc.gov).

## Additional Standard Notices

Please refer to [CDC's Template Repository](https://github.com/CDCgov/template) for more information about [contributing to this repository](https://github.com/CDCgov/template/blob/main/CONTRIBUTING.md), [public domain notices and disclaimers](https://github.com/CDCgov/template/blob/main/DISCLAIMER.md), and [code of conduct](https://github.com/CDCgov/template/blob/main/code-of-conduct.md).
