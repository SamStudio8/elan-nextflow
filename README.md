# elan
A nextflow pipeline for quality checking dispersed files and publishing them to a common location.
This repository is under development and you should expect breaking changes at any time.
This pipeline was designed to handle uploaded data for the COG-UK consortium.

## Execute

```
nextflow run elan.nf -c elan.config -resume --manifest <path> --publish <path> --schemegit <path>
```

### Params

* `manifest` location of TSV manifest
* `publish` location to publish artifacts
* `schemegit` location of artic-ncov scheme repository

