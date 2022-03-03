<div align="center">
<p align="center">
   <img src="/images/inbound-ops-patch.png" alt="Inbound Data Operations Badge" width="250"> 
</p>
<h1 align="center">Elan</h1>
</div>

Elan is the inbound distribution pipeline for CLIMB-COVID.
Elan is a Nextflow DSL2 pipeline for quality checking dispersed files and publishing them to a common location.

Elan was authored by [@samstudio8](https://github.com/samstudio8) as part of CLIMB-COVID's work with the COVID-19 Genomics UK Consortium.
Elan is now maintained by [@BioWilko](https://github.com/BioWilko).

## Parameters and environment variables

### Controlling Elan

#### go-full-elan.sh

| Name | Description |
| ---- | ----------- |
| `DATESTAMP` | YYYYMMDD datestamp to identify today's run |
| `UPLOADS_DIR_GLOB` | Uploads glob to expand as part of `resolve_uploads` |
| `ELAN_CONFIG` | Path to current Nextflow configuration |
| `ELAN_SOFTWARE_DIR` | Path to local clone of elan-nextflow |
| `ELAN_RUN_DIR` | Path to dir to run Elan from (scratch) |
| `ELAN_DIR` | Path to CLIMB-COVID staged artifacts root (nicholsz/), passed as `--publish` to elan-nf |
| `ARTIFACTS_ROOT` | Path to new CLIMB-COVID published artifact root (/artifacts/), passed as `--artifacts_root` to elan-nf |
| `NEXTFLOW_BIN` | Path to nextflow binary |
| `SLACK_MGMT_HOOK` | Slack HTTPS webhook for posting debug messages |
| `SLACK_REAL_HOOK` | Slack HTTPS webhook for posting inbound-dist messages |
| `MQTT_HOST` | IP for MQTT broker |
| `MQTT_ENV` | Root MQTT topic (`CLIMBDEV` or `COGUK`) |
| `CONDA_OCARINA` | conda prefix to `conda activate` when performing ocarina calls outside of Elan |
| `CONDA_IPC` | conda prefix to `conda activate` when sending MQTT messages with Tael |

`go-full-elan.sh` will immediately terminate with exit 64 (`EX_USAGE`) if any of the listed parameters are missing from the environment.

#### cog-publish.sh

Note these variables are checked inside go-full-elan.sh as it is the main entrypoint.
Additionally, variables defined above may be used in cog-publish.sh without listing them below.

| Name | Description |
| ---- | ----------- |
| `COG_PUBLISH_MODE` | Set to `local` or `slurm` to control how the daily consensus is generated |
| `CONDA_POSTELAN` | conda prefix to `conda activate` for publish related activities |


### Running Elan

#### elan-nextflow parameters

| Name | Description |
| ---- | ----------- |
| `--mode` | `inbound` or `ocarina` |
| `--ocarina_profile` | Ocarina profile to use for `save_manifest` (inbound) or `play_ocarina` (ocarina) |
| `--datestamp` (elan) | YYYYMMDD datestamp to identify today's run |
| `--uploads` (elan) | Glob path for CLIMB-COVID user uploads (ensure to quote appropriately to prevent premature glob expansion) |
| `--publish` (elan) | Path to CLIMB-COVID staged artifacts root (nicholsz/) |
| `--artifacts_root` (elan) | Path to new CLIMB-COVID published artifact root (/artifacts/) |
| `--minlen` (elan) | Minimum genome size required to pass the `screen_uploads` step [int] |
| `--manifest` (ocarina) | Path to Ocarina manifest created by Elan pipeline |

#### elan-nextflow environment variables

| Name | Description |
| ---- | ----------- |
| `ELAN_SLACK_MGMT_HOOK` | HTTPS hook for posting management and control messages to Slack |
| `ELAN_SLACK_INBOUND_HOOK` | HTTPS hook for posting counts and QC messages to Slack |
| `OCARINA_CONF_FILE` | Path to Ocarina JSON configuration |
| `OCARINA_PROFILE` | Profile to load from Ocarina JSON configuration |

Note that Elan will only error if `MAJORA_DOMAIN` is unset, all other `MAJORA_*` variables are not checked.

#### Nextflow environment variables

| Name | Description |
| ---- | ----------- |
| `NXF_WORK` | Path to NXF working dir (basename must be `nxf_work` on CLIMB-COVID) |
| `NXF_CONDA_CACHEDIR` | Path to conda cache dir (basename must be `.conda` on CLIMB-COVID) |
| `NXF_DEBUG` | One of: 1,2,3. See https://www.nextflow.io/docs/latest/config.html |


## Invocation

Add the following line to the execution node's `crontab`:

```
1 4 * * * export EAGLEOWL_CONF=/path/to/eagle-owl/config; export DATESTAMP=`date +\%Y\%m\%d`; /path/to/elan/repo/bin/control/go-full-elan.sh $DATESTAMP
```

## Etymology

Elan is named after an aqueduct whose source is in the Elan Valley in Powys, Mid Wales.
In the late 1800s, the Birmingham Corporation Water Department constructed a series of dams in the valley and an aqueduct stretching 73 miles (117 km) to carry water to Birmingham. The aqueduct served a vital supply of clean water to the West Midlands to improve public health in the region.

