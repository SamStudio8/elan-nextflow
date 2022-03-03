<div align="center">
<p align="center">
   <img src="/images/inbound-ops-patch.png" alt="Inbound Data Operations Badge" width="250"> 
</p>
<h1 align="center">Elan</h1>
</div>

Elan is the inbound distribution pipeline for CLIMB-COVID.
Elan is a Nextflow DSL1 pipeline for quality checking dispersed files and publishing them to a common location.

## Parameters and environment variables

### Controlling Elan

#### go-full-elan.sh

| Name | Description |
| ---- | ----------- |
| `DATESTAMP` | YYYYMMDD datestamp to identify today's run |
| `ELAN_CONFIG` | Path to current Nextflow configuration |
| `ELAN_SOFTWARE_DIR` | Path to local clone of elan-nextflow |
| `ELAN_RUN_DIR` | Path to dir to run Elan from (scratch) |
| `ELAN_DIR` | Path to CLIMB-COVID staged artifacts root (nicholsz/), passed as `--publish` to elan-nf |
| `ARTIFACTS_ROOT` | Path to new CLIMB-COVID published artifact root (/artifacts/), passed as `--artifacts_root` to elan-nf |
| `NEXTFLOW_BIN` | Path to nextflow binary |
| `SLACK_MGMT_HOOK` | Slack HTTPS webhook for posting debug messages |
| `SLACK_REAL_HOOK` | Slack HTTPS webhook for posting inbound-dist messages |
| `MQTT_HOST` | IP for MQTT broker |

`go-full-elan.sh` will immediately terminate with exit 64 (`EX_USAGE`) if any of the listed parameters are missing from the environment.

#### cog-publish.sh

Note these variables are checked inside go-full-elan.sh as it is the main entrypoint.
Additionally, variables defined above may be used in cog-publish.sh without listing them below.

| Name | Description |
| ---- | ----------- |
| `COG_PUBLISH_MODE` | Set to `local` or `slurm` to control how the daily consensus is generated |


### Running Elan

#### elan-nextflow parameters

| Name | Description |
| ---- | ----------- |
| `--datestamp` | YYYYMMDD datestamp to identify today's run |
| `--uploads` | Glob path for CLIMB-COVID user uploads |
| `--publish` | Path to CLIMB-COVID staged artifacts root (nicholsz/) |
| `--artifacts_root` | Path to new CLIMB-COVID published artifact root (/artifacts/) |
| `--minlen` | Minimum genome size required to pass the save_uploads step [int] |


#### elan-nextflow environment variables

| Name | Description |
| ---- | ----------- |
| `ELAN_SLACK_HOOK` | HTTPS hook for posting post-resolve messages to Slack |
| `MAJORA_DOMAIN` | `MAJORA_*` variables control authentication with Majora for Ocarina steps. See https://github.com/SamStudio8/ocarina#readme |
| `MAJORA_USER` | " |
| `MAJORA_TOKEN` | " |
| `MAJORA_CLIENT_ID` | " |
| `MAJORA_CLIENT_SECRET` | " |

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
31 5 * * * DATE=`date +\%Y\%m\%d`; /path/to/elan-nextflow/bin/control/go-full-elan.sh $DATE
```

## Etymology

Elan is named after an aqueduct whose source is in the Elan Valley in Powys, Mid Wales.
In the late 1800s, the Birmingham Corporation Water Department constructed a series of dams in the valley and an aqueduct stretching 73 miles (117 km) to carry water to Birmingham. The aqueduct served a vital supply of clean water to the West Midlands to improve public health in the region.

