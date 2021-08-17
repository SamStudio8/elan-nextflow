# elan

Elan is the inbound distribution pipeline for CLIMB-COVID.
Elan is a Nextflow DSL1 pipeline for quality checking dispersed files and publishing them to a common location.

## Parameters and environment variables

### Command line parameters

| Name | Description |
| ---- | ----------- |
| `--datestamp` | YYYYMMDD datestamp to identify today's run |
| `--uploads` | Glob path for CLIMB-COVID user uploads |
| `--publish` | Path to CLIMB-COVID staged artifacts root (nicholsz/) |
| `--cog_publish` | Path to CLIMB-COVID published artifact root (artifacts/) |
| `--minlen` | Minimum genome size required to pass the save_uploads step [int] |
| `--schemegit` | Path to local copy of https://github.com/artic-network/artic-ncov2019 repo |


### Environment variables
#### Elan

| Name | Description |
| ---- | ----------- |
| `ELAN_SLACK_HOOK` | HTTPS hook for posting post-resolve messages to Slack |
| `MAJORA_DOMAIN` | `MAJORA_*` variables control authentication with Majora for Ocarina steps. See https://github.com/SamStudio8/ocarina#readme |
| `MAJORA_USER` | " |
| `MAJORA_TOKEN` | " |
| `MAJORA_CLIENT_ID` | " |
| `MAJORA_CLIENT_SECRET` | " |

Note that Elan will only error if `MAJORA_DOMAIN` is unset, all other `MAJORA_*` variables are not checked.

#### Nextflow

| Name | Description |
| ---- | ----------- |
| `NXF_WORK` | Path to NXF working dir |
| `NXF_DEBUG` | One of: 1,2,3. See https://www.nextflow.io/docs/latest/config.html |


## Invocation

Consult `bin/control/go-full-elan.sh`
