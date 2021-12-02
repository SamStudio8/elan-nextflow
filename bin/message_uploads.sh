#!/bin/bash
latest="$1"
q="$2"
t="$3"
msg="$4"
slack="$5"
summary="$6"
DATESTAMP="$7"

echo "[ELAN]" `date` " - Handle messages"
COUNT_MAJORA=`wc -l $latest | cut -f1 -d' '`
COUNT_ELAN_NEW=`grep -c '^1' $q`
COUNT_ELAN_OLD=`awk '{sum+=$6} END {print sum}' $summary`
COUNT_ELAN_OLDANDNEW=`expr $COUNT_ELAN_NEW + $COUNT_ELAN_OLD`
SITE_COUNTS=`awk '$14=="SANG" {print $14 " ("$13")"; next}; {print $14}' $q | sort | uniq -c | sort -nr`
SITE_COUNTS_NEW=`grep '^1' $q | awk '($14=="SANG" || $14=="PHEC") {print $14 " ("$13")"; next}; {print $14}' | sort | uniq -c | sort -nr`

FILE_MISSING_META=`grep 'ORPHAN-NEW-DIRX' $t | grep -v jacksond | awk '$2 > 2 {print $2,$8}' | sort -nr | column -t`
OLD_FILE_MISSING_META=`grep 'ORPHAN-USER-DIRX' $t | grep -v jacksond | cut -f2,3 -d'|' | sed 's,|, ,' | column -t`

SANG_MISSING_META=`grep 'ORPHAN-FILE' $t | grep jacksond | cut -f2 -d' ' | awk -F'/' '{print $(NF-1)}' | cut -c1-4 | tr -d '[0-9]_' | sort | uniq -c | sort -nr`

INVALID_LC=`grep 'BAD-LINES' $t | cut -f2 -d'|'`
if [ $INVALID_LC -gt 0 ]; then
    INVALID_MSG='Additionally, <@U0100HRPMAR>, there are '${INVALID_LC}' bad lines of metadata for you to investigate.'
else
    INVALID_MSG=""
fi

RECENT_DAYS_DEF=`grep 'RECENT-DAYS' $t | cut -f2 -d'|'`
NEW_SITE_MISSING_FILE=`grep 'ORPHAN-NEW-SITE' $t | awk '$2 > 0 {print $6 " " $2}' | sort -k2nr`
OLD_SITE_MISSING_FILE=`grep 'ORPHAN-OLD-SITE' $t | awk '$2 > 0 {print $6 " " $2}' | sort -k2nr`

###############################################################################
if [ "$msg" = "LONG" ]; then

LONG='{
    "attachments": [
        {
            "blocks": [
                {
                    "type": "header",
                    "text": {
                        "type": "plain_text",
                        "text": ":star: COG-UK inbound-distribution pre-pipeline report",
                        "emoji": true
                    }
                },
                {
                    "type": "section",
                    "text": {
                        "type": "mrkdwn",
                        "text": "*'$COUNT_ELAN_NEW'* new sequences today"
                    },
                    "accessory": {
                        "type": "image",
                        "image_url": "https://avatars.slack-edge.com/2019-05-03/627972616934_a621b7d3a28c2b6a7bd1_512.jpg",
                        "alt_text": "Majora is watching."
                    }
                },
                {
                    "type": "divider"
                },
            ],
            "color": "#36C5F0",
        },
        {
            "blocks": [
                {
                    "type": "header",
                    "text": {
                        "type": "plain_text",
                        "text": "Recent samples with metadata but missing uploaded sequences on CLIMB, by sequencing centre",
                        "emoji": true
                    }
                },
                {
                    "type": "section",
                    "text": {
                        "type": "mrkdwn",
                        "text":"
Please check your upload directories...'"\`\`\`${NEW_SITE_MISSING_FILE}\`\`\`"'"
                    }
                },
                {
                    "type": "divider"
                },
                {
                    "type": "context",
                    "elements": [
                        {
                            "type": "mrkdwn",
                            "text":"
To inspect the barcodes with metadata but missing an uploaded sequence from your site, execute:
```grep ORPHAN-COGX '$COG_PUBLISHED_DIR'/elan/'$DATESTAMP'.missing.ls | grep '"'"'\\[BIRM\\]'"'"'```
_Replace BIRM with your site code from the table above. Ensure to maintain the brackets and quotes._
_Recent is defined as sequenced in the past '${RECENT_DAYS_DEF}' days._"
                        }
                    ]
                }
            ],
            "color": "#E01E5A",
        },
        {
            "blocks": [
                {
                    "type": "header",
                    "text": {
                        "type": "plain_text",
                        "text": "Recently uploaded sequences missing metadata, by secondary directory",
                        "emoji": true
                    }
                },
                {
                    "type": "section",
                    "text": {
                        "type": "mrkdwn",
                        "text":"
These directories contain one or more directories with recently uploaded sequences that do not have metadata in Majora.
Please check you have uploaded all your metadata..."
                    }
                }
            ],
            "color": "#E01E5A"
        }
    ]
}'
curl -X POST -H 'Content-type: application/json' --data "${!msg}" "${slack}"


# Message payload too large for attachment block
LONG='{"text": "'"\`\`\`${FILE_MISSING_META}\`\`\`"'"}'
curl -X POST -H 'Content-type: application/json' --data "${!msg}" "${slack}"


LONG='{
    "attachments": [
        {
            "blocks": [
                {
                    "type": "divider"
                },
                {
                    "type": "context",
                    "elements": [
                        {
                            "type": "mrkdwn",
                            "text":"
To inspect the uploaded sequences that are missing metadata from your site, execute:
```grep ORPHAN-FILE '$COG_PUBLISHED_DIR'/elan/'$DATESTAMP'.missing.ls | grep climb-covid19-nichollss```
_Replace climb-covid19-nichollss with the username of the uploader from the table above._
_Recent is defined as uploaded in the past '${RECENT_DAYS_DEF}' days._"
                        }
                    ]
                }
            ],
            "color": "#E01E5A",
        },
        {
            "blocks": [
                {
                    "type": "header",
                    "text": {
                        "type": "plain_text",
                        "text": "All samples sequenced by Sanger missing biosample metadata",
                        "emoji": true
                    }
                },
                {
                    "type": "section",
                    "text": {
                        "type": "mrkdwn",
                        "text":"
Metadata is missing for samples submitted from these sites, for sequencing at Sanger.
Local sites should ensure they have uploaded the biosample-only metadata for samples not sequenced locally.'"\`\`\`${SANG_MISSING_META}\`\`\`"'"
                    }
                },
                {
                    "type": "divider"
                },
                {
                    "type": "context",
                    "elements": [
                        {
                            "type": "mrkdwn",
                            "text":"
To inspect the sequences uploaded from Sanger that are missing metadata from your site, execute:
```grep ORPHAN-FILE '$COG_PUBLISHED_DIR'/elan/'$DATESTAMP'.missing.ls | grep climb-covid19-jacksond | grep BIRM```
_Replace BIRM with your barcode prefix from the table above, without changing the username._
_Sanger sequences are assumed to be uploaded by climb-covid19-jacksond. If that changes, contact Sam._"
                        }
                    ]
                }
            ],
            "color": "#E01E5A",
        }
    ]
}'
curl -X POST -H 'Content-type: application/json' --data "${!msg}" "${slack}"


LONG='{
    "attachments": [
        {
            "blocks": [
                {
                    "type": "header",
                    "text": {
                        "type": "plain_text",
                        "text": "Older uploaded sequences missing metadata by uploading user",
                        "emoji": true
                    }
                },
                {
                    "type": "section",
                    "text": {
                        "type": "mrkdwn",
                        "text":"
These users have orphaned sample sequences that are still missing metadata.
'"\`\`\`${OLD_FILE_MISSING_META}\`\`\`"'",
                    }
                },
                {
                    "type": "divider"
                },
                {
                    "type": "context",
                    "elements": [
                        {
                            "type": "mrkdwn",
                            "text":"
To inspect the uploaded sequences that are missing metadata from your site, execute:
```grep ORPHAN-FILE '$COG_PUBLISHED_DIR'/elan/'$DATESTAMP'.missing.ls | grep climb-covid19-nichollss```
_Replace climb-covid19-nichollss with the username of the uploader from the table above._
_Older is defined as uploaded more than '${RECENT_DAYS_DEF}' days ago._"
                        }
                    ]
                }
            ],
            "color": "#f2c744",
        },
        {
            "blocks": [
                {
                    "type": "header",
                    "text": {
                        "type": "plain_text",
                        "text": "New sequences by centre",
                        "emoji": true
                    }
                },
                {
                    "type": "divider"
                },
                {
                    "type": "section",
                    "text": {
                        "type": "mrkdwn",
                        "text": "_Files and metadata uploaded in the past 30 minutes may not have been seen in time to be counted here. Files and metadata must be uploaded before one minute past the next hour to appear in the end of day report. The pipeline will start autonomously tomorrow morning, not even Sam can stop the pipeline now..._"
                    }
                }
            ],
            "color": "#36a64f"
        }
    ]
}'
curl -X POST -H 'Content-type: application/json' --data "${!msg}" "${slack}"


LONG='{"text": "'"\`\`\`${SITE_COUNTS_NEW}\`\`\`"'"}'
curl -X POST -H 'Content-type: application/json' --data "${!msg}" "${slack}"
fi
###############################################################################
if [ "$msg" = "SHORT" ]; then

SHORT='{
    "attachments": [
        {
            "blocks": [
                {
                    "type": "header",
                    "text": {
                        "type": "plain_text",
                        "text": ":star2: COG-UK inbound-distribution pre-pipeline count",
                        "emoji": true
                    }
                },
                {
                    "type": "section",
                    "text": {
                        "type": "mrkdwn",
                        "text": "
*'$COUNT_MAJORA'* sample sequencing experiments in Majora
*'$COUNT_ELAN_NEW'* new sequences today
*'$COUNT_ELAN_OLDANDNEW'* sequences matched to Majora metadata"
                    },
                    "accessory": {
                        "type": "image",
                        "image_url": "https://avatars.slack-edge.com/2019-05-03/627972616934_a621b7d3a28c2b6a7bd1_512.jpg",
                        "alt_text": "Majora is watching."
                    }
                },
                {
                    "type": "context",
                    "elements": [
                        {
                            "type": "mrkdwn",
                            "text":"
_The pipeline will start tomorrow morning. Good night!_"
                        }
                    ]
                }
            ],
            "color": "#36C5F0",
        }
    ]
}'
curl -X POST -H 'Content-type: application/json' --data "${!msg}" "${slack}"
SHORT='{"text": "*New sequences by centre*\n'"\`\`\`${SITE_COUNTS_NEW}\`\`\`"'"}'
curl -X POST -H 'Content-type: application/json' --data "${!msg}" "${slack}"
fi
###############################################################################
if [ "$msg" = "SHORTSTART" ]; then

SHORTSTART='{
    "attachments": [
        {
            "blocks": [
                {
                    "type": "header",
                    "text": {
                        "type": "plain_text",
                        "text": ":sunrise: COG-UK inbound-distribution pipeline ready",
                        "emoji": true
                    }
                },
                {
                    "type": "section",
                    "text": {
                        "type": "mrkdwn",
                        "text": "
*'$COUNT_MAJORA'* sample sequencing experiments in Majora
*'$COUNT_ELAN_NEW'* new sequences today
*'$COUNT_ELAN_OLDANDNEW'* sequences matched to Majora metadata"
                    },
                    "accessory": {
                        "type": "image",
                        "image_url": "https://avatars.slack-edge.com/2019-05-03/627972616934_a621b7d3a28c2b6a7bd1_512.jpg",
                        "alt_text": "Majora is watching."
                    }
                },
                {
                    "type": "context",
                    "elements": [
                        {
                            "type": "mrkdwn",
                            "text":"
_Good morning! The pipeline will start shortly. Have a nice day!_"
                        }
                    ]
                }
            ],
            "color": "#36C5F0",
        }
    ]
}'
curl -X POST -H 'Content-type: application/json' --data "${!msg}" "${slack}"
SHORTSTART='{"text": "*New sequences by centre*\n'"\`\`\`${SITE_COUNTS_NEW}\`\`\`"'"}'
curl -X POST -H 'Content-type: application/json' --data "${!msg}" "${slack}"
fi
###############################################################################
echo "[ELAN]" `date` " - Pull complete"
