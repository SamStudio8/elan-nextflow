import paho.mqtt.client as mqtt
import json
import time
import argparse
import os
import requests
import toml

def on_connect(client, userdata, flags, rc):
    print("subbed to ", topic)
    client.subscribe(topic, qos=2)

def on_message_wrap(client, userdata, msg):
    try:
        return on_message(client, userdata, msg)
    except Exception as e:
        print(e)

def on_message(client, userdata, msg):
    status_code = -77 # EX_NOPERM
    text = {}
    payload = {}

    if msg.topic in allow_set:
        payload = json.loads(msg.payload)
        payload = make_payload(msg.topic, payload, payload["ts"] if payload.get("ts") else None)

        try:
            r = send_splunk(payload, profile)
            status_code = r.status_code
            text = r.text
        except Exception as e:
            status_code = -69 # EX_UNAVAILABLE
            text = {"exception": str(e)}

    with open(output, 'a') as out_fh:
        out_fh.write('\t'.join([str(x) for x in [
            msg.topic,
            time.time(),
            payload,
            status_code,
            json.dumps(text),
        ]]) + '\n')

def get_splunk_creds(profile):
    splunk_toml = os.getenv("SPLUNK_TOML")
    if not splunk_toml:
        raise Exception("Splunk credentials could not be loaded: SPLUNK_TOML unset")
    if not os.path.exists(splunk_toml):
        raise Exception("Splunk credentials could not be loaded: SPLUNK_TOML set but does not exist")

    try:
        config = toml.load(splunk_toml)
    except Exception as e:
        sys.stderr.write("Splunk credentials could not be loaded: SPLUNK_TOML set but encountered an error on open\n")
        raise e

    creds = config.get("profile", {}).get(profile, {})
    for k in ["splunk_auth", "splunk_endpoint"]:
        if k not in creds:
            raise Exception("Splunk credentials could not be loaded: SPUNK_TOML set and valid but required keys appear unset")

    return creds

def make_payload(topic, event, ts=None):
    splunk_payload = {
        "time": ts if ts else int(time.time()),
        "host": "climb-covid",
        "source": topic,
        "event": event,
    }
    return splunk_payload

def send_splunk(payload, profile):

    if len(payload) == 0:
        raise Exception("Payload cannot be empty")

    creds = get_splunk_creds(profile)
    splunk_url = creds["splunk_endpoint"]
    splunk_auth = creds["splunk_auth"]

    headers = {
        "Authorization": splunk_auth,
    }
    return requests.post(splunk_url, json=payload, headers=headers)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--topic', default="COGUK/#")
    parser.add_argument('--allow', help="new-line delimited list of topics to allow")
    parser.add_argument("--host", default="localhost")
    parser.add_argument("-o", required=True)
    parser.add_argument("--profile", required=True)
    args = parser.parse_args()

    output = args.o
    topic = args.topic
    profile = args.profile
    allow_set = set([x.strip() for x in open(args.allow).readlines()])
    client = mqtt.Client()
    client.on_connect = on_connect
    client.on_message = on_message_wrap

    client.connect(args.host, 1883, 60)
    client.loop_forever()
