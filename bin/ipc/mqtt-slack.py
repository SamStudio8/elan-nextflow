import paho.mqtt.client as mqtt
import json
import subprocess
import argparse
import os

from slack import WebClient
from slack.errors import SlackApiError

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--topic', default="COGUK/#")
parser.add_argument('-c', '--channel', default="#majora-test")
args = parser.parse_args()

def on_connect(client, userdata, flags, rc):
    print("subbed to ", args.topic)
    client.subscribe(args.topic, qos=2)

def on_message(client, userdata, msg):
    payload = json.loads(msg.payload)
    token = os.getenv('SLACK_TOKEN')
    if not token:
        print("boo no token")
        return

    print(payload)
    try:
        sclient = WebClient(token=token)
    except Exception as e:
        print(e)

    try:
        smsg = "\n".join([
                #"*COGUK Interprocess Communication System*",
                "`%s`" % msg.topic,
                "",
                "```",
                json.dumps(payload, indent=4, sort_keys=True),
                "```",
        ])
        response = sclient.chat_postMessage(
                channel=args.channel,
                text=smsg,
        )
    except SlackApiError as e:
        print("boo bad slack")
    except Exception as e:
        print(e)
        return

client = mqtt.Client()
client.on_connect = on_connect
client.on_message = on_message

client.connect("localhost", 1883, 60)
client.loop_forever()
