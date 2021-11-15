#!/usr/bin/env python
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
parser.add_argument('-d', '--drop', nargs='*', default=[])
parser.add_argument('--host', default="localhost")
args = parser.parse_args()

print("ignoring", args.drop)

def on_connect(client, userdata, flags, rc):
    print("subbed to ", args.topic)
    client.subscribe(args.topic, qos=2)

def on_message_wrap(client, userdata, msg):
    try:
        return on_message(client, userdata, msg)
    except Exception as e:
        print(e)

def on_message(client, userdata, msg):

    try:
        payload = json.loads(msg.payload)
    except Exception as e:
        print("invalid", msg.topic)
        return

    token = os.getenv('SLACK_TOKEN')
    if not token:
        print("boo no token")
        return

    for d in args.drop:
        if d in msg.topic:
            print('ignored (%s)' % d, msg.topic, payload)
            return

    print(msg.topic, payload)

    announce = payload.get("announce", False)

    try:
        sclient = WebClient(token=token)
    except Exception as e:
        print(e)

    try:
        smsg = "\n".join([
                #"*COGUK Interprocess Communication System*",
                "%s`%s`" % ("<!channel>\n" if announce else "", msg.topic),
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
client.on_message = on_message_wrap

client.connect(args.host, 1883, 60)
client.loop_forever()
