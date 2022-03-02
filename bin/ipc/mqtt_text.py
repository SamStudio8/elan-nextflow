#!/usr/bin/env python3
import paho.mqtt.client as mqtt
import json
import time
import argparse
import os

def on_connect(client, userdata, flags, rc):
    print("subbed to ", topic)
    client.subscribe(topic, qos=2)

def on_message_wrap(client, userdata, msg):
    try:
        return on_message(client, userdata, msg)
    except Exception as e:
        print(e)

def on_message(client, userdata, msg):
    payload = json.loads(msg.payload)
    with open(output, 'a') as out_fh:
        out_fh.write('\t'.join([str(x) for x in [
            msg.topic,
            time.time(),
            payload,
        ]]) + '\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--topic', default="COGUK/#")
    parser.add_argument("--host", default="localhost")
    parser.add_argument("-o", required=True)
    args = parser.parse_args()

    topic = args.topic
    output = args.o
    client = mqtt.Client()
    client.on_connect = on_connect
    client.on_message = on_message_wrap

    client.connect(args.host, 1883, 60)
    client.loop_forever()
