import paho.mqtt.client as mqtt
import json
import subprocess
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--topic', default="COGUK/infrastructure/pipelines/elan/status")
parser.add_argument('-c', '--cmd', required=True)
args = parser.parse_args()

def on_connect(client, userdata, flags, rc):
    client.subscribe(args.topic, qos=2)

def on_message(client, userdata, msg):
    payload = json.loads(msg.payload)
    status = payload.get("status")
    if status == "finished":
        print("elan finished")
        try:
            print("starting command")
            print("  %s" % args.cmd)
            subprocess.call(args.cmd, shell=True)
            print("finished command")
        except:
            print("unable to initialise subprocess")

client = mqtt.Client()
client.on_connect = on_connect
client.on_message = on_message

client.connect("localhost", 1883, 60)
client.loop_forever()
