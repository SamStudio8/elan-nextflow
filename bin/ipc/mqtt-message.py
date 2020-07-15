import paho.mqtt.client as mqtt
import paho.mqtt.publish as publish
import json

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--topic', required=True)
parser.add_argument("--attr", action='append', nargs=2, metavar=('key', 'value'))
args = parser.parse_args()


payload = json.dumps({x[0]: x[1] for x in args.attr})

publish.single(
    args.topic,
    payload=payload,
    hostname="localhost",
    transport="tcp",
    port=1883,
    qos=2,
    client_id="",
    keepalive=60,
    retain=False,
    will=None,
    auth=None,
    tls=None,
    protocol=mqtt.MQTTv311,
)
