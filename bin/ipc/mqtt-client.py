import paho.mqtt.client as mqtt
import paho.mqtt.publish as publish
import json
import subprocess
import argparse
import os
from datetime import datetime

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--topic', default="COGUK/infrastructure/pipelines/elan/status")
parser.add_argument('-c', '--cmd', required=True)
parser.add_argument('--who', required=True)
parser.add_argument('--envprefix')
parser.add_argument('--envreq', nargs='+')
args = parser.parse_args()

def emit(who, payload):
    payload["ts"] = int(datetime.now().strftime("%s"))
    publish.single(
        "COGUK/infrastructure/pipelines/%s/status" % who,
        payload=json.dumps(payload),
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

def on_connect(client, userdata, flags, rc):
    print("subbed to ", args.topic)
    client.subscribe(args.topic, qos=2)

def on_message(client, userdata, msg):
    payload = json.loads(msg.payload)
    status = payload.get("status")

    if args.envreq:
        upped_payload = [x.upper() for x in payload.keys()]
        for e in args.envreq:
            if e.upper() not in upped_payload:
                print("cowardly refusing to start command without '%s'. report this to the message sender." % e)
                return

    new_partial_env = {}
    if args.envprefix:
        new_partial_env.update( {"%s_%s" % (args.envprefix, k.upper()): v for k,v in payload.items()} )
    env = os.environ.copy()
    env.update(new_partial_env)

    if status == "finished":
        print("elan finished")
        try:
            print("starting command")
            emit(args.who, {"status": "started", "announce": False})

            print("[cmd] %s" % args.cmd)
            print("[env] %s" % str(new_partial_env))

            proc = subprocess.Popen(
                    args.cmd,
                    shell=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    env=env,
            )
            stdout, stderr = proc.communicate()
            rc = proc.returncode

            emit(args.who, {
                "status": "finished",
                "return_code": rc,
                "announce": True if rc > 0 else False,
            })
            print("finished command with return code %s" % str(rc))
        except Exception as e:
            print("unable to initialise subprocess")
            print(e)

client = mqtt.Client()
client.on_connect = on_connect
client.on_message = on_message

client.connect("localhost", 1883, 60)
client.loop_forever()
