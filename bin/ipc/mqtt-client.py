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
parser.add_argument('--payload-passthrough', nargs='+')
parser.add_argument("--host", default="localhost")
args = parser.parse_args()

def emit(who, payload):
    payload["ts"] = int(datetime.now().strftime("%s"))
    publish.single(
        "COGUK/infrastructure/pipelines/%s/status" % who,
        payload=json.dumps(payload),
        hostname=args.host,
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
    control_topic = "COGUK/infrastructure/pipelines/%s/control" % args.who
    print("subbed to", args.topic)
    client.subscribe(args.topic, qos=2)
    print("control topic", control_topic)
    client.subscribe(control_topic, qos=2)

def on_message(client, userdata, msg):
    payload = json.loads(msg.payload)
    status = action = None
    if msg.topic.endswith("status"):
        status = payload.get("status")
        if not status:
            print("received status message without status parameter")
            return
    elif msg.topic.endswith("control"):
        action = payload.get("action")
        if not action:
            print("received control message without action parameter")
            return

    envreq = []
    if args.envreq:
        envreq = [e.upper() for e in args.envreq]
        upped_payload = [x.upper() for x in payload.keys()]
        for e in envreq:
            if e not in upped_payload:
                print("cowardly refusing to start command without '%s'. report this to the message sender." % e)
                return

    new_partial_env = {}
    if args.envprefix:
        # Use the prefix to safely add all payload variables to the env as a form of namespacing
        new_partial_env.update( {"%s_%s" % (args.envprefix, k.upper()): v for k,v in payload.items()} )
    else:
        # Without envprefix, only push through payload variables that are required by envreq
        new_partial_env.update({"%s" % k.upper(): v for k,v in payload.items() if k.upper() in envreq})
    env = os.environ.copy()
    env.update(new_partial_env)

    start_cmd = False
    reason = "unknown"
    if status == "finished":
        start_cmd = True
        reason = "%s finished" % args.topic

    if action == "raise":
        start_cmd = True
        reason = "manually raised by control message"

    if start_cmd:

        print(reason)
        try:
            print("starting command")
            extend_payload = {}
            if args.payload_passthrough:
                for p in args.payload_passthrough:
                    if p not in env:
                        print("cannot passthrough environment variable '%s'. make sure to use uppercasing. if you are using the envreq prefix, make sure to use the full prefixed name." % p)
                    else:
                        extend_payload[p] = env[p]

            payload = {
                "status": "started",
                "announce": False,
                "reason": reason,
            }
            payload.update(extend_payload)
            emit(args.who, payload)

            print("[cmd] %s" % args.cmd)
            print("[env] %s" % str(new_partial_env))

            start_time = datetime.now()

            proc = subprocess.Popen(
                    args.cmd,
                    shell=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    env=env,
            )
            stdout, stderr = proc.communicate()
            rc = proc.returncode

            end_time = datetime.now()

            if rc == 0:
                status = "finished"
            else:
                status = "failed"

            payload = {
                "status": status,
                "return_code": rc,
                "announce": True if rc > 0 else False,
                "time_elapsed": str(end_time - start_time),
            }
            payload.update(extend_payload)
            emit(args.who, payload)
            print("finished command with return code %s" % str(rc))
        except Exception as e:
            print("unable to initialise subprocess")
            print(e)

client = mqtt.Client()
client.on_connect = on_connect
client.on_message = on_message

client.connect(args.host, 1883, 60)
client.loop_forever()
