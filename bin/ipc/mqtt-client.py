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
    print("[INFO] listener topic:", args.topic)
    client.subscribe(args.topic, qos=2)
    print("[INFO] control topic:", control_topic)
    client.subscribe(control_topic, qos=2)

def on_message(client, userdata, msg):
    payload = json.loads(msg.payload)
    status = action = None

    print("[RECV] %s" % msg.topic)

    if msg.topic.endswith("status"):
        status = payload.get("status")
        if not status:
            print("[SKIP] received status message without status parameter")
            return
    elif msg.topic.endswith("control"):
        action = payload.get("action")
        if not action:
            print("[SKIP] received control message without action parameter")
            return

    envreq = []
    if args.envreq:
        envreq = [e.upper() for e in args.envreq]
        upped_payload = [x.upper() for x in payload.keys()]
        for e in envreq:
            if e not in upped_payload:
                print("[SKIP] cowardly refusing to start command without required payload key '%s'. report this to the message sender." % e)
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

    if not start_cmd:
        print("[SKIP] cowardly refusing to start command")
        print("       * perhaps status was not 'finished'")
        print("       * perhaps action was not 'raise'")
        return
    else:
        try:
            print("[OKGO] starting command (%s)" % reason)
            extend_payload = {}
            if args.payload_passthrough:
                for p in args.payload_passthrough:
                    p = str(p).upper()
                    if p not in env:
                        print("[WARN] cannot passthrough environment variable '%s'..." % p)
                        print("       * if you are using --envprefix, make sure you have included the prefix to the variable to --payload-passthrough")
                        print("       * if you are not using --envprefix, you need to specifically allow variables for passthrough with --envreq")
                    else:
                        extend_payload[p] = env[p]

            payload = {
                "status": "started",
                "announce": False,
                "reason": reason,
            }
            payload.update(extend_payload)
            emit(args.who, payload)

            print("[OKGO][cmd] %s" % args.cmd)
            print("[OKGO][env] %s" % str(new_partial_env))

            start_time = datetime.now()

            # Ensure all env values are strings
            env = {str(k): str(v) for k,v in env.items()}

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
            print("[DONE] finished command with return code %s" % str(rc))
        except Exception as e:
            print("[FAIL] unable to initialise subprocess:")
            print(e)
            payload = {
                "status": "falsestart",
                "return_code": None,
                "announce": True,
                "time_elapsed": None,
            }
            payload.update(extend_payload)
            emit(args.who, payload)

client = mqtt.Client()
client.on_connect = on_connect
client.on_message = on_message

client.connect(args.host, 1883, 60)
client.loop_forever()
