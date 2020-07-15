# COGUK interpipeline comms

## Housekeeping

    conda create -n my-ipc python=python3.7
    conda activate my-ipc
    pip install paho-mqtt

## Send a message

    python mqtt-message.py -t COGUK/infrastructure/pipelines/elan/status --attr status finished

## Receive a message and run an action

    python mqtt-client.py -c '/path/to/my_script.sh'
